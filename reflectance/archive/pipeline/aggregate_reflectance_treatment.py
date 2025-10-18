#!/usr/bin/env python3
"""Aggregate sample-level reflectance spectra to the dose level.

Inputs:
    - aggregated_reflectance/reflectance_trimmed_stats.csv (sample-level means/stds)

Outputs (written to aggregated_reflectance/):
    - dose_reflectance_stats.csv: dose-level means and dispersion
    - dose_reflectance_composites.csv: Σ (angle average) and Δ (angle difference)
    - dose_reflectance_metadata.json: summary of parameters

Aggregation approach (applied per dose and angle):
    1. Stack the five sample-level mean spectra, compute a trimmed mean (default: 0.2),
       and capture the trimmed standard deviation alongside the median absolute deviation (MAD).
    2. Build Σ = (mean_12 + mean_6) / 2 and Δ = mean_12 - mean_6 spectra with
       dispersion approximated via simple propagation rules.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from dose_metadata import attach_dose_metadata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate sample-level reflectance spectra to dose-level summaries."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("aggregated_reflectance/reflectance_trimmed_stats.csv"),
        help="Sample-level reflectance stats produced by aggregate_reflectance.py",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("aggregated_reflectance"),
        help="Output directory for dose summaries",
    )
    parser.add_argument(
        "--trim-fraction",
        type=float,
        default=0.2,
        help="Symmetric trim fraction applied when aggregating sample means",
    )
    return parser.parse_args()


def trimmed_mean(data: np.ndarray, trim_fraction: float) -> Tuple[np.ndarray, np.ndarray, int]:
    """Compute trimmed mean and trimmed std across the first axis."""
    if data.ndim != 2:
        raise ValueError(f"Expected 2D array (samples x wavelengths), got {data.shape}")
    n_samples = data.shape[0]
    if not 0 <= trim_fraction < 0.5:
        raise ValueError("trim_fraction must be in [0, 0.5)")
    if n_samples == 0:
        raise ValueError("No spectra provided for aggregation.")

    if trim_fraction == 0 or n_samples <= 2:
        trimmed = data
    else:
        cut = int(np.floor(trim_fraction * n_samples))
        if cut * 2 >= n_samples:
            cut = max(0, (n_samples - 1) // 2)
        sorted_data = np.sort(data, axis=0)
        trimmed = sorted_data[cut : n_samples - cut]
        if trimmed.size == 0:
            trimmed = sorted_data

    mean = trimmed.mean(axis=0)
    if trimmed.shape[0] > 1:
        std = trimmed.std(axis=0, ddof=1)
    else:
        std = np.zeros_like(mean)
    return mean, std, int(trimmed.shape[0])


def mad(data: np.ndarray) -> np.ndarray:
    """Median absolute deviation along axis 0."""
    median = np.median(data, axis=0)
    return np.median(np.abs(data - median), axis=0)


def aggregate_doses(
    df: pd.DataFrame, trim_fraction: float
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate sample-level means per dose and angle."""
    wavelength_cols = [c for c in df.columns if c.startswith("mean_")]
    doses: List[Dict[str, object]] = []
    composites = []

    for (dose_id, angle), group in df.groupby(["dose_id", "angle"]):
        samples = group["sample_label"].tolist()
        spectra = group[wavelength_cols].to_numpy()
        mean, std, kept = trimmed_mean(spectra, trim_fraction)
        mad_vals = mad(spectra)
        dose_meta = attach_dose_metadata(dose_id)

        record = {
            "dose_id": dose_id,
            "angle": angle,
            "samples_total": len(samples),
            "samples_used": kept,
            "trim_fraction": trim_fraction,
            "uva_mw_cm2": dose_meta.uva_mw_cm2,
            "uvb_mw_cm2": dose_meta.uvb_mw_cm2,
        }
        for idx, col in enumerate(wavelength_cols):
            record[f"mean_{idx:03d}"] = float(mean[idx])
            record[f"std_{idx:03d}"] = float(std[idx])
            record[f"mad_{idx:03d}"] = float(mad_vals[idx])
        doses.append(record)

    dose_df = pd.DataFrame(doses)

    # Build composites (Σ and Δ) per dose
    for dose_id, group in dose_df.groupby("dose_id"):
        if group.shape[0] != 2 or set(group["angle"]) != {"12Oclock", "6Oclock"}:
            raise ValueError(f"Expected both angles for {dose_id}, found {group['angle'].tolist()}")

        group = group.set_index("angle")
        dose_meta = attach_dose_metadata(dose_id)
        mean_12 = group.loc["12Oclock", [c for c in dose_df.columns if c.startswith("mean_")]].to_numpy()
        mean_6 = group.loc["6Oclock", [c for c in dose_df.columns if c.startswith("mean_")]].to_numpy()

        std_12 = group.loc["12Oclock", [c for c in dose_df.columns if c.startswith("std_")]].to_numpy()
        std_6 = group.loc["6Oclock", [c for c in dose_df.columns if c.startswith("std_")]].to_numpy()

        mad_12 = group.loc["12Oclock", [c for c in dose_df.columns if c.startswith("mad_")]].to_numpy()
        mad_6 = group.loc["6Oclock", [c for c in dose_df.columns if c.startswith("mad_")]].to_numpy()

        # Convert to arrays of floats
        mean_12 = mean_12.astype(float)
        mean_6 = mean_6.astype(float)
        std_12 = std_12.astype(float)
        std_6 = std_6.astype(float)
        mad_12 = mad_12.astype(float)
        mad_6 = mad_6.astype(float)

        sigma = 0.5 * (mean_12 + mean_6)
        delta = mean_12 - mean_6

        # Rough propagation for dispersion
        sigma_std = 0.5 * np.sqrt(std_12**2 + std_6**2)
        delta_std = np.sqrt(std_12**2 + std_6**2)
        sigma_mad = 0.5 * np.sqrt(mad_12**2 + mad_6**2)
        delta_mad = np.sqrt(mad_12**2 + mad_6**2)

        base = {
            "dose_id": dose_id,
            "samples_12": int(group.loc["12Oclock", "samples_used"]),
            "samples_6": int(group.loc["6Oclock", "samples_used"]),
            "trim_fraction": trim_fraction,
            "uva_mw_cm2": dose_meta.uva_mw_cm2,
            "uvb_mw_cm2": dose_meta.uvb_mw_cm2,
        }
        sigma_record = base.copy()
        sigma_record["composite"] = "Sigma"
        delta_record = base.copy()
        delta_record["composite"] = "Delta"

        for idx, col in enumerate([c for c in dose_df.columns if c.startswith("mean_")]):
            sigma_record[f"mean_{idx:03d}"] = float(sigma[idx])
            sigma_record[f"std_{idx:03d}"] = float(sigma_std[idx])
            sigma_record[f"mad_{idx:03d}"] = float(sigma_mad[idx])

            delta_record[f"mean_{idx:03d}"] = float(delta[idx])
            delta_record[f"std_{idx:03d}"] = float(delta_std[idx])
            delta_record[f"mad_{idx:03d}"] = float(delta_mad[idx])

        composites.append(sigma_record)
        composites.append(delta_record)

    composite_df = pd.DataFrame(composites)
    return dose_df, composite_df


def save_outputs(
    dose_df: pd.DataFrame,
    composite_df: pd.DataFrame,
    outdir: Path,
    trim_fraction: float,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    treatment_path = outdir / "dose_reflectance_stats.csv"
    composite_path = outdir / "dose_reflectance_composites.csv"
    meta_path = outdir / "dose_reflectance_metadata.json"

    dose_df.to_csv(treatment_path, index=False)
    composite_df.to_csv(composite_path, index=False)

    meta = {
        "dose_ids": sorted(dose_df["dose_id"].unique()),
        "angles": sorted(dose_df["angle"].unique()),
        "wavelength_count": len([c for c in dose_df.columns if c.startswith("mean_")]),
        "trim_fraction": trim_fraction,
        "uva_uvb_by_dose": {
            dose_id: {
                "uva_mw_cm2": float(group["uva_mw_cm2"].iloc[0]),
                "uvb_mw_cm2": float(group["uvb_mw_cm2"].iloc[0]),
            }
            for dose_id, group in dose_df.groupby("dose_id")
        },
    }
    meta_path.write_text(json.dumps(meta, indent=2))


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    dose_df, composite_df = aggregate_doses(df, args.trim_fraction)
    save_outputs(dose_df, composite_df, args.outdir, args.trim_fraction)
    print(
        f"Dose summaries written to {args.outdir} "
        f"(stats: {len(dose_df)} rows; composites: {len(composite_df)} rows)."
    )


if __name__ == "__main__":
    main()
