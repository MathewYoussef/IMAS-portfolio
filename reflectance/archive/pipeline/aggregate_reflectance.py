#!/usr/bin/env python3
"""Aggregate denoised reflectance replicates into robust mean spectra.

The script expects the reflectance directory structure:
    root/
        treatment_*/
            sample_*/
                12Oclock/
                    rep_###_denoised.npy
                6Oclock/
                    rep_###_denoised.npy

For every directory (historically named `treatment_*`, but mapped to dose IDs),
sample, and angle combination, the script loads all replicate spectra, computes
a trimmed mean and trimmed standard deviation (default: 10 % trim from each tail),
and writes:
    - a CSV summary table with dose/sample/angle metadata, mean/std columns per wavelength
    - NPZ bundles containing the mean and standard deviation spectra

Usage example:
    python aggregate_reflectance.py \\
        --root denoised_full_run_staging \\
        --outdir aggregated_reflectance \\
        --trim-fraction 0.1
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

from dose_metadata import attach_dose_metadata

def iter_angle_files(angle_dir: Path) -> Iterable[Path]:
    """Yield replicate files in sorted order for a given angle directory."""
    return sorted(angle_dir.glob("rep_*_denoised.npy"))


def compute_trimmed_stats(
    data: np.ndarray, trim_fraction: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Return trimmed mean, trimmed std (ddof=1), MAD, and replicates retained."""
    if data.ndim != 2:
        raise ValueError(f"Expected 2D array, got shape {data.shape}")

    n_reps = data.shape[0]
    if not 0 <= trim_fraction < 0.5:
        raise ValueError("trim_fraction must be in [0, 0.5)")

    if trim_fraction == 0 or n_reps <= 2:
        trimmed = data
    else:
        cut = int(np.floor(trim_fraction * n_reps))
        if cut == 0:
            trimmed = data
        else:
            if cut * 2 >= n_reps:
                cut = max(0, (n_reps - 1) // 2)
            sorted_data = np.sort(data, axis=0)
            trimmed = sorted_data[cut : n_reps - cut]
            if trimmed.size == 0:
                trimmed = sorted_data

    mean = trimmed.mean(axis=0)
    if trimmed.shape[0] > 1:
        std = trimmed.std(axis=0, ddof=1)
        median = np.median(trimmed, axis=0)
        mad = np.median(np.abs(trimmed - median), axis=0)
    else:
        std = np.zeros_like(mean)
        mad = np.zeros_like(mean)
    return mean, std, mad, int(trimmed.shape[0])


def sanitize_key(*parts: str) -> str:
    """Generate a safe key for NPZ storage."""
    return "_".join(part.replace(" ", "").replace("/", "_") for part in parts)


def treatment_dir_to_dose(treatment_dir: Path) -> str:
    """Convert the on-disk treatment directory name into a canonical dose id."""
    name = treatment_dir.name
    if name.startswith("treatment_"):
        try:
            ordinal = int(name.split("_", 1)[1])
        except ValueError as exc:
            raise ValueError(f"Unparseable treatment directory name: {name}") from exc
        if not 1 <= ordinal <= 6:
            raise ValueError(f"Unexpected treatment ordinal {ordinal} in {name}")
        # The historical treatment labels run opposite the true dose order.
        dose_idx = 7 - ordinal  # treatment_1 -> dose_6 (highest UVA), treatment_6 -> dose_1.
        return f"dose_{dose_idx}"
    return name


def aggregate_reflectance(
    root: Path, trim_fraction: float
) -> Tuple[pd.DataFrame, Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Aggregate all spectra under root and return summary objects."""
    records: List[Dict[str, object]] = []
    mean_store: Dict[str, np.ndarray] = {}
    std_store: Dict[str, np.ndarray] = {}
    mad_store: Dict[str, np.ndarray] = {}
    wavelength_count: int | None = None

    for source_dir in sorted(root.iterdir()):
        if not source_dir.is_dir():
            continue
        dose_id = treatment_dir_to_dose(source_dir)
        for sample_dir in sorted(source_dir.iterdir()):
            if not sample_dir.is_dir():
                continue
            for angle_dir in sorted(sample_dir.iterdir()):
                if not angle_dir.is_dir():
                    continue
                files = list(iter_angle_files(angle_dir))
                if not files:
                    continue
                spectra = np.stack([np.load(path) for path in files], axis=0)
                if wavelength_count is None:
                    wavelength_count = spectra.shape[1]
                elif spectra.shape[1] != wavelength_count:
                    raise ValueError(
                        f"Inconsistent spectrum length in {angle_dir}: "
                        f"{spectra.shape[1]} vs expected {wavelength_count}"
                    )

                mean, std, mad, kept = compute_trimmed_stats(spectra, trim_fraction)
                key = sanitize_key(dose_id, sample_dir.name, angle_dir.name, "mean")
                mean_store[key] = mean
                std_store[sanitize_key(dose_id, sample_dir.name, angle_dir.name, "std")] = std
                mad_store[sanitize_key(dose_id, sample_dir.name, angle_dir.name, "mad")] = mad

                dose_meta = attach_dose_metadata(dose_id)

                record: Dict[str, object] = {
                    "dose_id": dose_id,
                    "source_dir": source_dir.name,
                    "sample_label": sample_dir.name,
                    "angle": angle_dir.name,
                    "replicates_total": spectra.shape[0],
                    "replicates_used": kept,
                    "trim_fraction": trim_fraction,
                    "uva_mw_cm2": dose_meta.uva_mw_cm2,
                    "uvb_mw_cm2": dose_meta.uvb_mw_cm2,
                }
                for idx in range(mean.shape[0]):
                    record[f"mean_{idx:03d}"] = float(mean[idx])
                    record[f"std_{idx:03d}"] = float(std[idx])
                    record[f"mad_{idx:03d}"] = float(mad[idx])
                records.append(record)

    if not records:
        raise RuntimeError(f"No spectra found under {root}")

    df = pd.DataFrame.from_records(records)
    return df, mean_store, std_store, mad_store


def save_outputs(
    df: pd.DataFrame,
    mean_store: Dict[str, np.ndarray],
    std_store: Dict[str, np.ndarray],
    mad_store: Dict[str, np.ndarray],
    outdir: Path,
) -> None:
    """Persist aggregated results to disk."""
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "reflectance_trimmed_stats.csv"
    df.to_csv(csv_path, index=False)

    mean_npz_path = outdir / "reflectance_trimmed_means.npz"
    np.savez(mean_npz_path, **mean_store)

    std_npz_path = outdir / "reflectance_trimmed_stds.npz"
    np.savez(std_npz_path, **std_store)

    mad_npz_path = outdir / "reflectance_trimmed_mads.npz"
    np.savez(mad_npz_path, **mad_store)

    meta_path = outdir / "aggregation_metadata.json"
    meta = {
        "records": int(df.shape[0]),
        "wavelength_count": len([c for c in df.columns if c.startswith("mean_")]),
        "trim_fraction": float(df["trim_fraction"].iloc[0]),
        "dose_ids": sorted(df["dose_id"].unique()),
        "uva_uvb_by_dose": {
            dose_id: {
                "uva_mw_cm2": float(group["uva_mw_cm2"].iloc[0]),
                "uvb_mw_cm2": float(group["uvb_mw_cm2"].iloc[0]),
            }
            for dose_id, group in df.groupby("dose_id")
        },
    }
    meta_path.write_text(json.dumps(meta, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate denoised reflectance replicates using trimmed means."
    )
    parser.add_argument(
        "--root",
        type=Path,
        required=True,
        help="Path to the reflectance directory (e.g., denoised_full_run_staging).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("aggregated_reflectance"),
        help="Directory for aggregated outputs.",
    )
    parser.add_argument(
        "--trim-fraction",
        type=float,
        default=0.1,
        help="Fraction (0–0.5) to trim from each tail when computing means.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df, mean_store, std_store, mad_store = aggregate_reflectance(args.root, args.trim_fraction)
    save_outputs(df, mean_store, std_store, mad_store, args.outdir)
    print(f"Aggregated {len(df)} records. Output saved to {args.outdir}")


if __name__ == "__main__":
    main()
