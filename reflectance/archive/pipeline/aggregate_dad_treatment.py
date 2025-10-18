#!/usr/bin/env python3
"""Aggregate DAD concentration metrics by dose with robust statistics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dose_metadata import attach_dose_metadata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute dose-level summaries for DAD concentrations."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("Combined_Scytonemin_Concentrations.csv"),
        help="CSV containing sample-level DAD concentrations (with sample_id).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("aggregated_reflectance"),
        help="Directory where dose-level summaries will be written.",
    )
    parser.add_argument(
        "--trim-fraction",
        type=float,
        default=0.2,
        help="Symmetric trim fraction for robust means/SDs per dose.",
    )
    return parser.parse_args()


def trimmed_stats(values: np.ndarray, trim_fraction: float) -> tuple[float, float, int]:
    """Return trimmed mean, trimmed std, and retained count."""
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return np.nan, np.nan, 0
    if trim_fraction < 0 or trim_fraction >= 0.5:
        raise ValueError("trim_fraction must be in [0, 0.5)")

    n = values.size
    if trim_fraction == 0 or n <= 2:
        trimmed = np.sort(values)
        kept = n
    else:
        cut = int(np.floor(trim_fraction * n))
        if cut * 2 >= n:
            cut = max(0, (n - 1) // 2)
        sorted_vals = np.sort(values)
        trimmed = sorted_vals[cut : n - cut]
        kept = trimmed.size
        if kept == 0:
            trimmed = sorted_vals
            kept = trimmed.size
    mean = float(trimmed.mean())
    std = float(trimmed.std(ddof=1)) if kept > 1 else float("nan")
    return mean, std, kept


def mad(values: np.ndarray) -> float:
    """Median absolute deviation."""
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return float("nan")
    median = np.median(values)
    return float(np.median(np.abs(values - median)))


def aggregate_dose_dad(
    df: pd.DataFrame, trim_fraction: float
) -> pd.DataFrame:
    """Compute robust dose-level stats for DAD concentrations."""
    metrics = [
        "dad_total_mg_per_gDW",
        "dad_oxidized_mg_per_gDW",
        "dad_reduced_mg_per_gDW",
    ]
    doses = []

    for dose_id, group in df.groupby("dose_id"):
        dose_meta = attach_dose_metadata(dose_id)
        record = {
            "dose_id": dose_id,
            "sample_count": len(group),
            "trim_fraction": trim_fraction,
            "source_sample_ids": ",".join(sorted(group["sample_id"].tolist())),
            "uva_mw_cm2": dose_meta.uva_mw_cm2,
            "uvb_mw_cm2": dose_meta.uvb_mw_cm2,
        }
        for metric in metrics:
            values = group[metric].dropna().to_numpy()
            mean_trimmed, std_trimmed, kept = trimmed_stats(values, trim_fraction)
            median_val = float(np.median(values)) if values.size else float("nan")
            mad_val = mad(values)
            record[f"{metric}_mean_trimmed"] = mean_trimmed
            record[f"{metric}_std_trimmed"] = std_trimmed
            record[f"{metric}_median"] = median_val
            record[f"{metric}_mad"] = mad_val
            record[f"{metric}_n_used"] = kept
        doses.append(record)

    dose_df = pd.DataFrame(doses).sort_values("dose_id")
    return dose_df


def save_outputs(df: pd.DataFrame, outdir: Path, trim_fraction: float) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / "dose_dad_concentrations.csv"
    meta_path = outdir / "dose_dad_metadata.json"

    df.to_csv(output_path, index=False)
    meta = {
        "dose_ids": sorted(df["dose_id"].unique()),
        "trim_fraction": trim_fraction,
        "metrics": ["dad_total_mg_per_gDW", "dad_oxidized_mg_per_gDW", "dad_reduced_mg_per_gDW"],
        "uva_uvb_by_dose": {
            dose_id: {
                "uva_mw_cm2": float(group["uva_mw_cm2"].iloc[0]),
                "uvb_mw_cm2": float(group["uvb_mw_cm2"].iloc[0]),
            }
            for dose_id, group in df.groupby("dose_id")
        },
    }
    meta_path.write_text(json.dumps(meta, indent=2))


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    if "sample_id" not in df.columns:
        raise ValueError("Input CSV must contain 'sample_id'.")
    treatment_labels = df["sample_id"].str[0].astype(int)
    df["dose_id"] = treatment_labels.apply(lambda x: f"dose_{7 - x}")

    dose_df = aggregate_dose_dad(df, args.trim_fraction)
    save_outputs(dose_df, args.outdir, args.trim_fraction)
    print(
        f"DAD dose summaries written to {args.outdir} "
        f"({dose_df.shape[0]} doses)."
    )


if __name__ == "__main__":
    main()
