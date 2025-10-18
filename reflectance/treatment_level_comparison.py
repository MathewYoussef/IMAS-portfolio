#!/usr/bin/env python3
"""Compare dose-level reflectance summaries with DAD concentrations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats

from dose_metadata import attach_dose_metadata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate dose-level features and compare reflectance vs DAD statistics."
    )
    parser.add_argument(
        "--reflectance-stats",
        type=Path,
        default=Path("canonical_dataset/dose_reflectance_stats.csv"),
        help="Dose-level per-angle reflectance summary (from canonical_dataset).",
    )
    parser.add_argument(
        "--reflectance-composites",
        type=Path,
        default=Path("canonical_dataset/dose_reflectance_composites.csv"),
        help="Dose-level composite spectra (Σ/Δ) in canonical_dataset.",
    )
    parser.add_argument(
        "--dad-stats",
        type=Path,
        default=Path("canonical_dataset/dose_dad_concentrations.csv"),
        help="Dose-level DAD concentration summaries (canonical).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("canonical_dataset"),
        help="Directory to write comparison outputs (will co-locate with canonical tables).",
    )
    parser.add_argument(
        "--wavelength-start",
        type=float,
        default=350.0,
        help="Wavelength corresponding to index 0 of the spectra (nm).",
    )
    parser.add_argument(
        "--wavelength-step",
        type=float,
        default=1.0,
        help="Increment between wavelength samples (nm).",
    )
    parser.add_argument(
        "--subset-min",
        type=float,
        default=320.0,
        help="Lower bound (inclusive) of wavelength subset used for feature extraction (nm).",
    )
    parser.add_argument(
        "--subset-max",
        type=float,
        default=480.0,
        help="Upper bound (inclusive) of wavelength subset used for feature extraction (nm).",
    )
    return parser.parse_args()


def extract_matrix(row: pd.Series, prefix: str) -> np.ndarray:
    cols = sorted([col for col in row.index if col.startswith(prefix)])
    values = row[cols].to_numpy(dtype=float)
    return values


def summarise_spectrum(mean_vals: np.ndarray) -> Dict[str, float]:
    return {
        "mean_reflectance": float(mean_vals.mean()),
        "median_reflectance": float(np.median(mean_vals)),
        "min_reflectance": float(mean_vals.min()),
        "max_reflectance": float(mean_vals.max()),
        "area_reflectance": float(mean_vals.sum()),
        "std_reflectance": float(mean_vals.std(ddof=1)),
        "range_reflectance": float(mean_vals.max() - mean_vals.min()),
    }


def summarise_dispersion(std_vals: np.ndarray, mad_vals: np.ndarray) -> Dict[str, float]:
    return {
        "mean_std": float(std_vals.mean()),
        "median_std": float(np.median(std_vals)),
        "mean_mad": float(mad_vals.mean()),
        "median_mad": float(np.median(mad_vals)),
    }


def build_feature_table(
    df: pd.DataFrame,
    label_column: str,
    mask: np.ndarray,
) -> pd.DataFrame:
    records: List[Dict[str, object]] = []
    for _, row in df.iterrows():
        mean_vals = extract_matrix(row, "mean_")
        std_vals = extract_matrix(row, "std_")
        mad_vals = extract_matrix(row, "mad_")
        mean_vals = mean_vals[mask]
        std_vals = std_vals[mask]
        mad_vals = mad_vals[mask]
        summary = summarise_spectrum(mean_vals)
        dispersion = summarise_dispersion(std_vals, mad_vals)
        dose_meta = attach_dose_metadata(row["dose_id"])
        record = {
            "dose_id": row["dose_id"],
            "kind": row[label_column],
            "uva_mw_cm2": float(row.get("uva_mw_cm2", dose_meta.uva_mw_cm2)),
            "uvb_mw_cm2": float(row.get("uvb_mw_cm2", dose_meta.uvb_mw_cm2)),
        }
        record.update(summary)
        record.update(dispersion)
        record["samples"] = row.get("samples_used", row.get("samples_12", None))
        record["trim_fraction"] = row.get("trim_fraction")
        records.append(record)
    return pd.DataFrame(records)


def compute_correlations(ref_df: pd.DataFrame, dad_df: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    """Return Pearson/Spearman correlations between reflectance metrics and DAD metrics."""
    merged = ref_df.merge(
        dad_df,
        on=["uva_mw_cm2", "uvb_mw_cm2"],
        how="inner",
        suffixes=("_reflectance", "_dad"),
    )
    if merged.empty:
        raise ValueError("No overlapping dose records between reflectance and DAD tables.")
    results: Dict[str, Dict[str, float]] = {}
    for metric in ["dad_total_mg_per_gDW_mean_trimmed", "dad_oxidized_mg_per_gDW_mean_trimmed", "dad_reduced_mg_per_gDW_mean_trimmed"]:
        for refl_metric in [
            "mean_reflectance",
            "median_reflectance",
            "area_reflectance",
            "min_reflectance",
            "max_reflectance",
            "range_reflectance",
        ]:
            x = merged[refl_metric].to_numpy()
            y = merged[metric].to_numpy()
            pearson_r, pearson_p = stats.pearsonr(x, y)
            spearman_rho, spearman_p = stats.spearmanr(x, y)
            key = f"{refl_metric} vs {metric}"
            results[key] = {
                "pearson_r": float(pearson_r),
                "pearson_p": float(pearson_p),
                "spearman_rho": float(spearman_rho),
                "spearman_p": float(spearman_p),
            }
    return results


def main() -> None:
    args = parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    angle_df = pd.read_csv(args.reflectance_stats)
    composite_df = pd.read_csv(args.reflectance_composites)
    dad_df = pd.read_csv(args.dad_stats)

    for frame in (angle_df, composite_df, dad_df):
        if "dose_id" in frame.columns:
            if "uva_mw_cm2" not in frame.columns or frame["uva_mw_cm2"].isnull().any():
                frame["uva_mw_cm2"] = frame["dose_id"].map(
                    lambda label: attach_dose_metadata(label).uva_mw_cm2
                )
            if "uvb_mw_cm2" not in frame.columns or frame["uvb_mw_cm2"].isnull().any():
                frame["uvb_mw_cm2"] = frame["dose_id"].map(
                    lambda label: attach_dose_metadata(label).uvb_mw_cm2
                )

    n_wavelengths = len([col for col in angle_df.columns if col.startswith("mean_")])
    wavelengths = args.wavelength_start + args.wavelength_step * np.arange(n_wavelengths)
    mask = (wavelengths >= args.subset_min) & (wavelengths <= args.subset_max)
    if not mask.any():
        raise ValueError(
            "Wavelength subset resulted in empty selection. "
            "Check start/step/subset arguments."
        )

    # Build feature tables
    angle_features = build_feature_table(angle_df, "angle", mask)
    composite_features = build_feature_table(composite_df, "composite", mask)

    angle_features.to_csv(outdir / "dose_reflectance_features_angles.csv", index=False)
    composite_features.to_csv(outdir / "dose_reflectance_features_composites.csv", index=False)

    # Join Sigma composite with DAD stats for direct comparison
    sigma_features = composite_features[composite_features["kind"] == "Sigma"].copy()
    comparison_df = sigma_features.merge(
        dad_df,
        on=["uva_mw_cm2", "uvb_mw_cm2"],
        how="inner",
        suffixes=("_reflectance", "_dad"),
    )
    comparison_df.to_csv(outdir / "dose_reflectance_dad_summary.csv", index=False)

    # Compute correlations
    corr_results = compute_correlations(sigma_features, dad_df)
    corr_path = outdir / "dose_reflectance_dad_correlations.json"
    corr_path.write_text(json.dumps(corr_results, indent=2))

    print("Angle features saved to", angle_features.shape)
    print("Composite features saved to", composite_features.shape)
    print("Joined Sigma + DAD summary saved with shape", comparison_df.shape)
    print("Correlations written to", corr_path)


if __name__ == "__main__":
    main()
