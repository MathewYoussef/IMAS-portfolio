#!/usr/bin/env python3
"""
Derive Phase 3 reflectance features over specified wavelength windows.

The script reads the canonical dose-level reflectance tables via the Phase 3
loader manifest, computes additional summary metrics for configurable
wavelength bands, and writes the results to analysis/phase3/.

Outputs:
    - CSV files for per-angle and composite (Σ/Δ) metrics per window.
    - A manifest JSON describing window bounds and feature columns.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


def load_manifest(manifest_path: Path) -> Dict[str, object]:
    data = json.loads(manifest_path.read_text())
    return data


def extract_wavelength_axis(manifest: Dict[str, object]) -> np.ndarray:
    axis_info = manifest.get("wavelength_axis", {})
    start = axis_info.get("start_nm")
    step = axis_info.get("step_nm")
    count = axis_info.get("count")
    if start is None or step is None or count is None:
        raise ValueError("Manifest missing wavelength axis metadata; rerun loader with wavelength grid.")
    return start + step * np.arange(count, dtype=float)


def load_reflectance_tables(canonical_root: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    angles = pd.read_csv(canonical_root / "dose_reflectance_stats.csv")
    composites = pd.read_csv(canonical_root / "dose_reflectance_composites.csv")
    return angles, composites


def build_window_mask(wavelengths: np.ndarray, lower: float, upper: float) -> np.ndarray:
    mask = (wavelengths >= lower) & (wavelengths <= upper)
    if not mask.any():
        raise ValueError(f"Wavelength window {lower}-{upper} nm selects no points.")
    return mask


def compute_features(
    df: pd.DataFrame,
    wavelengths: np.ndarray,
    mask: np.ndarray,
    group_columns: List[str],
) -> pd.DataFrame:
    spectral_cols = [c for c in df.columns if c.startswith("mean_")]
    mean_matrix = df[spectral_cols].to_numpy(dtype=float)
    std_cols = [c for c in df.columns if c.startswith("std_")]
    mad_cols = [c for c in df.columns if c.startswith("mad_")]

    masked_means = mean_matrix[:, mask]
    features = {
        "mean_reflectance": masked_means.mean(axis=1),
        "median_reflectance": np.median(masked_means, axis=1),
        "min_reflectance": masked_means.min(axis=1),
        "max_reflectance": masked_means.max(axis=1),
        "range_reflectance": masked_means.max(axis=1) - masked_means.min(axis=1),
        "area_reflectance": np.trapz(masked_means, x=wavelengths[mask], axis=1),
    }

    if std_cols:
        std_matrix = df[std_cols].to_numpy(dtype=float)[:, mask]
        features["mean_std"] = std_matrix.mean(axis=1)
        features["median_std"] = np.median(std_matrix, axis=1)

    if mad_cols:
        mad_matrix = df[mad_cols].to_numpy(dtype=float)[:, mask]
        features["mean_mad"] = mad_matrix.mean(axis=1)
        features["median_mad"] = np.median(mad_matrix, axis=1)

    feature_df = df[group_columns].copy()
    for key, values in features.items():
        feature_df[key] = values
    feature_df["wavelength_lower_nm"] = float(wavelengths[mask].min())
    feature_df["wavelength_upper_nm"] = float(wavelengths[mask].max())
    feature_df["wavelength_step_nm"] = float(wavelengths[1] - wavelengths[0])
    feature_df["wavelength_point_count"] = int(mask.sum())
    return feature_df


def write_outputs(
    angles_df: pd.DataFrame,
    composites_df: pd.DataFrame,
    window_label: str,
    output_dir: Path,
) -> Tuple[Path, Path]:
    angle_path = output_dir / f"phase3_reflectance_features_angles_{window_label}.csv"
    composite_path = output_dir / f"phase3_reflectance_features_composites_{window_label}.csv"
    angles_df.to_csv(angle_path, index=False)
    composites_df.to_csv(composite_path, index=False)
    return angle_path, composite_path


def parse_window_arg(arg: str) -> Tuple[float, float, str]:
    bounds = arg.split(":")
    if len(bounds) != 2:
        raise ValueError(f"Window specification '{arg}' must be in the form lower:upper")
    lower = float(bounds[0])
    upper = float(bounds[1])
    if lower >= upper:
        raise ValueError(f"Window lower bound must be < upper bound (received {lower} >= {upper}).")
    label = f"{int(lower)}_{int(upper)}nm"
    return lower, upper, label


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Derive reflectance features for Phase 3 diagnostics.")
    parser.add_argument(
        "--manifest",
        default="analysis/phase3/phase3_data_manifest.json",
        help="Path to the Phase 3 loader manifest.",
    )
    parser.add_argument(
        "--canonical-root",
        default="canonical_dataset",
        help="Directory containing canonical reflectance tables.",
    )
    parser.add_argument(
        "--window",
        action="append",
        default=["320:480"],
        help="Wavelength window (lower:upper, nm). May be provided multiple times.",
    )
    parser.add_argument(
        "--output-dir",
        default="analysis/phase3",
        help="Directory to write derived feature tables.",
    )
    parser.add_argument(
        "--manifest-json",
        default="analysis/phase3/phase3_reflectance_feature_manifest.json",
        help="Path to write a manifest describing generated feature tables.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.manifest).resolve()
    manifest = load_manifest(manifest_path)
    wavelengths = extract_wavelength_axis(manifest)

    canonical_root = Path(args.canonical_root).resolve()
    angles, composites = load_reflectance_tables(canonical_root)

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    feature_manifest = {
        "manifest_source": str(manifest_path),
        "generated": [],
    }

    for window_spec in args.window:
        lower, upper, label = parse_window_arg(window_spec)
        mask = build_window_mask(wavelengths, lower, upper)
        angle_features = compute_features(
            angles,
            wavelengths,
            mask,
            group_columns=[
                "dose_id",
                "angle",
                "uva_mw_cm2",
                "uvb_mw_cm2",
                "samples_total",
                "samples_used",
                "trim_fraction",
            ],
        )
        composite_features = compute_features(
            composites,
            wavelengths,
            mask,
            group_columns=[
                "dose_id",
                "composite",
                "samples_12",
                "samples_6",
                "trim_fraction",
                "uva_mw_cm2",
                "uvb_mw_cm2",
            ],
        )
        angle_path, composite_path = write_outputs(angle_features, composite_features, label, output_dir)
        feature_manifest["generated"].append(
            {
                "window": {
                    "lower_nm": lower,
                    "upper_nm": upper,
                    "step_nm": float(wavelengths[1] - wavelengths[0]),
                    "point_count": int(mask.sum()),
                },
                "angles_csv": str(angle_path),
                "composites_csv": str(composite_path),
                "feature_columns": [c for c in angle_features.columns if c not in {"dose_id", "angle", "uva_mw_cm2", "uvb_mw_cm2", "samples_total", "samples_used", "trim_fraction", "wavelength_lower_nm", "wavelength_upper_nm", "wavelength_step_nm", "wavelength_point_count"}],
            }
        )
        print(f"[INFO] Generated features for {label}: angles -> {angle_path.name}, composites -> {composite_path.name}")

    manifest_path_out = Path(args.manifest_json).resolve()
    manifest_path_out.write_text(json.dumps(feature_manifest, indent=2))
    print(f"[INFO] Feature manifest written to {manifest_path_out}")


if __name__ == "__main__":
    main()
