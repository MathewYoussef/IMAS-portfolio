#!/usr/bin/env python3
"""
Assemble Phase 3 concentration targets aligned with reflectance features.

The script reads the Phase 3 manifests, canonical concentration tables, and
derived reflectance feature files, then produces long-form tables that join
each wavelength window to chromatogram, DAD, and latent concentration metrics.

Outputs (per window):
    - phase3_join_reflectance_chrom_{label}.csv
    - phase3_join_reflectance_dad_{label}.csv
    - phase3_join_reflectance_latent_{label}.csv
    - A manifest JSON listing sources, renames, and output paths.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd


@dataclass
class ManifestPaths:
    phase3_manifest: Path
    feature_manifest: Path
    canonical_root: Path
    output_dir: Path
    output_manifest: Path


CONCENTRATION_RENAMES = {
    "chrom_total_mg_per_gDW_trimmed_mean": "chrom_total_mean",
    "chrom_total_mg_per_gDW_trimmed_sd": "chrom_total_sd",
    "chrom_total_mg_per_gDW_ci_low": "chrom_total_ci_low",
    "chrom_total_mg_per_gDW_ci_high": "chrom_total_ci_high",
    "chrom_oxidized_mg_per_gDW_trimmed_mean": "chrom_oxidized_mean",
    "chrom_oxidized_mg_per_gDW_trimmed_sd": "chrom_oxidized_sd",
    "chrom_oxidized_mg_per_gDW_ci_low": "chrom_oxidized_ci_low",
    "chrom_oxidized_mg_per_gDW_ci_high": "chrom_oxidized_ci_high",
    "chrom_reduced_mg_per_gDW_trimmed_mean": "chrom_reduced_mean",
    "chrom_reduced_mg_per_gDW_trimmed_sd": "chrom_reduced_sd",
    "chrom_reduced_mg_per_gDW_ci_low": "chrom_reduced_ci_low",
    "chrom_reduced_mg_per_gDW_ci_high": "chrom_reduced_ci_high",
    "dad_total_mg_per_gDW_trimmed_mean": "dad_total_mean",
    "dad_total_mg_per_gDW_trimmed_sd": "dad_total_sd",
    "dad_total_mg_per_gDW_ci_low": "dad_total_ci_low",
    "dad_total_mg_per_gDW_ci_high": "dad_total_ci_high",
    "dad_oxidized_mg_per_gDW_trimmed_mean": "dad_oxidized_mean",
    "dad_oxidized_mg_per_gDW_trimmed_sd": "dad_oxidized_sd",
    "dad_oxidized_mg_per_gDW_ci_low": "dad_oxidized_ci_low",
    "dad_oxidized_mg_per_gDW_ci_high": "dad_oxidized_ci_high",
    "dad_reduced_mg_per_gDW_trimmed_mean": "dad_reduced_mean",
    "dad_reduced_mg_per_gDW_trimmed_sd": "dad_reduced_sd",
    "dad_reduced_mg_per_gDW_ci_low": "dad_reduced_ci_low",
    "dad_reduced_mg_per_gDW_ci_high": "dad_reduced_ci_high",
}


def load_json(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text())


def load_canonical_tables(root: Path) -> Dict[str, pd.DataFrame]:
    tables = {
        "chrom": pd.read_csv(root / "dose_summary.csv"),
        "dad": pd.read_csv(root / "dose_dad_concentrations.csv"),
        "latent": pd.read_csv(root / "precision_weighted_concentrations_treatment.csv"),
        "summary": pd.read_csv(root / "dose_level_canonical_summary.csv"),
    }
    return tables


def harmonise_chrom(df: pd.DataFrame) -> pd.DataFrame:
    subset = df[
        [
            "dose_id",
            "uva_mw_cm2",
            "uvb_mw_cm2",
            "chrom_total_mg_per_gDW_trimmed_mean",
            "chrom_total_mg_per_gDW_trimmed_sd",
            "chrom_total_mg_per_gDW_ci_low",
            "chrom_total_mg_per_gDW_ci_high",
            "chrom_oxidized_mg_per_gDW_trimmed_mean",
            "chrom_oxidized_mg_per_gDW_trimmed_sd",
            "chrom_oxidized_mg_per_gDW_ci_low",
            "chrom_oxidized_mg_per_gDW_ci_high",
            "chrom_reduced_mg_per_gDW_trimmed_mean",
            "chrom_reduced_mg_per_gDW_trimmed_sd",
            "chrom_reduced_mg_per_gDW_ci_low",
            "chrom_reduced_mg_per_gDW_ci_high",
        ]
    ].copy()
    subset = subset.rename(columns=CONCENTRATION_RENAMES)
    return subset


def harmonise_dad(df: pd.DataFrame) -> pd.DataFrame:
    subset = df[
        [
            "dose_id",
            "uva_mw_cm2",
            "uvb_mw_cm2",
            "dad_total_mg_per_gDW_mean_trimmed",
            "dad_total_mg_per_gDW_std_trimmed",
            "dad_total_mg_per_gDW_median",
            "dad_total_mg_per_gDW_mad",
            "dad_total_mg_per_gDW_n_used",
            "dad_oxidized_mg_per_gDW_mean_trimmed",
            "dad_oxidized_mg_per_gDW_std_trimmed",
            "dad_oxidized_mg_per_gDW_median",
            "dad_oxidized_mg_per_gDW_mad",
            "dad_oxidized_mg_per_gDW_n_used",
            "dad_reduced_mg_per_gDW_mean_trimmed",
            "dad_reduced_mg_per_gDW_std_trimmed",
            "dad_reduced_mg_per_gDW_median",
            "dad_reduced_mg_per_gDW_mad",
            "dad_reduced_mg_per_gDW_n_used",
        ]
    ].copy()
    subset = subset.rename(
        columns={
            "dad_total_mg_per_gDW_mean_trimmed": "dad_total_mean",
            "dad_total_mg_per_gDW_std_trimmed": "dad_total_sd",
            "dad_total_mg_per_gDW_median": "dad_total_median",
            "dad_total_mg_per_gDW_mad": "dad_total_mad",
            "dad_total_mg_per_gDW_n_used": "dad_total_n_used",
            "dad_oxidized_mg_per_gDW_mean_trimmed": "dad_oxidized_mean",
            "dad_oxidized_mg_per_gDW_std_trimmed": "dad_oxidized_sd",
            "dad_oxidized_mg_per_gDW_median": "dad_oxidized_median",
            "dad_oxidized_mg_per_gDW_mad": "dad_oxidized_mad",
            "dad_oxidized_mg_per_gDW_n_used": "dad_oxidized_n_used",
            "dad_reduced_mg_per_gDW_mean_trimmed": "dad_reduced_mean",
            "dad_reduced_mg_per_gDW_std_trimmed": "dad_reduced_sd",
            "dad_reduced_mg_per_gDW_median": "dad_reduced_median",
            "dad_reduced_mg_per_gDW_mad": "dad_reduced_mad",
            "dad_reduced_mg_per_gDW_n_used": "dad_reduced_n_used",
        }
    )
    return subset


def harmonise_latent(df: pd.DataFrame) -> pd.DataFrame:
    subset = df[
        [
            "dose_id",
            "uva_mw_cm2",
            "uvb_mw_cm2",
            "treatment",
            "trim_fraction",
            "sample_count",
            "total_latent_mean_trimmed",
            "total_latent_std_trimmed",
            "total_latent_n_used",
            "total_latent_se_median",
            "oxidized_latent_mean_trimmed",
            "oxidized_latent_std_trimmed",
            "oxidized_latent_n_used",
            "oxidized_latent_se_median",
            "reduced_latent_mean_trimmed",
            "reduced_latent_std_trimmed",
            "reduced_latent_n_used",
            "reduced_latent_se_median",
        ]
    ].copy()
    subset = subset.rename(
        columns={
            "total_latent_mean_trimmed": "latent_total_mean",
            "total_latent_std_trimmed": "latent_total_sd",
            "total_latent_n_used": "latent_total_n_used",
            "total_latent_se_median": "latent_total_se_median",
            "oxidized_latent_mean_trimmed": "latent_oxidized_mean",
            "oxidized_latent_std_trimmed": "latent_oxidized_sd",
            "oxidized_latent_n_used": "latent_oxidized_n_used",
            "oxidized_latent_se_median": "latent_oxidized_se_median",
            "reduced_latent_mean_trimmed": "latent_reduced_mean",
            "reduced_latent_std_trimmed": "latent_reduced_sd",
            "reduced_latent_n_used": "latent_reduced_n_used",
            "reduced_latent_se_median": "latent_reduced_se_median",
        }
    )
    return subset


def reshape_reflectance(angle_df: pd.DataFrame, composite_df: pd.DataFrame) -> pd.DataFrame:
    angle_subset = angle_df.rename(columns={"angle": "kind"})
    composite_subset = composite_df.rename(columns={"composite": "kind"})
    composite_subset["samples_total"] = composite_subset["samples_12"]
    composite_subset["samples_used"] = composite_subset["samples_12"]

    shared_cols = [
        "dose_id",
        "kind",
        "uva_mw_cm2",
        "uvb_mw_cm2",
        "mean_reflectance",
        "median_reflectance",
        "min_reflectance",
        "max_reflectance",
        "range_reflectance",
        "area_reflectance",
        "mean_std",
        "median_std",
        "mean_mad",
        "median_mad",
        "wavelength_lower_nm",
        "wavelength_upper_nm",
        "wavelength_step_nm",
        "wavelength_point_count",
    ]

    angle_keep = shared_cols + ["samples_total", "samples_used", "trim_fraction"]
    composite_keep = shared_cols + ["samples_12", "samples_6", "trim_fraction", "samples_total", "samples_used"]

    angle_trimmed = angle_subset[[c for c in angle_keep if c in angle_subset.columns]].copy()
    composite_trimmed = composite_subset[[c for c in composite_keep if c in composite_subset.columns]].copy()
    composite_trimmed["samples_total"] = composite_subset["samples_12"]
    composite_trimmed["samples_used"] = composite_subset["samples_12"]

    combined = pd.concat([angle_trimmed, composite_trimmed], ignore_index=True)
    return combined


def join_reflectance_with_target(reflectance: pd.DataFrame, target: pd.DataFrame) -> pd.DataFrame:
    merged = reflectance.merge(
        target,
        on=["dose_id", "uva_mw_cm2", "uvb_mw_cm2"],
        how="inner",
        validate="many_to_one",
    )
    return merged


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assemble Phase 3 reflectance vs concentration join tables.")
    parser.add_argument(
        "--phase3-manifest",
        default="analysis/phase3/phase3_data_manifest.json",
        help="Path to the Phase 3 loader manifest.",
    )
    parser.add_argument(
        "--feature-manifest",
        default="analysis/phase3/phase3_reflectance_feature_manifest.json",
        help="Manifest describing derived reflectance feature tables.",
    )
    parser.add_argument(
        "--canonical-root",
        default="canonical_dataset",
        help="Directory containing canonical concentration tables.",
    )
    parser.add_argument(
        "--output-dir",
        default="analysis/phase3",
        help="Directory for join outputs.",
    )
    parser.add_argument(
        "--output-manifest",
        default="analysis/phase3/phase3_concentration_target_manifest.json",
        help="Path to write manifest describing join outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = ManifestPaths(
        phase3_manifest=Path(args.phase3_manifest).resolve(),
        feature_manifest=Path(args.feature_manifest).resolve(),
        canonical_root=Path(args.canonical_root).resolve(),
        output_dir=Path(args.output_dir).resolve(),
        output_manifest=Path(args.output_manifest).resolve(),
    )

    phase3_manifest = load_json(paths.phase3_manifest)
    feature_manifest = load_json(paths.feature_manifest)
    canonical_tables = load_canonical_tables(paths.canonical_root)

    chrom_df = harmonise_chrom(canonical_tables["chrom"])
    dad_df = harmonise_dad(canonical_tables["dad"])
    latent_df = harmonise_latent(canonical_tables["latent"])

    paths.output_dir.mkdir(parents=True, exist_ok=True)

    output_records: List[Dict[str, object]] = []

    for entry in feature_manifest.get("generated", []):
        angles_path = Path(entry["angles_csv"]).resolve()
        composites_path = Path(entry["composites_csv"]).resolve()
        if not angles_path.exists() or not composites_path.exists():
            print(f"[WARN] Reflectance feature files missing for window entry {entry}; skipping.")
            continue
        window = entry["window"]
        label = f"{int(window['lower_nm'])}_{int(window['upper_nm'])}nm"

        angle_df = pd.read_csv(angles_path)
        composite_df = pd.read_csv(composites_path)
        reflectance_long = reshape_reflectance(angle_df, composite_df)

        chrom_join = join_reflectance_with_target(reflectance_long, chrom_df)
        dad_join = join_reflectance_with_target(reflectance_long, dad_df)
        latent_join = join_reflectance_with_target(reflectance_long, latent_df)

        chrom_out = paths.output_dir / f"phase3_join_reflectance_chrom_{label}.csv"
        dad_out = paths.output_dir / f"phase3_join_reflectance_dad_{label}.csv"
        latent_out = paths.output_dir / f"phase3_join_reflectance_latent_{label}.csv"

        chrom_join.to_csv(chrom_out, index=False)
        dad_join.to_csv(dad_out, index=False)
        latent_join.to_csv(latent_out, index=False)

        output_records.append(
            {
                "window": window,
                "angles_csv": str(angles_path),
                "composites_csv": str(composites_path),
                "chrom_join": str(chrom_out),
                "dad_join": str(dad_out),
                "latent_join": str(latent_out),
            }
        )

        print(f"[INFO] Joined window {label}: chrom -> {chrom_out.name}, dad -> {dad_out.name}, latent -> {latent_out.name}")

    manifest_out = {
        "phase3_manifest": str(paths.phase3_manifest),
        "feature_manifest": str(paths.feature_manifest),
        "canonical_root": str(paths.canonical_root),
        "outputs": output_records,
    }
    paths.output_manifest.write_text(json.dumps(manifest_out, indent=2))
    print(f"[INFO] Join manifest written to {paths.output_manifest}")


if __name__ == "__main__":
    main()
