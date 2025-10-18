#!/usr/bin/env python3
"""
Build consolidated summary tables for the Phase 3–6 comparison workspace.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import pandas as pd

MANIFEST_PATH = Path("analysis/phase7/inputs.json")


def load_manifest() -> Dict[str, object]:
    return json.loads(MANIFEST_PATH.read_text())


def build_modality_summary(manifest: Dict[str, object]) -> pd.DataFrame:
    corr_path = Path(manifest["phase3"]["correlations"]["path"])
    reg_path = Path(manifest["phase3"]["regression_summary"]["path"])

    corr_df = pd.read_json(corr_path)
    corr_df = corr_df[corr_df["metric"].str.contains("_total")]
    corr_df = corr_df[["modality", "pearson_perm_p", "spearman_perm_p", "kendall_perm_p"]]

    reg_df = pd.read_csv(reg_path)
    reg_df = reg_df[reg_df["metric"].str.contains("_total")]
    reg_df = reg_df[["modality", "metric", "slope", "p_value"]]

    modality_summary = corr_df.merge(reg_df, on="modality", how="left")
    modality_summary = modality_summary.rename(columns={"metric": "phase3_metric"})
    return modality_summary


def build_kind_modality_summary(manifest: Dict[str, object]) -> pd.DataFrame:
    dcor_path = Path(manifest["phase4"]["dcor_summary"]["path"])
    congruence_path = Path(manifest["phase4"]["congruence_summary"]["path"])
    multiblock_path = Path(manifest["phase5"]["multiblock_summary"]["path"])
    dcor_df = pd.read_csv(dcor_path)
    if "window" in dcor_df.columns:
        dcor_df = dcor_df.rename(columns={"statistic": "dcor_statistic", "permutation_p": "dcor_permutation_p"})
    else:
        dcor_df["window"] = "320_480nm"
        dcor_df = dcor_df.rename(columns={"statistic": "dcor_statistic", "permutation_p": "dcor_permutation_p"})

    congruence_df = pd.read_csv(congruence_path)
    congruence_cols = [
        "window",
        "kind",
        "modality",
        "procrustes_permutation_p",
        "rv_value",
        "rv_permutation_p",
    ]
    congruence_df = congruence_df[congruence_cols]

    merged = dcor_df.merge(congruence_df, on=["window", "kind", "modality"], how="outer")

    multiblock_df = pd.read_csv(multiblock_path)
    multiblock_df = multiblock_df.rename(
        columns={
            "blocks_used": "multiblock_blocks_used",
            "used_joint_components": "multiblock_joint_components",
            "warning": "multiblock_warning",
        }
    )[["window", "modality", "multiblock_blocks_used", "multiblock_joint_components", "multiblock_warning"]]

    merged = merged.merge(multiblock_df, on=["window", "modality"], how="left")

    merged["component"] = "total"

    baseline_loocv_path = Path(manifest["phase6"]["baseline_loocv"]["path"])
    baseline_summary_path = Path(manifest["phase6"]["baseline_summary"]["path"])

    baseline_loocv_df = pd.read_csv(baseline_loocv_path)
    baseline_summary_df = pd.read_csv(baseline_summary_path)[["target", "kind", "band_label"]]
    baseline_df = baseline_loocv_df.merge(baseline_summary_df, on=["target", "kind"], how="left")
    target_modality_component = {
        "chrom_total": ("chrom", "total"),
        "chrom_oxidized": ("chrom", "oxidized"),
        "chrom_reduced": ("chrom", "reduced"),
        "dad_total": ("dad", "total"),
        "dad_oxidized": ("dad", "oxidized"),
        "dad_reduced": ("dad", "reduced"),
        "latent_total": ("latent", "total"),
    }
    baseline_df["modality"] = baseline_df["target"].map(lambda t: target_modality_component.get(t, (None, None))[0])
    baseline_df["component"] = baseline_df["target"].map(lambda t: target_modality_component.get(t, (None, None))[1])
    baseline_df = baseline_df.rename(
        columns={"rmse_loocv": "baseline_rmse_loocv", "mae_loocv": "baseline_mae_loocv"}
    )[["modality", "kind", "band_label", "component", "baseline_rmse_loocv", "baseline_mae_loocv"]]

    baseline_keys = baseline_df[["modality", "kind", "band_label", "component"]].drop_duplicates()
    merged = merged.merge(baseline_keys, on=["modality", "kind", "component"], how="left")
    merged = merged.merge(baseline_df, on=["modality", "kind", "band_label", "component"], how="left")

    merged["feature_set"] = "phase6"
    merged["model_type"] = "linear"
    merged["band_lower_nm"] = pd.NA
    merged["band_upper_nm"] = pd.NA

    # Attach advanced baselines if available.
    phase9_config = manifest.get("phase9", {})
    advanced_path = phase9_config.get("advanced_baseline_summary", {}).get("path")
    advanced_rows_df: pd.DataFrame | None = None
    if advanced_path and Path(advanced_path).exists():
        advanced_df = pd.read_csv(advanced_path)
        if not advanced_df.empty:
            preferred = advanced_df[advanced_df["model_type"] == "ridge"].copy()
            preferred["modality"] = preferred["target"].map(lambda t: target_modality_component.get(t, (None, None))[0])
            preferred["component"] = preferred["target"].map(lambda t: target_modality_component.get(t, (None, None))[1])
            preferred = preferred.rename(
                columns={
                    "rmse_loocv": "baseline_rmse_loocv",
                    "mae_loocv": "baseline_mae_loocv",
                }
            )
            preferred["feature_set"] = "phase9"
            preferred = preferred[
                [
                    "modality",
                    "kind",
                    "band_label",
                    "band_lower_nm",
                    "band_upper_nm",
                    "component",
                    "model_type",
                    "feature_set",
                    "baseline_rmse_loocv",
                    "baseline_mae_loocv",
                ]
            ]

            metrics_base = (
                merged.drop(columns=["band_label", "baseline_rmse_loocv", "baseline_mae_loocv", "feature_set", "model_type", "band_lower_nm", "band_upper_nm"])
                .drop_duplicates()
            )
            metrics_base = metrics_base.drop(columns=["component"], errors="ignore")
            advanced_rows_df = metrics_base.merge(preferred, on=["modality", "kind"], how="inner")

    if advanced_rows_df is not None and not advanced_rows_df.empty:
        merged = pd.concat([merged, advanced_rows_df], ignore_index=True, sort=False)

    return merged


def main() -> None:
    manifest = load_manifest()

    modality_summary = build_modality_summary(manifest)
    modality_summary_path = Path("analysis/phase7/modality_summary.csv")
    modality_summary.to_csv(modality_summary_path, index=False)

    kind_modality_summary = build_kind_modality_summary(manifest)
    kind_modality_summary_path = Path("analysis/phase7/kind_modality_summary.csv")
    kind_modality_summary.to_csv(kind_modality_summary_path, index=False)

    manifest_update = json.loads(MANIFEST_PATH.read_text())
    manifest_update.setdefault("phase7", {})
    manifest_update["phase7"]["modality_summary"] = {
        "path": str(modality_summary_path),
        "source": "analysis/phase7/build_summary.py",
    }
    manifest_update["phase7"]["kind_modality_summary"] = {
        "path": str(kind_modality_summary_path),
        "source": "analysis/phase7/build_summary.py",
    }
    MANIFEST_PATH.write_text(json.dumps(manifest_update, indent=2))

    print(f"[INFO] Wrote {modality_summary_path}")
    print(f"[INFO] Wrote {kind_modality_summary_path}")


if __name__ == "__main__":
    main()
