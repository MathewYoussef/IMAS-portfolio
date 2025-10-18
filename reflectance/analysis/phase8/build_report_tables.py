#!/usr/bin/env python3
"""Generate aggregated tables for the Phase 8 report."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

MANIFEST_PATH = Path("analysis/phase8/inputs/phase8_input_manifest.json")
OUTPUT_TABLE = Path("analysis/phase8/tables/comparison_overview.csv")


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text())


def build_table() -> pd.DataFrame:
    manifest = load_manifest()
    kind_modality_summary = pd.read_csv(manifest["phase7"]["kind_modality_summary"]["path"])
    modality_summary = pd.read_csv(manifest["phase7"]["modality_summary"]["path"])

    key_cols = [
        "window",
        "kind",
        "modality",
        "dcor_permutation_p",
        "rv_permutation_p",
        "procrustes_permutation_p",
        "baseline_rmse_loocv",
        "baseline_mae_loocv",
        "multiblock_joint_components",
    ]
    extra_cols = [
        "feature_set",
        "model_type",
        "band_label",
        "band_lower_nm",
        "band_upper_nm",
        "component",
    ]
    available_cols = key_cols + [col for col in extra_cols if col in kind_modality_summary.columns]
    comparison = kind_modality_summary[available_cols].copy()
    if "feature_set" in comparison.columns:
        comparison = comparison[comparison["feature_set"] == "phase9"]

    modality_pvals = modality_summary[[
        "modality",
        "pearson_perm_p",
        "spearman_perm_p",
        "kendall_perm_p",
    ]]

    comparison = comparison.merge(modality_pvals, on="modality", how="left")
    comparison = comparison.sort_values(["modality", "kind"])
    return comparison


def main() -> None:
    table = build_table()
    OUTPUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUTPUT_TABLE, index=False)

    manifest = load_manifest()
    manifest.setdefault("phase8", {})
    manifest["phase8"]["comparison_overview"] = {
        "path": str(OUTPUT_TABLE),
        "source": "analysis/phase8/build_report_tables.py",
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2))
    print(f"[INFO] Wrote {OUTPUT_TABLE}")


if __name__ == "__main__":
    main()
