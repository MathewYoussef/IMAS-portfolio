#!/usr/bin/env python3
"""
Phase 3 data loader and validator.

This script centralises ingestion of the canonical dose-level tables that drive
the Phase 3 diagnostics.  It enforces a consistent dose ordering, checks the
expected Σ/Δ and angle metadata, and records wavelength axis metadata so
downstream notebooks can operate from a single source of truth.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd
import numpy as np

# Ensure project root (two levels up) is on the Python path so we can import
# shared modules when running from the analysis directory.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
import sys

if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from dose_metadata import dose_mapping

# Canonical dose ordering (low → high UVA)
DOSE_ORDER = {label: idx + 1 for idx, label in enumerate(sorted(dose_mapping().keys()))}

# Tables we rely on for Phase 3 diagnostics.
CANONICAL_TABLES = {
    "dose_reflectance_stats": "dose_reflectance_stats.csv",
    "dose_reflectance_composites": "dose_reflectance_composites.csv",
    "dose_reflectance_features_angles": "dose_reflectance_features_angles.csv",
    "dose_reflectance_features_composites": "dose_reflectance_features_composites.csv",
    "dose_dad_concentrations": "dose_dad_concentrations.csv",
    "dose_summary_chrom": "dose_summary.csv",
    "precision_weighted_treatment": "precision_weighted_concentrations_treatment.csv",
    "dose_level_canonical_summary": "dose_level_canonical_summary.csv",
}


@dataclass
class LoaderConfig:
    canonical_root: Path
    wavelength_grid: np.ndarray | None = None


def _sorted_by_dose(df: pd.DataFrame, dose_col: str = "dose_id") -> pd.DataFrame:
    missing = set(df[dose_col].unique()) - DOSE_ORDER.keys()
    if missing:
        raise ValueError(f"Unknown dose labels encountered: {sorted(missing)}")
    return df.sort_values(dose_col, key=lambda s: s.map(DOSE_ORDER))


def _validate_angles(df: pd.DataFrame) -> None:
    expected_angles = {"12Oclock", "6Oclock"}
    grouped = df.groupby("dose_id")["angle"].apply(set)
    for dose_id, seen in grouped.items():
        if seen != expected_angles:
            raise ValueError(f"{dose_id}: expected angles {expected_angles}, found {seen}")


def _validate_composites(df: pd.DataFrame) -> None:
    expected = {"Sigma", "Delta"}
    grouped = df.groupby("dose_id")["composite"].apply(set)
    for dose_id, seen in grouped.items():
        if seen != expected:
            raise ValueError(f"{dose_id}: expected composites {expected}, found {seen}")


def _load_wavelength_grid(path: Path) -> np.ndarray:
    grid = np.load(path)
    if grid.ndim != 1:
        raise ValueError(f"Wavelength grid must be 1-D array, got shape {grid.shape} from {path}")
    if grid.size < 2:
        raise ValueError("Wavelength grid must contain at least two points to infer spacing.")
    return grid


def load_tables(config: LoaderConfig) -> Dict[str, pd.DataFrame]:
    tables: Dict[str, pd.DataFrame] = {}
    for key, filename in CANONICAL_TABLES.items():
        path = config.canonical_root / filename
        if not path.exists():
            raise FileNotFoundError(f"Expected canonical table missing: {path}")
        df = pd.read_csv(path)
        if "dose_id" in df.columns:
            df = _sorted_by_dose(df, "dose_id")
        tables[key] = df

    # Validate reflectance metadata expectations.
    _validate_angles(tables["dose_reflectance_stats"])
    _validate_composites(tables["dose_reflectance_composites"])

    return tables


def summarise_loader(tables: Dict[str, pd.DataFrame], config: LoaderConfig) -> Dict[str, object]:
    summary: Dict[str, object] = {
        "canonical_root": str(config.canonical_root),
        "dose_order": list(DOSE_ORDER.keys()),
        "tables": {},
    }

    for key, df in tables.items():
        summary["tables"][key] = {
            "rows": int(df.shape[0]),
            "columns": list(df.columns),
        }

    reflectance_df = tables["dose_reflectance_stats"]
    spectral_cols: List[str] = [c for c in reflectance_df.columns if c.startswith("mean_")]
    if config.wavelength_grid is not None:
        grid = config.wavelength_grid
        if len(spectral_cols) != grid.size:
            raise ValueError(
                f"Wavelength grid length ({grid.size}) does not match spectral column count ({len(spectral_cols)})."
            )
        summary["wavelength_axis"] = {
            "start_nm": float(grid[0]),
            "step_nm": float(grid[1] - grid[0]),
            "count": int(grid.size),
        }
    else:
        summary["wavelength_axis"] = {
            "start_nm": None,
            "step_nm": None,
            "count": len(spectral_cols),
        }
    summary["angles"] = sorted(reflectance_df["angle"].unique())
    summary["composites"] = sorted(tables["dose_reflectance_composites"]["composite"].unique())

    return summary


def write_summary(summary: Dict[str, object], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Load and validate canonical tables for Phase 3 analysis.")
    parser.add_argument(
        "--canonical-root",
        default="canonical_dataset",
        help="Directory containing the canonical tables (default: canonical_dataset).",
    )
    parser.add_argument(
        "--wavelength-grid",
        type=Path,
        default=Path("denoised_full_run_staging") / "wavelength_grid.npy",
        help="Path to the wavelength grid (.npy) matching the reflectance spectra.",
    )
    parser.add_argument(
        "--summary-json",
        default="analysis/phase3/phase3_data_manifest.json",
        help="Path to write loader summary metadata (default: analysis/phase3/phase3_data_manifest.json).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    canonical_root = Path(args.canonical_root).resolve()
    wavelength_grid_path = args.wavelength_grid.resolve()
    wavelength_grid = None
    if wavelength_grid_path.exists():
        wavelength_grid = _load_wavelength_grid(wavelength_grid_path)
    else:
        print(f"[WARN] Wavelength grid not found at {wavelength_grid_path}; summary will omit axis metadata.")
    config = LoaderConfig(
        canonical_root=canonical_root,
        wavelength_grid=wavelength_grid,
    )
    tables = load_tables(config)
    summary = summarise_loader(tables, config)
    write_summary(summary, Path(args.summary_json).resolve())

    print("Loaded tables:")
    for key, info in summary["tables"].items():
        print(f"  {key}: {info['rows']} rows, {len(info['columns'])} columns")
    print(f"Wavelength axis: start={summary['wavelength_axis']['start_nm']} nm, "
          f"step={summary['wavelength_axis']['step_nm']} nm, "
          f"count={summary['wavelength_axis']['count']}")


if __name__ == "__main__":
    main()
