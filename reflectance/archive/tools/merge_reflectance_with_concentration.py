#!/usr/bin/env python3
"""Prepare an auditable crosswalk between reflectance folders and concentration IDs.

The script expects:
  1. An aggregated reflectance CSV produced by aggregate_reflectance.py.
  2. The combined scytonemin concentration CSV (with sample_id metadata).
  3. A mapping CSV that links (treatment, sample_label) pairs to sample_id.

The mapping CSV must have columns:
    treatment,sample_label,sample_id
where:
    - treatment matches the reflectance naming (e.g., "treatment_1")
    - sample_label matches entries like "sample_A"
    - sample_id matches identifiers in Combined_Scytonemin_Concentrations.csv

If the mapping file is missing, the script will generate a template listing all
required (treatment, sample_label) combinations with blank sample_id values.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd


def load_reflectance(reflectance_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(reflectance_csv)
    required_cols = {"treatment", "sample_label", "angle"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Reflectance CSV missing columns: {missing}")
    return df


def load_concentrations(concentration_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(concentration_csv)
    if "sample_id" not in df.columns:
        raise ValueError("Concentration CSV must contain 'sample_id'.")
    return df


def load_mapping(mapping_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(mapping_csv)
    required = {"treatment", "sample_label", "sample_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Mapping file missing columns: {missing}")
    if df["sample_id"].isna().any():
        raise ValueError(
            "Mapping file contains blank sample_id entries. "
            "Fill them before merging."
        )
    duplicates = df.duplicated(subset=["treatment", "sample_label"])
    if duplicates.any():
        bad = df.loc[duplicates, ["treatment", "sample_label", "sample_id"]]
        raise ValueError(
            "Duplicate entries in mapping for treatment/sample_label:\n"
            f"{bad}"
        )
    return df


def generate_mapping_template(
    reflectance_df: pd.DataFrame, template_path: Path
) -> None:
    unique_pairs = (
        reflectance_df[["treatment", "sample_label"]]
        .drop_duplicates()
        .sort_values(["treatment", "sample_label"])
    )
    unique_pairs["sample_id"] = ""
    template_path.write_text(unique_pairs.to_csv(index=False))
    print(
        f"Mapping template written to {template_path}. "
        "Fill 'sample_id' values and rerun."
    )


def build_crosswalk(
    reflectance_df: pd.DataFrame,
    concentration_df: pd.DataFrame,
    mapping_df: pd.DataFrame,
) -> pd.DataFrame:
    unique_pairs = (
        reflectance_df[["treatment", "sample_label"]]
        .drop_duplicates()
        .sort_values(["treatment", "sample_label"])
    )
    crosswalk = unique_pairs.merge(mapping_df, on=["treatment", "sample_label"], how="left")
    if crosswalk["sample_id"].isna().any():
        missing = crosswalk[crosswalk["sample_id"].isna()][["treatment", "sample_label"]]
        raise ValueError(
            "Mapping incomplete. Missing sample_id assignments for:\n"
            f"{missing.to_string(index=False)}"
        )

    concentration_ids = set(concentration_df["sample_id"].unique())
    crosswalk["has_concentration"] = crosswalk["sample_id"].isin(concentration_ids)
    crosswalk["note"] = np.where(
        crosswalk["has_concentration"],
        "",
        "no matching concentration ID available (destructive sampling; treatment-level only)",
    )
    return crosswalk


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge aggregated reflectance spectra with concentration data."
    )
    parser.add_argument(
        "--reflectance-csv",
        type=Path,
        default=Path("aggregated_reflectance/reflectance_trimmed_stats.csv"),
        help="Path to aggregated reflectance CSV.",
    )
    parser.add_argument(
        "--concentration-csv",
        type=Path,
        default=Path("Combined_Scytonemin_Concentrations.csv"),
        help="Path to concentration CSV with sample_id column.",
    )
    parser.add_argument(
        "--mapping-csv",
        type=Path,
        default=Path("reflectance_sample_mapping.csv"),
        help="CSV mapping treatment/sample_label to sample_id.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("reflectance_concentration_crosswalk.csv"),
        help="Output CSV path for the reflectance-to-concentration crosswalk.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reflectance_df = load_reflectance(args.reflectance_csv)
    concentration_df = load_concentrations(args.concentration_csv)

    if not args.mapping_csv.exists():
        generate_mapping_template(reflectance_df, args.mapping_csv)
        return

    mapping_df = load_mapping(args.mapping_csv)
    crosswalk = build_crosswalk(reflectance_df, concentration_df, mapping_df)
    crosswalk.to_csv(args.output_csv, index=False)
    matches = crosswalk["has_concentration"].sum()
    total = len(crosswalk)
    print(
        f"Crosswalk written to {args.output_csv} "
        f"({total} rows; concentration coverage {matches}/{total} = {matches/total:.1%})."
    )
    missing_records = crosswalk[~crosswalk["has_concentration"]]
    if not missing_records.empty:
        missing_path = args.output_csv.with_name("crosswalk_missing_concentration_ids.csv")
        missing_records.to_csv(missing_path, index=False)
        print(
            f"WARNING: {len(missing_records)} entries lack concentration data. "
            f"Details saved to {missing_path}"
        )


if __name__ == "__main__":
    main()
