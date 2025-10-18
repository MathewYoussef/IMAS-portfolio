#!/usr/bin/env python3
"""Back-fill UVA/UVB columns into existing aggregation outputs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import pandas as pd

from dose_metadata import DoseRecord, attach_dose_metadata, iter_dose_records


ROOT = Path("aggregated_reflectance")


def iter_records() -> Iterable[DoseRecord]:
    return iter_dose_records()


def assign_by_dose(df: pd.DataFrame, dose_col: str) -> pd.DataFrame:
    if dose_col not in df.columns:
        raise KeyError(f"Column '{dose_col}' not found in dataframe.")
    df["uva_mw_cm2"] = df[dose_col].map(lambda label: attach_dose_metadata(label).uva_mw_cm2)
    df["uvb_mw_cm2"] = df[dose_col].map(lambda label: attach_dose_metadata(label).uvb_mw_cm2)
    return df


def treatment_to_dose(label: str) -> str:
    if not label.startswith("treatment_"):
        raise ValueError(f"Unexpected treatment label {label!r}")
    ordinal = int(label.split("_", 1)[1])
    return f"dose_{7 - ordinal}"


def sample_to_dose(sample_id: str) -> str:
    ordinal = int(sample_id[0])
    return f"dose_{7 - ordinal}"


def update_csv(path: Path, dose_col: str | None = "dose_id", allow_missing: bool = False) -> None:
    df = pd.read_csv(path)
    if dose_col and dose_col not in df.columns:
        if allow_missing:
            return
        raise KeyError(f"{dose_col} missing from {path}")

    if dose_col is None:
        raise ValueError("dose_col must be provided for assignment.")

    assign_by_dose(df, dose_col)
    df.to_csv(path, index=False)


def update_precision_weighted_sample(path: Path) -> None:
    df = pd.read_csv(path)
    if "dose_id" not in df.columns:
        df["dose_id"] = df["sample_id"].map(sample_to_dose)
    assign_by_dose(df, "dose_id")
    df.to_csv(path, index=False)


def update_precision_weighted_treatment(path: Path) -> None:
    df = pd.read_csv(path)
    if "dose_id" not in df.columns:
        df["dose_id"] = df["treatment"].map(treatment_to_dose)
    assign_by_dose(df, "dose_id")
    df.to_csv(path, index=False)


def write_meta(path: Path) -> None:
    data = json.loads(path.read_text())
    data["uva_uvb_by_dose"] = {
        record.label: {
            "uva_mw_cm2": record.uva_mw_cm2,
            "uvb_mw_cm2": record.uvb_mw_cm2,
        }
        for record in iter_records()
    }
    path.write_text(json.dumps(data, indent=2))


def main() -> None:
    update_csv(ROOT / "reflectance_trimmed_stats.csv")
    update_csv(ROOT / "dose_reflectance_stats.csv")
    update_csv(ROOT / "dose_reflectance_composites.csv")
    update_csv(ROOT / "dose_dad_concentrations.csv")
    update_precision_weighted_sample(ROOT / "precision_weighted_concentrations.csv")
    update_precision_weighted_treatment(ROOT / "precision_weighted_concentrations_treatment.csv")

    write_meta(ROOT / "aggregation_metadata.json")
    write_meta(ROOT / "dose_reflectance_metadata.json")
    write_meta(ROOT / "dose_dad_metadata.json")
    write_meta(ROOT / "precision_weighted_concentration_metadata.json")


if __name__ == "__main__":
    main()
