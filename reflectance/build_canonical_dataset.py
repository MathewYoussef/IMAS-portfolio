#!/usr/bin/env python3
"""
Build a harmonised, canonical dataset from legacy aggregation outputs.

The script reads the existing tables under ``archive/aggregated_reflectance/``
without modifying them, normalises identifiers and dose metadata, and writes a
clean bundle of CSV/JSON artefacts to the directory specified via ``--output``
(default ``canonical_dataset``). Downstream analyses should use the canonical
outputs.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from dose_metadata import DoseRecord, attach_dose_metadata, iter_dose_records


# Canonical inputs now live under archive/aggregated_reflectance after the
# harmonisation. Keep the path centralised so audit scripts/tests share it.
SOURCE_ROOT = Path("archive") / "aggregated_reflectance"


def iter_doses() -> Iterable[DoseRecord]:
    """Yield canonical dose records in ascending UVA order."""
    return iter_dose_records()


def dose_order(labels: Sequence[str]) -> list[int]:
    """Return numerical ordering for strings like ``dose_1``."""
    return [int(label.split("_", 1)[1]) for label in labels]


def sort_by_dose(df: pd.DataFrame, column: str = "dose_id") -> pd.DataFrame:
    """Sort dataframe by the ordinal encoded in ``column`` (dose_1 … dose_6)."""
    return df.sort_values(column, key=lambda s: pd.Series(dose_order(s)))


def add_physical_columns(df: pd.DataFrame, dose_col: str) -> pd.DataFrame:
    """Attach UVA/UVB metadata based on the provided dose label column."""
    df = df.copy()
    if "dose_id" not in df.columns:
        df["dose_id"] = df[dose_col]
    df["uva_mw_cm2"] = df["dose_id"].map(lambda label: attach_dose_metadata(label).uva_mw_cm2)
    df["uvb_mw_cm2"] = df["dose_id"].map(lambda label: attach_dose_metadata(label).uvb_mw_cm2)
    return df


def label_from_uv(uva: float, uvb: float, tol: float = 1e-6) -> str:
    """Identify canonical dose label for the given UVA/UVB intensities."""
    for record in iter_doses():
        if math.isclose(uva, record.uva_mw_cm2, abs_tol=tol) and math.isclose(uvb, record.uvb_mw_cm2, abs_tol=tol):
            return record.label
    raise ValueError(f"No canonical dose matches UVA={uva}, UVB={uvb}")


def canonical_dose_summary(output_dir: Path) -> None:
    """Normalise legacy dose_summary.csv (numeric IDs) to canonical labels."""
    src = SOURCE_ROOT / "dose_summary.csv"
    df = pd.read_csv(src)
    df = df.rename(columns={"uva": "uva_mw_cm2", "uvb": "uvb_mw_cm2"})
    df["uva_mw_cm2"] = df["uva_mw_cm2"].astype(float)
    df["uvb_mw_cm2"] = df["uvb_mw_cm2"].astype(float)
    df["dose_id"] = df.apply(lambda row: label_from_uv(row["uva_mw_cm2"], row["uvb_mw_cm2"]), axis=1)
    columns = ["dose_id", "uva_mw_cm2", "uvb_mw_cm2"] + [
        c for c in df.columns if c not in {"dose_id", "uva_mw_cm2", "uvb_mw_cm2"}
    ]
    df = sort_by_dose(df)
    df[columns].to_csv(output_dir / "dose_summary.csv", index=False)


def canonical_table(
    name: str,
    output_dir: Path,
    dose_col: str = "dose_id",
    extra_drop: Iterable[str] | None = None,
    rename: dict[str, str] | None = None,
) -> None:
    """Copy a CSV to the canonical directory with harmonised dose metadata."""
    src = SOURCE_ROOT / name
    df = pd.read_csv(src)
    if rename:
        df = df.rename(columns=rename)
    df = add_physical_columns(df, dose_col)
    drop = set(extra_drop or [])
    df = df[[c for c in df.columns if c not in drop]]
    columns = ["dose_id", "uva_mw_cm2", "uvb_mw_cm2"] + [
        c for c in df.columns if c not in {"dose_id", "uva_mw_cm2", "uvb_mw_cm2"}
    ]
    df = sort_by_dose(df)
    df[columns].to_csv(output_dir / name, index=False)


def canonical_precision_weighted_sample(output_dir: Path) -> None:
    src = SOURCE_ROOT / "precision_weighted_concentrations.csv"
    df = pd.read_csv(src)
    if "dose_id" not in df.columns:
        df["dose_id"] = df["sample_id"].map(lambda sid: f"dose_{7 - int(sid[0])}")
    df = add_physical_columns(df, "dose_id")
    columns = ["sample_id", "dose_id", "uva_mw_cm2", "uvb_mw_cm2"] + [
        c for c in df.columns if c not in {"sample_id", "dose_id", "uva_mw_cm2", "uvb_mw_cm2"}
    ]
    df[columns].to_csv(output_dir / "precision_weighted_concentrations.csv", index=False)


def canonical_precision_weighted_treatment(output_dir: Path) -> None:
    src = SOURCE_ROOT / "precision_weighted_concentrations_treatment.csv"
    df = pd.read_csv(src)
    if "dose_id" not in df.columns:
        df["dose_id"] = df["treatment"].map(lambda label: f"dose_{7 - int(label.split('_', 1)[1])}")
    df = add_physical_columns(df, "dose_id")
    columns = ["dose_id", "uva_mw_cm2", "uvb_mw_cm2"] + [
        c for c in df.columns if c not in {"dose_id", "uva_mw_cm2", "uvb_mw_cm2"}
    ]
    df = sort_by_dose(df)
    df[columns].to_csv(output_dir / "precision_weighted_concentrations_treatment.csv", index=False)


def canonical_reflectance_dad_summary(output_dir: Path) -> None:
    src = SOURCE_ROOT / "dose_reflectance_dad_summary.csv"
    df = pd.read_csv(src)
    if not (df["dose_id_reflectance"] == df["dose_id_dad"]).all():
        raise ValueError("Mismatch between reflectance and DAD dose IDs in summary.")
    df = df.rename(columns={"dose_id_reflectance": "dose_id"}).drop(columns=["dose_id_dad"])
    df = add_physical_columns(df, "dose_id")
    columns = ["dose_id", "kind", "uva_mw_cm2", "uvb_mw_cm2"] + [
        c for c in df.columns if c not in {"dose_id", "kind", "uva_mw_cm2", "uvb_mw_cm2"}
    ]
    df = sort_by_dose(df)
    df[columns].to_csv(output_dir / "dose_reflectance_dad_summary.csv", index=False)


def copy_json(src: Path, dest: Path) -> None:
    data = json.loads(src.read_text())
    dest.write_text(json.dumps(data, indent=2))


def build_canonical_dataset(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    canonical_dose_summary(output_dir)
    canonical_table("reflectance_trimmed_stats.csv", output_dir)
    canonical_table("dose_reflectance_stats.csv", output_dir)
    canonical_table("dose_reflectance_composites.csv", output_dir)
    canonical_table("dose_reflectance_features_angles.csv", output_dir)
    canonical_table("dose_reflectance_features_composites.csv", output_dir)
    canonical_table("dose_dad_concentrations.csv", output_dir)
    canonical_precision_weighted_sample(output_dir)
    canonical_precision_weighted_treatment(output_dir)
    canonical_reflectance_dad_summary(output_dir)

    copy_json(
        SOURCE_ROOT / "dose_reflectance_dad_correlations.json",
        output_dir / "dose_reflectance_dad_correlations.json",
    )

    merged = (
        pd.read_csv(output_dir / "dose_summary.csv")
        .merge(
            pd.read_csv(output_dir / "dose_dad_concentrations.csv")[
                ["dose_id", "dad_total_mg_per_gDW_mean_trimmed"]
            ],
            on="dose_id",
            how="left",
        )
        .merge(
            pd.read_csv(output_dir / "dose_reflectance_features_composites.csv")
            .query("kind == 'Sigma'")[["dose_id", "mean_reflectance", "area_reflectance"]],
            on="dose_id",
            how="left",
        )
        .merge(
            pd.read_csv(output_dir / "precision_weighted_concentrations_treatment.csv")[
                ["dose_id", "total_latent_mean_trimmed", "total_latent_std_trimmed"]
            ],
            on="dose_id",
            how="left",
        )
    )
    merged = sort_by_dose(merged)
    merged.to_csv(output_dir / "dose_level_canonical_summary.csv", index=False)

    metadata = {
        "source_root": str(SOURCE_ROOT.resolve()),
        "output_root": str(output_dir),
        "files": sorted(p.name for p in output_dir.glob("*") if p.is_file()),
        "uva_uvb_lookup": {
            record.label: {
                "uva_mw_cm2": record.uva_mw_cm2,
                "uvb_mw_cm2": record.uvb_mw_cm2,
            }
            for record in iter_doses()
        },
    }
    (output_dir / "README.json").write_text(json.dumps(metadata, indent=2))


def archive_legacy(output_dir: Path, archive_root: Path) -> None:
    """Move legacy aggregation directory into an archive for safekeeping."""
    archive_root.mkdir(parents=True, exist_ok=True)
    target = archive_root / SOURCE_ROOT.name
    if target.exists():
        raise FileExistsError(f"Archive target {target} already exists.")
    SOURCE_ROOT.rename(target)
    print(f"Legacy data moved to {target}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build canonical dataset from legacy aggregations.")
    parser.add_argument("--output", default="canonical_dataset", help="Directory for canonical artefacts.")
    parser.add_argument(
        "--archive-dir",
        default="archive",
        help="Directory where the legacy aggregated_reflectance directory will be moved.",
    )
    parser.add_argument(
        "--archive",
        action="store_true",
        help="Move the original aggregated_reflectance/ directory into --archive-dir after building.",
    )
    args = parser.parse_args()

    output_dir = Path(args.output).resolve()
    build_canonical_dataset(output_dir)
    print(f"Canonical dataset written to {output_dir}")

    if args.archive:
        archive_root = Path(args.archive_dir).resolve()
        archive_legacy(output_dir, archive_root)
        print(f"Archive complete: {archive_root}")


if __name__ == "__main__":
    main()
