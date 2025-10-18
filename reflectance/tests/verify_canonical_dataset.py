#!/usr/bin/env python3
"""Regression check that canonical outputs still mirror the archived source.

The script compares every CSV/JSON artefact in ``canonical_dataset/`` against
its counterpart under ``archive/aggregated_reflectance/``.  Any divergence
greater than a tiny floating-point tolerance triggers a non-zero exit so the
CI/build can halt before publishing inconsistent data.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Iterable, Sequence


DEFAULT_ARCHIVE = Path("archive") / "aggregated_reflectance"
DEFAULT_CANONICAL = Path("canonical_dataset")
FLOAT_TOLERANCE = 1e-10


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _sorted_rows(rows: Iterable[dict[str, str]], keys: Sequence[str]) -> list[dict[str, str]]:
    return sorted(rows, key=lambda row: tuple(row[key] for key in keys))


def _is_float(value: str | None) -> bool:
    if value is None:
        return False
    stripped = value.strip()
    if stripped == "":
        return False
    try:
        float(stripped)
    except ValueError:
        return False
    return True


def _equal(value_a: str, value_b: str) -> bool:
    if _is_float(value_a) and _is_float(value_b):
        a = float(value_a)
        b = float(value_b)
        if math.isnan(a) and math.isnan(b):
            return True
        return abs(a - b) <= FLOAT_TOLERANCE
    return value_a == value_b


def _assert_same_rows(
    archive_rows: list[dict[str, str]],
    canonical_rows: list[dict[str, str]],
    sort_keys: Sequence[str],
    label: str,
    allowed_extra: set[str] | None = None,
) -> None:
    if len(archive_rows) != len(canonical_rows):
        raise AssertionError(f"{label}: row count differs ({len(archive_rows)} vs {len(canonical_rows)})")

    archive_sorted = _sorted_rows(archive_rows, sort_keys)
    canonical_sorted = _sorted_rows(canonical_rows, sort_keys)

    archive_columns = set(archive_rows[0].keys()) if archive_rows else set()
    canonical_columns = set(canonical_rows[0].keys()) if canonical_rows else set()
    allowed_extra = allowed_extra or set()
    extra_columns = canonical_columns - archive_columns
    unexpected = extra_columns - allowed_extra
    if unexpected:
        raise AssertionError(f"{label}: unexpected canonical columns {sorted(unexpected)}")

    compare_columns = sorted(archive_columns)

    for idx, (archive_row, canonical_row) in enumerate(zip(archive_sorted, canonical_sorted), start=1):
        for column in compare_columns:
            if not _equal(archive_row[column], canonical_row[column]):
                raise AssertionError(
                    f"{label}: difference at row {idx}, column '{column}' "
                    f"(archive={archive_row[column]!r}, canonical={canonical_row[column]!r})"
                )


def check_dose_summary(archive_root: Path, canonical_root: Path) -> None:
    archive_rows = _read_csv(archive_root / "dose_summary.csv")
    canonical_rows = _read_csv(canonical_root / "dose_summary.csv")
    canonical_lookup = {
        (float(row["uva_mw_cm2"]), float(row["uvb_mw_cm2"])): row for row in canonical_rows
    }

    if len(archive_rows) != len(canonical_rows):
        raise AssertionError("dose_summary.csv: differing row counts")

    for archive_row in archive_rows:
        key = (float(archive_row["uva"]), float(archive_row["uvb"]))
        canonical_row = canonical_lookup.get(key)
        if canonical_row is None:
            raise AssertionError(f"dose_summary.csv: UAV/UVB pair {key} missing in canonical dataset")

        for column, value in archive_row.items():
            if column in {"dose_id", "uva", "uvb"}:
                continue
            canonical_value = canonical_row[column]
            if not _equal(value, canonical_value):
                raise AssertionError(
                    "dose_summary.csv: mismatch for UVA/UVB "
                    f"{key} column '{column}' (archive={value!r}, canonical={canonical_value!r})"
                )


def check_json(archive_root: Path, canonical_root: Path) -> None:
    archive_json = json.loads((archive_root / "dose_reflectance_dad_correlations.json").read_text())
    canonical_json = json.loads((canonical_root / "dose_reflectance_dad_correlations.json").read_text())
    if archive_json != canonical_json:
        raise AssertionError("dose_reflectance_dad_correlations.json: JSON content differs")


def compare_all(archive_root: Path, canonical_root: Path) -> None:
    required_files = [
        "reflectance_trimmed_stats.csv",
        "dose_reflectance_stats.csv",
        "dose_reflectance_composites.csv",
        "dose_reflectance_features_angles.csv",
        "dose_reflectance_features_composites.csv",
        "dose_dad_concentrations.csv",
        "precision_weighted_concentrations.csv",
        "precision_weighted_concentrations_treatment.csv",
        "dose_reflectance_dad_summary.csv",
    ]

    for path in required_files:
        archive_path = archive_root / path
        canonical_path = canonical_root / path
        if not archive_path.exists():
            raise FileNotFoundError(f"Missing archive file: {archive_path}")
        if not canonical_path.exists():
            raise FileNotFoundError(f"Missing canonical file: {canonical_path}")

    check_dose_summary(archive_root, canonical_root)
    check_json(archive_root, canonical_root)

    compare_matrix = [
        ("reflectance_trimmed_stats.csv", ("dose_id", "source_dir", "sample_label", "angle")),
        ("dose_reflectance_stats.csv", ("dose_id", "angle")),
        ("dose_reflectance_composites.csv", ("dose_id", "composite")),
        ("dose_reflectance_features_angles.csv", ("dose_id", "kind")),
        ("dose_reflectance_features_composites.csv", ("dose_id", "kind")),
        ("dose_dad_concentrations.csv", ("dose_id",)),
        ("precision_weighted_concentrations.csv", ("sample_id",)),
        ("precision_weighted_concentrations_treatment.csv", ("dose_id",)),
    ]

    extra_columns = {"uva_mw_cm2", "uvb_mw_cm2"}

    for filename, keys in compare_matrix:
        archive_rows = _read_csv(archive_root / filename)
        canonical_rows = _read_csv(canonical_root / filename)
        _assert_same_rows(archive_rows, canonical_rows, keys, filename, allowed_extra=extra_columns)

    check_reflectance_dad_summary(archive_root, canonical_root)


def check_reflectance_dad_summary(archive_root: Path, canonical_root: Path) -> None:
    archive_rows = _read_csv(archive_root / "dose_reflectance_dad_summary.csv")
    canonical_rows = _read_csv(canonical_root / "dose_reflectance_dad_summary.csv")

    for row in archive_rows:
        if row["dose_id_reflectance"] != row["dose_id_dad"]:
            raise AssertionError(
                "dose_reflectance_dad_summary.csv: archive mismatch between reflectance and DAD dose IDs"
            )
        row["dose_id"] = row["dose_id_reflectance"]
        del row["dose_id_reflectance"]
        del row["dose_id_dad"]

    compare_columns = ("dose_id", "kind")
    _assert_same_rows(
        archive_rows,
        canonical_rows,
        compare_columns,
        "dose_reflectance_dad_summary.csv",
        allowed_extra={"uva_mw_cm2", "uvb_mw_cm2"},
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify canonical_dataset/ matches archive/aggregated_reflectance outputs."
    )
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE, help="Path to the archived source data")
    parser.add_argument("--canonical", type=Path, default=DEFAULT_CANONICAL, help="Path to the canonical dataset")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    compare_all(args.archive, args.canonical)
    print("Canonical dataset matches archive source within tolerance.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # pragma: no cover - simple CLI guard
        print(f"[verify_canonical_dataset] ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
