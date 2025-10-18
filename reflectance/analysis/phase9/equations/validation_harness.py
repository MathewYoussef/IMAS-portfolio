#!/usr/bin/env python3
"""Held-out validation harness for Phase 9 ridge equations."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

import numpy as np
import pandas as pd


DEFAULT_INVENTORY = Path("analysis/phase9/equations/model_inventory.json")
DEFAULT_FEATURES = Path("analysis/phase9/tables/advanced_band_features.csv")
DEFAULT_TARGETS_FILE = Path("canonical_dataset/dose_level_canonical_summary.csv")
DEFAULT_OUTPUT_SUMMARY = Path("analysis/phase9/equations/outputs/validation_summary.csv")
DEFAULT_OUTPUT_FOLDS = Path("analysis/phase9/equations/outputs/validation_folds.csv")

TARGET_TO_COLUMN = {
    "chrom_total": "chrom_total_mg_per_gDW_trimmed_mean",
    "chrom_oxidized": "chrom_oxidized_mg_per_gDW_trimmed_mean",
    "chrom_reduced": "chrom_reduced_mg_per_gDW_trimmed_mean",
    "dad_total": "dad_total_mg_per_gDW_trimmed_mean",
    "dad_oxidized": "dad_oxidized_mg_per_gDW_trimmed_mean",
    "dad_reduced": "dad_reduced_mg_per_gDW_trimmed_mean",
    "latent_total": "total_latent_mean_trimmed",
}


@dataclass
class ModelSpec:
    target: str
    modality: str
    component: str
    kind: str
    band_label: str
    model_type: str
    estimator: str
    reporting_role: str
    rmse_loocv: float
    mae_loocv: float
    rmse: float
    mae: float
    r2: float
    coefficients: np.ndarray
    intercept: float
    feature_columns: List[str]


def load_inventory(path: Path) -> List[ModelSpec]:
    payload = json.loads(path.read_text())
    models: List[ModelSpec] = []
    for entry in payload.get("models", []):
        models.append(
            ModelSpec(
                target=entry["target"],
                modality=entry["modality"],
                component=entry["component"],
                kind=entry["kind"],
                band_label=entry["band_label"],
                model_type=entry["model_type"],
                estimator=entry["estimator"],
                reporting_role=entry.get("reporting_role", "unknown"),
                rmse_loocv=float(entry["rmse_loocv"]),
                mae_loocv=float(entry["mae_loocv"]),
                rmse=float(entry["rmse"]),
                mae=float(entry["mae"]),
                r2=float(entry["r2"]),
                coefficients=np.array(entry["coefficients"], dtype=float),
                intercept=float(entry["intercept"]),
                feature_columns=list(entry["feature_columns"]),
            )
        )
    return models


def select_models(models: Iterable[ModelSpec], targets: Sequence[str] | None, roles: Sequence[str] | None) -> List[ModelSpec]:
    selected: List[ModelSpec] = []
    for model in models:
        if targets and model.target not in targets:
            continue
        if roles and model.reporting_role not in roles:
            continue
        selected.append(model)
    return selected


def prepare_dataset(spec: ModelSpec, features_df: pd.DataFrame, targets_df: pd.DataFrame) -> pd.DataFrame:
    subset = features_df[(features_df["kind"] == spec.kind) & (features_df["band_label"] == spec.band_label)].copy()
    if subset.empty:
        raise ValueError(f"No feature rows for {spec.target} ({spec.kind}, {spec.band_label})")

    target_column = TARGET_TO_COLUMN.get(spec.target)
    if target_column is None:
        raise ValueError(f"Target mapping missing for {spec.target}")

    merged = subset.merge(targets_df[["dose_id", target_column]], on="dose_id", how="left")
    if merged[target_column].isna().any():
        missing = merged[merged[target_column].isna()]["dose_id"].unique()
        raise ValueError(f"Missing concentrations for {spec.target}; dose_ids={missing}")

    merged = merged.rename(columns={target_column: "target_value"})
    return merged


def evaluate_spec(
    spec: ModelSpec,
    dataset: pd.DataFrame,
) -> tuple[dict[str, object], List[dict[str, object]]]:
    X = dataset[spec.feature_columns].fillna(0.0).to_numpy(dtype=float)
    y = dataset["target_value"].to_numpy(dtype=float)

    predictions = X @ spec.coefficients + spec.intercept
    residuals = predictions - y

    rmse = float(math.sqrt(np.mean(residuals ** 2)))
    mae = float(np.mean(np.abs(residuals)))

    summary = {
        "target": spec.target,
        "kind": spec.kind,
        "band_label": spec.band_label,
        "role": spec.reporting_role,
        "rmse_in_sample_replay": rmse,
        "mae_in_sample_replay": mae,
        "rmse_catalogue": spec.rmse,
        "mae_catalogue": spec.mae,
        "rmse_loocv_catalogue": spec.rmse_loocv,
        "mae_loocv_catalogue": spec.mae_loocv,
        "max_abs_residual": float(np.max(np.abs(residuals))),
    }

    fold_rows: List[dict[str, object]] = []
    for dose_id, pred, actual in zip(dataset["dose_id"], predictions, y):
        fold_rows.append(
            {
                "target": spec.target,
                "kind": spec.kind,
                "band_label": spec.band_label,
                "dose_id": dose_id,
                "predicted": float(pred),
                "actual": float(actual),
                "residual": float(pred - actual),
            }
        )

    return summary, fold_rows


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validation harness for Phase 9 ridge equations")
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--targets-file", type=Path, default=DEFAULT_TARGETS_FILE, help="Canonical concentration table")
    parser.add_argument("--role", action="append", dest="roles", help="Filter models by reporting role")
    parser.add_argument("--target", action="append", dest="filter_targets", help="Filter models by target id")
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_OUTPUT_SUMMARY)
    parser.add_argument("--fold-output", type=Path, default=DEFAULT_OUTPUT_FOLDS)
    args = parser.parse_args(argv)

    models = load_inventory(args.inventory)
    selected = select_models(models, args.filter_targets, args.roles)
    if not selected:
        raise SystemExit("No models selected for validation")

    features_df = pd.read_csv(args.features)
    targets_df = pd.read_csv(args.targets_file)

    summary_records: List[dict[str, object]] = []
    fold_records: List[dict[str, object]] = []

    for spec in selected:
        dataset = prepare_dataset(spec, features_df, targets_df)
        summary, folds = evaluate_spec(spec, dataset)
        summary_records.append(summary)
        fold_records.extend(folds)

    summary_df = pd.DataFrame(summary_records)
    summary_df.sort_values(["target"], inplace=True)
    fold_df = pd.DataFrame(fold_records)
    fold_df.sort_values(["target", "dose_id"], inplace=True)

    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(args.summary_output, index=False)

    args.fold_output.parent.mkdir(parents=True, exist_ok=True)
    fold_df.to_csv(args.fold_output, index=False)

    print(summary_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
