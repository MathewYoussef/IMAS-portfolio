#!/usr/bin/env python3
"""Re-evaluate Phase 9 ridge equations and compare against stored outputs."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

import numpy as np
import pandas as pd


DEFAULT_INVENTORY = Path("analysis/phase9/equations/model_inventory.json")
DEFAULT_FEATURES = Path("analysis/phase9/tables/advanced_band_features.csv")
DEFAULT_TARGETS_FILE = Path("canonical_dataset/dose_level_canonical_summary.csv")
DEFAULT_PREDICTIONS = Path("analysis/phase9/tables/advanced_baseline_predictions.csv")

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


def evaluate_model(
    spec: ModelSpec,
    features_df: pd.DataFrame,
    targets_df: pd.DataFrame,
    predictions_df: pd.DataFrame,
    tolerance: float,
) -> dict:
    subset = features_df[(features_df["kind"] == spec.kind) & (features_df["band_label"] == spec.band_label)].copy()
    if subset.empty:
        raise ValueError(f"No feature rows for {spec.target} ({spec.kind}, {spec.band_label})")

    missing = [col for col in spec.feature_columns if col not in subset.columns]
    if missing:
        raise ValueError(f"Missing feature columns for {spec.target}: {missing}")

    X = subset[spec.feature_columns].fillna(0.0).to_numpy(dtype=float)
    preds = X @ spec.coefficients + spec.intercept

    target_column = TARGET_TO_COLUMN.get(spec.target)
    if target_column is None:
        raise ValueError(f"Target mapping missing for {spec.target}")
    merged = subset[["dose_id"]].merge(targets_df[["dose_id", target_column]], on="dose_id", how="left")
    if merged[target_column].isna().any():
        raise ValueError(f"Missing concentration values for {spec.target}")
    actual = merged[target_column].to_numpy(dtype=float)

    residuals = preds - actual
    rmse = float(np.sqrt(np.mean(residuals ** 2)))
    mae = float(np.mean(np.abs(residuals)))

    diffs: dict[str, float | bool | str] = {}
    if not predictions_df.empty:
        mask = (
            (predictions_df["target"] == spec.target)
            & (predictions_df["kind"] == spec.kind)
            & (predictions_df["band_label"] == spec.band_label)
            & (predictions_df["model_type"] == spec.model_type)
        )
        baseline = predictions_df[mask]
        if not baseline.empty:
            baseline = baseline.sort_values("dose_id")
            predicted = baseline["predicted"].to_numpy(dtype=float)
            if predicted.shape == preds.shape:
                delta = np.abs(predicted - preds)
                diffs["max_pred_delta"] = float(delta.max())
                diffs["within_tolerance"] = bool((delta <= tolerance + 1e-15).all())
            diffs["loocv_available"] = True
        else:
            diffs["note"] = "No baseline predictions found"

    return {
        "target": spec.target,
        "kind": spec.kind,
        "band_label": spec.band_label,
        "role": spec.reporting_role,
        "rmse_recomputed": rmse,
        "mae_recomputed": mae,
        "rmse_catalogue": spec.rmse,
        "rmse_loocv_catalogue": spec.rmse_loocv,
        **diffs,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Evaluate Phase 9 ridge equations against canonical data")
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--targets-file", type=Path, default=DEFAULT_TARGETS_FILE)
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--role", action="append", dest="roles")
    parser.add_argument("--target", action="append", dest="filter_targets")
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)

    models = load_inventory(args.inventory)
    selected = select_models(models, getattr(args, "filter_targets", None), args.roles)
    if not selected:
        print("[WARN] No models selected", file=sys.stderr)
        return 1

    features_df = pd.read_csv(args.features)
    targets_df = pd.read_csv(args.targets_file)
    predictions_df = pd.read_csv(args.predictions) if args.predictions.exists() else pd.DataFrame()

    records = [
        evaluate_model(spec, features_df, targets_df, predictions_df, args.tolerance)
        for spec in selected
    ]
    summary = pd.DataFrame(records)
    print(summary.to_string(index=False))

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.output, index=False)

    if "within_tolerance" in summary.columns and not summary["within_tolerance"].fillna(True).all():
        print("[WARN] Deviations exceeded tolerance", file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
