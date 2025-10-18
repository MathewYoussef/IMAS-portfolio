#!/usr/bin/env python3
"""Run predefined feature-subset LOOCV sweeps for Phase 9 ridge equations."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable, List, Sequence

import numpy as np
import pandas as pd


DEFAULT_INVENTORY = Path("analysis/phase9/equations/model_inventory.json")
DEFAULT_FEATURES = Path("analysis/phase9/tables/advanced_band_features.csv")
DEFAULT_TARGETS_FILE = Path("canonical_dataset/dose_level_canonical_summary.csv")
DEFAULT_OUTPUT_SUMMARY = Path("analysis/phase9/equations/outputs/feature_subset_sweep_summary.csv")
DEFAULT_OUTPUT_FOLDS = Path("analysis/phase9/equations/outputs/feature_subset_sweep_folds.csv")

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


def default_subsets(spec: ModelSpec, include_full: bool) -> List[tuple[str, List[str]]]:
    # Core features we frequently ablate around
    preferred = [
        ("auc_only", ["auc"]),
        ("fwhm_only", ["fwhm_nm"]),
        ("auc_fwhm", ["auc", "fwhm_nm"]),
        ("peak_auc", ["peak_reflectance", "auc"]),
        ("peak_only", ["peak_reflectance"]),
    ]

    subsets: List[tuple[str, List[str]]] = []
    available = set(spec.feature_columns)

    # Always include a "fixed-term drop" variant when fixed peaks are present.
    no_fixed = [f for f in spec.feature_columns if not f.startswith("fixed_")]
    if len(no_fixed) != len(spec.feature_columns):
        subsets.append(("drop_fixed_terms", no_fixed))

    for name, cols in preferred:
        cols_present = [c for c in cols if c in available]
        if len(cols_present) == len(cols) and cols_present:
            subsets.append((name, cols_present))

    # Single-feature sweeps for remaining descriptors (beyond the preferred list)
    for col in spec.feature_columns:
        if col.startswith("fixed_"):
            continue
        if any(col in cols for _, cols in subsets):
            continue
        subsets.append((f"{col}_only", [col]))

    # Minimal two-feature combos (to gauge interactions) – limit to first 5 features to keep noise manageable.
    base_cols = [c for c in spec.feature_columns if not c.startswith("fixed_")][:5]
    for combo in combinations(base_cols, 2):
        subsets.append(("+".join(combo), list(combo)))

    if include_full:
        subsets.append(("full", list(spec.feature_columns)))

    # Ensure uniqueness while preserving order.
    seen = set()
    unique_subsets: List[tuple[str, List[str]]] = []
    for name, cols in subsets:
        key = tuple(cols)
        if key in seen:
            continue
        seen.add(key)
        unique_subsets.append((name, cols))
    return unique_subsets


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


def ridge_fit(X: np.ndarray, y: np.ndarray, alpha: float) -> tuple[np.ndarray, float]:
    if X.ndim != 2:
        raise ValueError("Expected 2D design matrix")
    n_samples, n_features = X.shape
    if n_samples == 0:
        raise ValueError("Empty training set")

    X_aug = np.hstack([X, np.ones((n_samples, 1), dtype=float)])
    penalty = np.eye(n_features + 1, dtype=float)
    penalty[-1, -1] = 0.0  # do not regularise intercept
    try:
        params = np.linalg.solve(X_aug.T @ X_aug + alpha * penalty, X_aug.T @ y)
    except np.linalg.LinAlgError:
        params, *_ = np.linalg.lstsq(X_aug.T @ X_aug + alpha * penalty, X_aug.T @ y, rcond=None)
    coefficients = params[:-1]
    intercept = params[-1]
    return coefficients, float(intercept)


def loocv_metrics(
    spec: ModelSpec,
    dataset: pd.DataFrame,
    feature_names: Sequence[str],
    alpha: float,
) -> tuple[dict[str, float], list[dict[str, object]]]:
    feature_names = list(feature_names)
    missing = [col for col in feature_names if col not in dataset.columns]
    if missing:
        raise ValueError(f"{spec.target}: missing columns {missing}")

    fold_rows: List[dict[str, object]] = []
    preds: List[float] = []
    actuals: List[float] = []
    dose_ids = dataset["dose_id"].unique()
    for dose_id in dose_ids:
        train = dataset[dataset["dose_id"] != dose_id]
        test = dataset[dataset["dose_id"] == dose_id]

        X_train = train[feature_names].fillna(0.0).to_numpy(dtype=float)
        y_train = train["target_value"].to_numpy(dtype=float)

        X_test = test[feature_names].fillna(0.0).to_numpy(dtype=float)
        y_test = test["target_value"].to_numpy(dtype=float)

        coeffs, intercept = ridge_fit(X_train, y_train, alpha)
        y_pred = X_test @ coeffs + intercept

        preds.extend(y_pred.tolist())
        actuals.extend(y_test.tolist())

        for d_id, pred_val, actual_val in zip(test["dose_id"], y_pred, y_test):
            fold_rows.append(
                {
                    "target": spec.target,
                    "kind": spec.kind,
                    "band_label": spec.band_label,
                    "dose_id": d_id,
                    "predicted": float(pred_val),
                    "actual": float(actual_val),
                }
            )

    residuals = np.array(preds) - np.array(actuals)
    rmse = float(math.sqrt(np.mean(residuals ** 2)))
    mae = float(np.mean(np.abs(residuals)))
    return {"rmse_loocv": rmse, "mae_loocv": mae}, fold_rows


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Feature-subset LOOCV sweeps for Phase 9 equations")
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--targets-file", type=Path, default=DEFAULT_TARGETS_FILE, help="Canonical concentrations CSV")
    parser.add_argument("--role", action="append", dest="roles", help="Filter models by reporting role")
    parser.add_argument("--target", action="append", dest="filter_targets", help="Filter models by target id")
    parser.add_argument("--alpha", type=float, default=0.1, help="Ridge regularisation strength")
    parser.add_argument("--include-full", action="store_true", help="Include full feature set in sweep output")
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_OUTPUT_SUMMARY)
    parser.add_argument("--fold-output", type=Path, default=DEFAULT_OUTPUT_FOLDS)
    args = parser.parse_args(argv)

    models = load_inventory(args.inventory)
    selected = select_models(models, args.filter_targets, args.roles)
    if not selected:
        raise SystemExit("No models selected for sweep")

    features_df = pd.read_csv(args.features)
    targets_df = pd.read_csv(args.targets_file)

    summary_records: List[dict[str, object]] = []
    fold_records: List[dict[str, object]] = []

    for spec in selected:
        dataset = prepare_dataset(spec, features_df, targets_df)
        subsets = default_subsets(spec, include_full=args.include_full)
        for name, columns in subsets:
            metrics, fold_rows = loocv_metrics(spec, dataset, columns, args.alpha)
            summary_records.append(
                {
                    "target": spec.target,
                    "kind": spec.kind,
                    "band_label": spec.band_label,
                    "role": spec.reporting_role,
                    "subset_name": name,
                    "feature_count": len(columns),
                    "features": ",".join(columns),
                    "alpha": args.alpha,
                    "rmse_loocv_subset": metrics["rmse_loocv"],
                    "mae_loocv_subset": metrics["mae_loocv"],
                    "rmse_loocv_catalogue": spec.rmse_loocv,
                    "rmse_in_sample_catalogue": spec.rmse,
                    "delta_rmse_loocv": metrics["rmse_loocv"] - spec.rmse_loocv,
                }
            )
            for fold in fold_rows:
                fold["subset_label"] = name
                fold_records.append(fold)

    summary_df = pd.DataFrame(summary_records)
    summary_df.sort_values(["target", "subset_name"], inplace=True)

    fold_df = pd.DataFrame(fold_records)
    if not fold_df.empty:
        fold_df.sort_values(["target", "subset_label", "dose_id"], inplace=True)

    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(args.summary_output, index=False)

    args.fold_output.parent.mkdir(parents=True, exist_ok=True)
    fold_df.to_csv(args.fold_output, index=False)

    print(summary_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
