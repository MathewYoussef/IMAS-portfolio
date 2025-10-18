#!/usr/bin/env python3
"""
Fit advanced reflectance baselines using Phase 9 spectral descriptors.

The script trains linear, ridge, and lasso models for each modality/kind/
band combination, computes leave-one-dose-out diagnostics, and records
predictions plus provenance for downstream consolidation.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


DEFAULT_MANIFEST = Path("analysis/phase9/inputs/phase9_input_manifest.json")

TARGET_SPECS: Dict[str, Dict[str, str]] = {
    "chrom_total": {
        "column": "chrom_total_mg_per_gDW_trimmed_mean",
        "modality": "chrom",
        "component": "total",
    },
    "chrom_oxidized": {
        "column": "chrom_oxidized_mg_per_gDW_trimmed_mean",
        "modality": "chrom",
        "component": "oxidized",
    },
    "chrom_reduced": {
        "column": "chrom_reduced_mg_per_gDW_trimmed_mean",
        "modality": "chrom",
        "component": "reduced",
    },
    "dad_total": {
        "column": "dad_total_mg_per_gDW_trimmed_mean",
        "modality": "dad",
        "component": "total",
    },
    "dad_oxidized": {
        "column": "dad_oxidized_mg_per_gDW_trimmed_mean",
        "modality": "dad",
        "component": "oxidized",
    },
    "dad_reduced": {
        "column": "dad_reduced_mg_per_gDW_trimmed_mean",
        "modality": "dad",
        "component": "reduced",
    },
    "latent_total": {
        "column": "total_latent_mean_trimmed",
        "modality": "latent",
        "component": "total",
    },
}

FEATURE_COLUMNS = [
    "peak_reflectance",
    "peak_wavelength_nm",
    "auc",
    "continuum_depth",
    "fwhm_nm",
    "max_curvature",
    "fixed_peak_reflectance",
    "fixed_continuum_depth",
]


@dataclass
class ModelSpec:
    name: str
    estimator_type: str
    alpha: float | None = None


@dataclass
class AdvancedBaselineConfig:
    manifest_path: Path
    features_path: Path
    dose_summary_path: Path
    output_summary_path: Path
    output_predictions_path: Path
    plots_dir: Path
    models: List[ModelSpec]
    random_seed: int


def load_manifest(path: Path) -> Dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run advanced baseline models for Phase 9 features.")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="Path to the Phase 9 manifest.")
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("analysis/phase9/tables/advanced_band_features.csv"),
        help="Advanced feature matrix.",
    )
    parser.add_argument(
        "--dose-summary",
        type=Path,
        default=Path("canonical_dataset/dose_level_canonical_summary.csv"),
        help="Canonical dose-level summary with concentration targets.",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path("analysis/phase9/tables/advanced_baseline_summary.csv"),
        help="Destination CSV for model diagnostics.",
    )
    parser.add_argument(
        "--predictions-output",
        type=Path,
        default=Path("analysis/phase9/tables/advanced_baseline_predictions.csv"),
        help="Destination CSV for per-dose predictions.",
    )
    return parser.parse_args()


def fit_linear_model(X: np.ndarray, y: np.ndarray, alpha: float = 0.0) -> Tuple[float, np.ndarray]:
    n_samples, n_features = X.shape
    ones = np.ones((n_samples, 1), dtype=float)
    design = np.hstack([ones, X])
    regularisation = alpha * np.eye(n_features + 1, dtype=float)
    regularisation[0, 0] = 0.0  # Do not penalise intercept
    gram = design.T @ design + regularisation
    target_vec = design.T @ y
    try:
        coefs = np.linalg.solve(gram, target_vec)
    except np.linalg.LinAlgError:
        coefs = np.linalg.pinv(gram) @ target_vec
    intercept = float(coefs[0])
    weights = np.asarray(coefs[1:], dtype=float)
    return intercept, weights


def predict_with_model(intercept: float, weights: np.ndarray, X: np.ndarray) -> np.ndarray:
    return X @ weights + intercept


def loocv_predictions_manual(X: np.ndarray, y: np.ndarray, alpha: float) -> np.ndarray:
    n = X.shape[0]
    preds = np.zeros(n, dtype=float)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        intercept, weights = fit_linear_model(X[mask], y[mask], alpha)
        preds[i] = predict_with_model(intercept, weights, X[[i]])[0]
    return preds


def compute_predictions(
    X: np.ndarray,
    y: np.ndarray,
    model_spec: ModelSpec,
    random_seed: int,
) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    estimator = model_spec.estimator_type
    if estimator == "ridge":
        alpha = model_spec.alpha if model_spec.alpha is not None else 1.0
        intercept, weights = fit_linear_model(X, y, alpha)
        predictions = predict_with_model(intercept, weights, X)
        loo_predictions = loocv_predictions_manual(X, y, alpha)
    elif estimator == "linear":
        alpha = 0.0
        intercept, weights = fit_linear_model(X, y, alpha)
        predictions = predict_with_model(intercept, weights, X)
        loo_predictions = loocv_predictions_manual(X, y, alpha)
    elif estimator == "lasso":
        alpha = model_spec.alpha if model_spec.alpha is not None else 0.1
        trained = Lasso(alpha=alpha, fit_intercept=True, max_iter=10000, random_state=random_seed)
        trained.fit(X, y)
        predictions = trained.predict(X)
        intercept = float(trained.intercept_)
        weights = trained.coef_.astype(float)

        loo_predictions = np.zeros(len(y), dtype=float)
        for i in range(len(y)):
            mask = np.ones(len(y), dtype=bool)
            mask[i] = False
            fold_model = Lasso(alpha=alpha, fit_intercept=True, max_iter=10000, random_state=random_seed)
            fold_model.fit(X[mask], y[mask])
            loo_predictions[i] = fold_model.predict(X[[i]])[0]
    else:
        raise ValueError(f"Unsupported estimator type: {estimator}")

    return intercept, weights, predictions, loo_predictions


def plot_predictions(
    output_dir: Path,
    target: str,
    kind: str,
    band_label: str,
    model_name: str,
    actual: np.ndarray,
    predicted: np.ndarray,
) -> str:
    output_dir.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(5, 5))
    axis.scatter(actual, predicted, color="#2ca02c")
    min_val = min(actual.min(), predicted.min())
    max_val = max(actual.max(), predicted.max())
    axis.plot([min_val, max_val], [min_val, max_val], linestyle="--", color="black")
    for dose, a, p in zip(range(1, len(actual) + 1), actual, predicted):
        axis.annotate(str(dose), (a, p), fontsize=8)
    axis.set_xlabel("Actual")
    axis.set_ylabel("Predicted")
    axis.set_title(f"{model_name} – {target} ({kind}, {band_label})")
    axis.grid(alpha=0.3)
    figure.tight_layout()
    filename = output_dir / f"{model_name}_{target}_{kind}_{band_label}.png"
    figure.savefig(filename, dpi=300)
    plt.close(figure)
    return filename.name


def collect_models(config: Dict[str, Any]) -> Tuple[List[ModelSpec], int]:
    modeling = config.get("modeling", {})
    model_entries = modeling.get("models", [])
    if not model_entries:
        model_entries = [
            {"name": "linear", "type": "linear"},
            {"name": "ridge", "type": "ridge", "alpha": 1.0},
        ]
    models = [ModelSpec(name=entry["name"], estimator_type=entry["type"], alpha=entry.get("alpha")) for entry in model_entries]
    random_seed = int(modeling.get("random_seed", 42))
    return models, random_seed


def run_advanced_baseline(config: AdvancedBaselineConfig) -> None:
    manifest = load_manifest(config.manifest_path)
    features_df = pd.read_csv(config.features_path)
    dose_summary = pd.read_csv(config.dose_summary_path)

    models = config.models

    target_columns = {key: spec["column"] for key, spec in TARGET_SPECS.items()}

    merge_columns = ["dose_id"] + list(target_columns.values())
    merge_targets = dose_summary[merge_columns].copy()

    summary_records: List[Dict[str, Any]] = []
    prediction_records: List[Dict[str, Any]] = []
    plots: List[str] = []

    for (kind, band_label), subset in features_df.groupby(["kind", "band_label"]):
        subset = subset.merge(merge_targets, on="dose_id", how="inner")
        if subset.empty:
            continue

        X = subset[FEATURE_COLUMNS].fillna(0.0).to_numpy(dtype=float)
        metadata = subset.iloc[0]
        band_lower = float(subset["actual_band_lower_nm"].iloc[0])
        band_upper = float(subset["actual_band_upper_nm"].iloc[0])

        for target_name, target_spec in TARGET_SPECS.items():
            target_column = target_spec["column"]
            if target_column not in subset:
                continue

            y = subset[target_column].to_numpy(dtype=float)
            dose_ids = subset["dose_id"].tolist()
            component = target_spec["component"]
            modality = target_spec["modality"]

            for model_spec in models:
                intercept, weights, predictions, loo_predictions = compute_predictions(
                    X=X,
                    y=y,
                    model_spec=model_spec,
                    random_seed=config.random_seed,
                )

                rmse = float(np.sqrt(mean_squared_error(y, predictions)))
                mae = float(mean_absolute_error(y, predictions))
                r2 = float(r2_score(y, predictions))

                rmse_loo = float(np.sqrt(mean_squared_error(y, loo_predictions)))
                mae_loo = float(mean_absolute_error(y, loo_predictions))

                coefficients = weights.tolist()

                summary_records.append(
                    {
                        "target": target_name,
                        "modality": modality,
                        "kind": kind,
                        "band_label": band_label,
                        "component": component,
                        "band_lower_nm": band_lower,
                        "band_upper_nm": band_upper,
                        "model_type": model_spec.name,
                        "estimator": model_spec.estimator_type,
                        "n_samples": int(len(subset)),
                        "rmse": rmse,
                        "mae": mae,
                        "r2": r2,
                        "rmse_loocv": rmse_loo,
                        "mae_loocv": mae_loo,
                        "coefficients": coefficients,
                        "intercept": intercept,
                        "feature_columns": FEATURE_COLUMNS,
                    }
                )

                for dose_id, actual, pred, pred_loo in zip(dose_ids, y, predictions, loo_predictions):
                    prediction_records.append(
                        {
                            "dose_id": dose_id,
                            "target": target_name,
                            "modality": modality,
                            "kind": kind,
                            "band_label": band_label,
                            "model_type": model_spec.name,
                            "component": component,
                            "actual": float(actual),
                            "predicted": float(pred),
                            "predicted_loocv": float(pred_loo),
                        }
                    )

                if model_spec.name == "ridge":
                    plot_name = plot_predictions(
                        output_dir=config.plots_dir,
                        target=target_name,
                        kind=kind,
                        band_label=band_label,
                        model_name=model_spec.name,
                        actual=y,
                        predicted=predictions,
                    )
                    plots.append(plot_name)

    summary_df = pd.DataFrame(summary_records)
    predictions_df = pd.DataFrame(prediction_records)

    config.output_summary_path.parent.mkdir(parents=True, exist_ok=True)
    config.output_predictions_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(config.output_summary_path, index=False)
    predictions_df.to_csv(config.output_predictions_path, index=False)

    output_manifest = {
        "source_manifest": str(config.manifest_path),
        "advanced_features": str(config.features_path),
        "dose_summary": str(config.dose_summary_path),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "models": [model.__dict__ for model in models],
        "summary_table": str(config.output_summary_path),
        "predictions_table": str(config.output_predictions_path),
        "plots_dir": str(config.plots_dir),
        "plots": plots,
        "targets": list(TARGET_SPECS.keys()),
        "components": sorted({spec["component"] for spec in TARGET_SPECS.values()}),
        "feature_columns": FEATURE_COLUMNS,
        "random_seed": config.random_seed,
        "notes": (
            "Advanced baselines fit per kind/band with leave-one-dose-out diagnostics. "
            "Linear and LASSO variants remain for QA; downstream reports highlight the ridge runs."
        ),
    }

    manifest_path = config.output_summary_path.with_name("advanced_baseline_manifest.json")
    manifest_path.write_text(json.dumps(output_manifest, indent=2))

    print(f"[INFO] Advanced baseline summary written to {config.output_summary_path}")
    print(f"[INFO] Advanced baseline predictions written to {config.output_predictions_path}")
    print(f"[INFO] Advanced baseline manifest written to {manifest_path}")


def main() -> None:
    args = parse_args()
    manifest = load_manifest(args.manifest)
    models, random_seed = collect_models(manifest)

    config = AdvancedBaselineConfig(
        manifest_path=args.manifest,
        features_path=args.features,
        dose_summary_path=args.dose_summary,
        output_summary_path=args.summary_output,
        output_predictions_path=args.predictions_output,
        plots_dir=Path(manifest["output_paths"]["plots_dir"]).resolve(),
        models=models,
        random_seed=random_seed,
    )
    run_advanced_baseline(config)


if __name__ == "__main__":
    main()
