#!/usr/bin/env python3
"""
Fit simple baseline regressions between band features and concentration targets.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


TARGET_COLUMNS = {
    "chrom_total": "chrom_total_mg_per_gDW_trimmed_mean",
    "dad_total": "dad_total_mg_per_gDW_trimmed_mean",
    "latent_total": "total_latent_mean_trimmed",
}

FEATURE_COLUMNS = ["band_mean", "band_median", "band_min", "band_max", "band_range", "band_area"]


@dataclass
class BaselineConfig:
    manifest_path: Path
    band_features_path: Path
    dose_summary_path: Path
    output_table_path: Path
    predictions_path: Path
    manifest_output_path: Path
    plots_dir: Path
    random_seed: int


def load_manifest(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text())


def run_loo_predictions(model: LinearRegression, X: np.ndarray, y: np.ndarray) -> np.ndarray:
    predictions = np.zeros_like(y)
    n = X.shape[0]
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        model_clone = LinearRegression()
        model_clone.fit(X[mask], y[mask])
        predictions[i] = model_clone.predict(X[[i]])[0]
    return predictions


def plot_baseline(dose_ids: List[str], actual: np.ndarray, predicted: np.ndarray, target: str, plots_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(actual, predicted, color="#1f77b4")
    min_val = min(actual.min(), predicted.min())
    max_val = max(actual.max(), predicted.max())
    ax.plot([min_val, max_val], [min_val, max_val], color="black", linestyle="--")
    for dose, act, pred in zip(dose_ids, actual, predicted):
        ax.annotate(dose, (act, pred))
    ax.set_xlabel("Actual")
    ax.set_ylabel("Predicted")
    ax.set_title(f"Baseline Fit – {target}")
    ax.grid(alpha=0.3)
    plots_dir.mkdir(parents=True, exist_ok=True)
    outfile = plots_dir / f"baseline_fit_{target}.png"
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)
    return outfile.name


def run_baseline(config: BaselineConfig) -> None:
    manifest = load_manifest(config.manifest_path)
    band_df = pd.read_csv(config.band_features_path)
    dose_summary = pd.read_csv(config.dose_summary_path)

    merge_df = dose_summary[
        [
            "dose_id",
            "chrom_total_mg_per_gDW_trimmed_mean",
            "dad_total_mg_per_gDW_trimmed_mean",
            "total_latent_mean_trimmed",
        ]
    ]

    summary_records = []
    prediction_records = []
    plot_files = []

    for target_name, target_column in TARGET_COLUMNS.items():
        for kind in band_df["kind"].unique():
            subset = band_df[band_df["kind"] == kind].merge(merge_df, on="dose_id", how="inner")
            if subset.empty:
                continue

            X = subset[FEATURE_COLUMNS].to_numpy(dtype=float)
            y = subset[target_column].to_numpy(dtype=float)

            model = LinearRegression()
            model.fit(X, y)
            predictions = model.predict(X)
            loo_predictions = run_loo_predictions(model, X, y)

            rmse = float(np.sqrt(mean_squared_error(y, predictions)))
            mae = float(mean_absolute_error(y, predictions))
            r2 = float(r2_score(y, predictions))
            rmse_loo = float(np.sqrt(mean_squared_error(y, loo_predictions)))
            mae_loo = float(mean_absolute_error(y, loo_predictions))

            summary_records.append(
                {
                    "target": target_name,
                    "kind": kind,
                    "band_label": subset["band_label"].iloc[0],
                    "n_samples": int(len(subset)),
                    "rmse": rmse,
                    "mae": mae,
                    "r2": r2,
                    "rmse_loocv": rmse_loo,
                    "mae_loocv": mae_loo,
                    "coefficients": model.coef_.tolist(),
                    "intercept": float(model.intercept_),
                    "features": FEATURE_COLUMNS,
                }
            )

            for dose_id, actual, pred, pred_loo in zip(subset["dose_id"], y, predictions, loo_predictions):
                prediction_records.append(
                    {
                        "dose_id": dose_id,
                        "target": target_name,
                        "kind": kind,
                        "band_label": subset["band_label"].iloc[0],
                        "actual": float(actual),
                        "predicted": float(pred),
                        "predicted_loocv": float(pred_loo),
                    }
                )

            plot_name = plot_baseline(
                subset["dose_id"].tolist(), y, predictions, f"{target_name}_{kind}", config.plots_dir
            )
            plot_files.append(plot_name)

    summary_df = pd.DataFrame(summary_records)
    summary_df.to_csv(config.output_table_path, index=False)

    predictions_df = pd.DataFrame(prediction_records)
    predictions_df.to_csv(config.predictions_path, index=False)

    manifest_output = {
        "band_features": str(config.band_features_path),
        "dose_summary": str(config.dose_summary_path),
        "summary_table": str(config.output_table_path),
        "predictions_table": str(config.predictions_path),
        "plots_dir": str(config.plots_dir),
        "plots": plot_files,
        "targets": list(TARGET_COLUMNS.keys()),
        "features": FEATURE_COLUMNS,
        "random_seed": config.random_seed,
        "notes": "Linear regression baseline with leave-one-dose-out predictions.",
    }
    config.manifest_output_path.write_text(json.dumps(manifest_output, indent=2))
    print(f"[INFO] Baseline summary written to {config.output_table_path}")
    print(f"[INFO] Baseline manifest written to {config.manifest_output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run minimal baseline models for Phase 6.")
    parser.add_argument(
        "--manifest",
        default="analysis/phase6/inputs/phase6_input_manifest.json",
        help="Phase 6 input manifest.",
    )
    parser.add_argument(
        "--band-features",
        default="analysis/phase6/tables/band_features.csv",
        help="Path to band features CSV.",
    )
    parser.add_argument(
        "--dose-summary",
        default="canonical_dataset/dose_level_canonical_summary.csv",
        help="Dose-level summary table containing chrom/dad/latent targets.",
    )
    parser.add_argument(
        "--summary-output",
        default="analysis/phase6/tables/baseline_summary.csv",
        help="Output CSV path for baseline metrics.",
    )
    parser.add_argument(
        "--predictions-output",
        default="analysis/phase6/tables/baseline_predictions.csv",
        help="Output CSV path for predictions.",
    )
    parser.add_argument(
        "--manifest-output",
        default="analysis/phase6/tables/phase6_baseline_manifest.json",
        help="Output JSON manifest describing baseline run.",
    )
    parser.add_argument(
        "--plots-dir",
        default="analysis/phase6/plots",
        help="Directory for baseline plots.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = load_manifest(Path(args.manifest))
    random_seed = manifest.get("random_seed", 42)
    config = BaselineConfig(
        manifest_path=Path(args.manifest).resolve(),
        band_features_path=Path(args.band_features).resolve(),
        dose_summary_path=Path(args.dose_summary).resolve(),
        output_table_path=Path(args.summary_output).resolve(),
        predictions_path=Path(args.predictions_output).resolve(),
        manifest_output_path=Path(args.manifest_output).resolve(),
        plots_dir=Path(args.plots_dir).resolve(),
        random_seed=random_seed,
    )
    run_baseline(config)


if __name__ == "__main__":
    main()
