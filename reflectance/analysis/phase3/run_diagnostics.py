#!/usr/bin/env python3
"""
Phase 3 diagnostics runner.

Generates scatter/regression plots, Bland–Altman plots, dispersion comparisons,
and correlation/regression tables for each reflectance–concentration pairing
described in phase3_concentration_target_manifest.json.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

PLOT_STYLE = {
    "Sigma": {"color": "#1f77b4", "marker": "o"},
    "Delta": {"color": "#ff7f0e", "marker": "s"},
    "12Oclock": {"color": "#2ca02c", "marker": "^"},
    "6Oclock": {"color": "#d62728", "marker": "v"},
}


@dataclass
class DiagnosticsConfig:
    join_manifest: Path
    output_dir: Path
    plots_dir: Path
    tables_dir: Path
    permutations: int
    random_seed: int


def load_manifest(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text())


def ensure_dirs(config: DiagnosticsConfig) -> None:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    config.plots_dir.mkdir(parents=True, exist_ok=True)
    config.tables_dir.mkdir(parents=True, exist_ok=True)


def scatter_plot(df: pd.DataFrame, x_col: str, y_col: str, title: str, output_path: Path) -> None:
    plt.figure(figsize=(6, 5))
    for kind, group in df.groupby("kind"):
        style = PLOT_STYLE.get(kind, {"color": "gray", "marker": "o"})
        plt.scatter(group[x_col], group[y_col], label=kind, color=style["color"], marker=style["marker"])
        for _, row in group.iterrows():
            plt.annotate(row["dose_id"], (row[x_col], row[y_col]), fontsize=8, alpha=0.7)
    slope, intercept, r_value, p_value, _ = stats.linregress(df[x_col], df[y_col])
    x_vals = np.linspace(df[x_col].min(), df[x_col].max(), 100)
    plt.plot(x_vals, intercept + slope * x_vals, color="black", linestyle="--", linewidth=1, alpha=0.7)
    plt.title(title)
    plt.xlabel(x_col.replace("_", " ").title())
    plt.ylabel(y_col.replace("_", " ").title())
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def bland_altman_plot(df: pd.DataFrame, x_col: str, y_col: str, title: str, output_path: Path) -> Tuple[float, float, float]:
    mean_values = 0.5 * (df[x_col] + df[y_col])
    diff_values = df[y_col] - df[x_col]
    mean_bias = float(diff_values.mean())
    std_diff = float(diff_values.std(ddof=1))
    loa_upper = mean_bias + 1.96 * std_diff
    loa_lower = mean_bias - 1.96 * std_diff

    plt.figure(figsize=(6, 5))
    plt.axhline(mean_bias, color="black", linestyle="-", linewidth=1, label=f"Bias {mean_bias:.3f}")
    plt.axhline(loa_upper, color="red", linestyle="--", linewidth=1, label=f"+1.96 SD {loa_upper:.3f}")
    plt.axhline(loa_lower, color="red", linestyle="--", linewidth=1, label=f"-1.96 SD {loa_lower:.3f}")

    for kind, group in df.groupby("kind"):
        style = PLOT_STYLE.get(kind, {"color": "gray", "marker": "o"})
        plt.scatter(mean_values[group.index], diff_values[group.index], label=kind, color=style["color"], marker=style["marker"])
        for _, row in group.iterrows():
            idx = row.name
            plt.annotate(row["dose_id"], (mean_values[idx], diff_values[idx]), fontsize=8, alpha=0.7)

    plt.title(title)
    plt.xlabel("Mean of measurements")
    plt.ylabel(f"{y_col} - {x_col}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    return mean_bias, loa_lower, loa_upper


def compute_correlations(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    permutations: int,
    rng: random.Random,
) -> Dict[str, float]:
    x = df[x_col].to_numpy()
    y = df[y_col].to_numpy()
    pearson_r, pearson_p = stats.pearsonr(x, y)
    spearman_rho, spearman_p = stats.spearmanr(x, y)
    kendall_tau, kendall_p = stats.kendalltau(x, y)

    def permute_stat(stat_func):
        observed = stat_func(x, y)
        count = 0
        for _ in range(permutations):
            y_perm = rng.sample(list(y), len(y))
            permuted = stat_func(x, y_perm)
            if abs(permuted) >= abs(observed):
                count += 1
        return observed, (count + 1) / (permutations + 1)

    pearson_stat, pearson_perm_p = permute_stat(lambda a, b: stats.pearsonr(a, b)[0])
    spearman_stat, spearman_perm_p = permute_stat(lambda a, b: stats.spearmanr(a, b)[0])
    kendall_stat, kendall_perm_p = permute_stat(lambda a, b: stats.kendalltau(a, b)[0])

    return {
        "pearson_r": float(pearson_r),
        "pearson_p": float(pearson_p),
        "pearson_perm_p": float(pearson_perm_p),
        "spearman_rho": float(spearman_rho),
        "spearman_p": float(spearman_p),
        "spearman_perm_p": float(spearman_perm_p),
        "kendall_tau": float(kendall_tau),
        "kendall_p": float(kendall_p),
        "kendall_perm_p": float(kendall_perm_p),
    }


def dispersion_table(
    df: pd.DataFrame,
    modality: str,
    columns: List[Tuple[str, str]],
) -> pd.DataFrame:
    records = []
    for (refl_col, target_col) in columns:
        for _, row in df.iterrows():
            records.append(
                {
                    "dose_id": row["dose_id"],
                    "kind": row["kind"],
                    "modality": modality,
                    "reflectance_metric": refl_col,
                    "reflectance_value": row.get(refl_col),
                    "target_metric": target_col,
                    "target_value": row.get(target_col),
                }
            )
    return pd.DataFrame(records)


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    df.to_csv(output_path, index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 3 diagnostics (plots, correlations, dispersion checks).")
    parser.add_argument(
        "--join-manifest",
        default="analysis/phase3/phase3_concentration_target_manifest.json",
        help="Manifest describing join tables for reflectance vs concentration.",
    )
    parser.add_argument(
        "--plots-dir",
        default="analysis/phase3/diagnostics/plots",
        help="Directory to store generated plots.",
    )
    parser.add_argument(
        "--tables-dir",
        default="analysis/phase3/diagnostics/tables",
        help="Directory to store statistical tables.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=10000,
        help="Number of permutations for correlation significance.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=42,
        help="Seed for permutation reproducibility.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = DiagnosticsConfig(
        join_manifest=Path(args.join_manifest).resolve(),
        output_dir=Path("analysis/phase3/diagnostics").resolve(),
        plots_dir=Path(args.plots_dir).resolve(),
        tables_dir=Path(args.tables_dir).resolve(),
        permutations=args.permutations,
        random_seed=args.random_seed,
    )

    ensure_dirs(config)

    rng = random.Random(config.random_seed)
    manifest = load_manifest(config.join_manifest)

    correlation_records: List[Dict[str, object]] = []
    regression_records: List[Dict[str, object]] = []
    bland_altman_records: List[Dict[str, object]] = []
    dispersion_records: List[pd.DataFrame] = []

    modality_info = [
        ("chrom", "chrom_join", ["chrom_total_mean", "chrom_oxidized_mean", "chrom_reduced_mean"]),
        ("dad", "dad_join", ["dad_total_mean", "dad_oxidized_mean", "dad_reduced_mean"]),
        ("latent", "latent_join", ["latent_total_mean", "latent_oxidized_mean", "latent_reduced_mean"]),
    ]

    for entry in manifest.get("outputs", []):
        window = entry["window"]
        window_label = f"{int(window['lower_nm'])}_{int(window['upper_nm'])}nm"

        for modality, key, metrics in modality_info:
            join_path = Path(entry[key]).resolve()
            if not join_path.exists():
                print(f"[WARN] Join file missing: {join_path}")
                continue
            df = pd.read_csv(join_path)

            for metric in metrics:
                scatter_path = config.plots_dir / f"{modality}_scatter_{metric}_{window_label}.png"
                title = f"{modality.upper()} {metric} vs Reflectance ({window_label})"
                scatter_plot(df, "mean_reflectance", metric, title, scatter_path)

                bland_path = config.plots_dir / f"{modality}_bland_altman_{metric}_{window_label}.png"
                bias, loa_lower, loa_upper = bland_altman_plot(df, "mean_reflectance", metric, title, bland_path)
                bland_altman_records.append(
                    {
                        "window": window_label,
                        "modality": modality,
                        "metric": metric,
                        "bias": bias,
                        "loa_lower": loa_lower,
                        "loa_upper": loa_upper,
                    }
                )

                corr_stats = compute_correlations(df, "mean_reflectance", metric, config.permutations, rng)
                corr_stats.update(
                    {
                        "window": window_label,
                        "modality": modality,
                        "metric": metric,
                    }
                )
                correlation_records.append(corr_stats)

                slope, intercept, r_value, p_value, stderr = stats.linregress(df["mean_reflectance"], df[metric])
                regression_records.append(
                    {
                        "window": window_label,
                        "modality": modality,
                        "metric": metric,
                        "slope": float(slope),
                        "intercept": float(intercept),
                        "r_value": float(r_value),
                        "p_value": float(p_value),
                        "stderr": float(stderr),
                    }
                )

            if modality == "latent":
                dispersion_pairs = [
                    ("mean_std", "latent_total_sd"),
                    ("mean_mad", None),
                ]
            elif modality == "dad":
                dispersion_pairs = [
                    ("mean_std", "dad_total_sd"),
                    ("mean_mad", "dad_total_mad"),
                ]
            else:  # chrom
                dispersion_pairs = [
                    ("mean_std", "chrom_total_sd"),
                    ("mean_mad", None),
                ]

            disp_df = dispersion_table(df, modality, dispersion_pairs)
            disp_df["window"] = window_label
            dispersion_records.append(disp_df)

    corr_table = pd.DataFrame(correlation_records)
    corr_path = config.tables_dir / "correlations.csv"
    save_table(corr_table, corr_path)

    regression_table = pd.DataFrame(regression_records)
    regression_path = config.tables_dir / "regression_results.csv"
    save_table(regression_table, regression_path)

    bland_table = pd.DataFrame(bland_altman_records)
    bland_path = config.tables_dir / "bland_altman_summary.csv"
    save_table(bland_table, bland_path)

    dispersion_combined = pd.concat(dispersion_records, ignore_index=True) if dispersion_records else pd.DataFrame()
    dispersion_path = config.tables_dir / "dispersion_comparisons.csv"
    save_table(dispersion_combined, dispersion_path)

    diagnostics_manifest = {
        "join_manifest": str(config.join_manifest),
        "plots_dir": str(config.plots_dir),
        "tables": {
            "correlations": str(corr_path),
            "regression": str(regression_path),
            "bland_altman": str(bland_path),
            "dispersion": str(dispersion_path),
        },
        "plots_generated": sorted([p.name for p in config.plots_dir.glob("*.png")]),
        "settings": {
            "permutations": config.permutations,
            "random_seed": config.random_seed,
        },
    }

    manifest_path = config.output_dir / "phase3_diagnostics_manifest.json"
    manifest_path.write_text(json.dumps(diagnostics_manifest, indent=2))
    print(f"[INFO] Diagnostics manifest written to {manifest_path}")


if __name__ == "__main__":
    main()
