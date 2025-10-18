#!/usr/bin/env python3
"""
Phase 4 dependence metrics (distance correlation and HSIC).

Consumes the cached Phase 4 matrices and computes dCor/HSIC for each
window/kind/modality combination with permutation-backed p-values.
"""

from __future__ import annotations

import argparse
import json
import math
import pickle
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist, squareform


@dataclass
class DependenceConfig:
    cache_path: Path
    metadata_path: Path
    input_manifest: Path
    tables_dir: Path
    plots_dir: Path
    permutations: int
    random_seed: int


def ensure_dirs(config: DependenceConfig) -> None:
    config.tables_dir.mkdir(parents=True, exist_ok=True)
    config.plots_dir.mkdir(parents=True, exist_ok=True)


def load_cache(cache_path: Path) -> Dict[str, object]:
    with cache_path.open("rb") as handle:
        return pickle.load(handle)


def extract_matrix(info: Dict[str, object]) -> Tuple[np.ndarray, List[str], List[str]]:
    matrix = np.array(info["matrix"])
    columns = info["columns"]
    zero_cols = info["scaling"]["zero_variance_columns"]
    keep_mask = [col not in zero_cols for col in columns]
    kept_cols = [col for col, keep in zip(columns, keep_mask) if keep]
    dropped = [col for col, keep in zip(columns, keep_mask) if not keep]
    if kept_cols:
        matrix = matrix[:, [i for i, keep in enumerate(keep_mask) if keep]]
    else:
        matrix = np.empty((matrix.shape[0], 0))
    return matrix, kept_cols, dropped


def distance_correlation(x: np.ndarray, y: np.ndarray) -> float:
    if x.shape[0] != y.shape[0]:
        raise ValueError("distance_correlation requires matrices with the same number of rows.")
    n = x.shape[0]
    if n <= 1:
        return 0.0

    dist_x = cdist(x, x, metric="euclidean")
    dist_y = cdist(y, y, metric="euclidean")

    def double_center(dist: np.ndarray) -> np.ndarray:
        row_mean = dist.mean(axis=1, keepdims=True)
        col_mean = dist.mean(axis=0, keepdims=True)
        total_mean = dist.mean()
        return dist - row_mean - col_mean + total_mean

    A = double_center(dist_x)
    B = double_center(dist_y)

    d_cov = (A * B).mean()
    d_var_x = (A * A).mean()
    d_var_y = (B * B).mean()

    if d_var_x <= 0 or d_var_y <= 0:
        return 0.0
    return float(d_cov / math.sqrt(d_var_x * d_var_y))


def hsic_rbf(x: np.ndarray, y: np.ndarray, sigma_x: float, sigma_y: float) -> float:
    if x.shape[0] != y.shape[0]:
        raise ValueError("HSIC requires matrices with the same number of rows.")
    n = x.shape[0]
    if n <= 1:
        return 0.0

    def rbf_kernel(data: np.ndarray, sigma: float) -> np.ndarray:
        if sigma <= 0:
            sigma = 1.0
        gamma = 1.0 / (2.0 * sigma ** 2)
        dist_sq = cdist(data, data, metric="sqeuclidean")
        return np.exp(-gamma * dist_sq)

    K = rbf_kernel(x, sigma_x)
    L = rbf_kernel(y, sigma_y)

    H = np.eye(n) - np.ones((n, n)) / n
    Kc = H @ K @ H
    Lc = H @ L @ H
    return float(np.trace(Kc @ Lc) / ((n - 1) ** 2))


def median_bandwidth(data: np.ndarray) -> float:
    if data.shape[0] <= 1:
        return 1.0
    distances = pdist(data, metric="euclidean")
    median = np.median(distances)
    if median <= 0:
        return 1.0
    return float(median)


def run_dependence(config: DependenceConfig) -> None:
    ensure_dirs(config)
    cache = load_cache(config.cache_path)
    rng = random.Random(config.random_seed)

    records: List[Dict[str, object]] = []

    for window, payload in cache["windows"].items():
        reflectance = payload["reflectance"]
        concentration = payload["concentration"]

        for kind, refl_info in reflectance.items():
            refl_matrix, refl_cols, refl_dropped = extract_matrix(refl_info)
            if refl_matrix.shape[1] == 0:
                records.append(
                    {
                        "window": window,
                        "kind": kind,
                        "modality": None,
                        "metric": "dCor",
                        "statistic": None,
                        "permutations": config.permutations,
                        "permutation_exceedances": None,
                        "permutation_p": None,
                        "kernel_bandwidth_x": None,
                        "kernel_bandwidth_y": None,
                        "seed": config.random_seed,
                        "reflectance_features": [],
                        "concentration_features": [],
                        "dropped_reflectance_columns": refl_dropped,
                        "dropped_concentration_columns": [],
                        "warning": "No usable reflectance features after zero-variance filtering.",
                    }
                )
                continue

            sigma_x = median_bandwidth(refl_matrix)

            for modality, conc_info in concentration.items():
                conc_matrix, conc_cols, conc_dropped = extract_matrix(conc_info)
                if conc_matrix.shape[1] == 0:
                    records.append(
                        {
                            "window": window,
                            "kind": kind,
                            "modality": modality,
                            "metric": "dCor",
                            "statistic": None,
                            "permutations": config.permutations,
                            "permutation_exceedances": None,
                            "permutation_p": None,
                            "kernel_bandwidth_x": sigma_x,
                            "kernel_bandwidth_y": None,
                            "seed": config.random_seed,
                            "reflectance_features": refl_cols,
                            "concentration_features": conc_cols,
                            "dropped_reflectance_columns": refl_dropped,
                            "dropped_concentration_columns": conc_dropped,
                            "warning": "No usable concentration features after zero-variance filtering.",
                        }
                    )
                    continue

                if refl_matrix.shape[0] != conc_matrix.shape[0]:
                    min_rows = min(refl_matrix.shape[0], conc_matrix.shape[0])
                    refl_matrix_used = refl_matrix[:min_rows]
                    conc_matrix_used = conc_matrix[:min_rows]
                    shape_warning = f"Adjusted rows from {refl_matrix.shape[0]} and {conc_matrix.shape[0]} to {min_rows}."
                else:
                    refl_matrix_used = refl_matrix
                    conc_matrix_used = conc_matrix
                    shape_warning = None

                obs_dcor = distance_correlation(refl_matrix_used, conc_matrix_used)
                exceed_dcor = 0
                for _ in range(config.permutations):
                    idx = rng.sample(range(conc_matrix_used.shape[0]), conc_matrix_used.shape[0])
                    perm_stat = distance_correlation(refl_matrix_used, conc_matrix_used[idx])
                    if perm_stat >= obs_dcor:
                        exceed_dcor += 1
                p_dcor = (exceed_dcor + 1) / (config.permutations + 1)

                sigma_y = median_bandwidth(conc_matrix_used)
                obs_hsic = hsic_rbf(refl_matrix_used, conc_matrix_used, sigma_x, sigma_y)
                exceed_hsic = 0
                for _ in range(config.permutations):
                    idx = rng.sample(range(conc_matrix_used.shape[0]), conc_matrix_used.shape[0])
                    perm_stat = hsic_rbf(refl_matrix_used, conc_matrix_used[idx], sigma_x, sigma_y)
                    if perm_stat >= obs_hsic:
                        exceed_hsic += 1
                p_hsic = (exceed_hsic + 1) / (config.permutations + 1)

                for metric, stat, exceed, p_val, bw_x, bw_y in [
                    ("dCor", obs_dcor, exceed_dcor, p_dcor, None, None),
                    ("HSIC", obs_hsic, exceed_hsic, p_hsic, sigma_x, sigma_y),
                ]:
                    records.append(
                        {
                            "window": window,
                            "kind": kind,
                            "modality": modality,
                            "metric": metric,
                            "statistic": stat,
                            "permutations": config.permutations,
                            "permutation_exceedances": exceed,
                            "permutation_p": p_val,
                            "kernel_bandwidth_x": bw_x,
                            "kernel_bandwidth_y": bw_y,
                            "seed": config.random_seed,
                            "reflectance_features": refl_cols,
                            "concentration_features": conc_cols,
                            "dropped_reflectance_columns": refl_dropped,
                            "dropped_concentration_columns": conc_dropped,
                            "warning": shape_warning,
                        }
                    )

    dependence_df = pd.DataFrame(records)
    table_path = config.tables_dir / "dependence_metrics.csv"
    dependence_df.to_csv(table_path, index=False)

    manifest = {
        "cache_path": str(config.cache_path),
        "metadata_path": str(config.metadata_path),
        "tables_dir": str(config.tables_dir),
        "plots_dir": str(config.plots_dir),
        "settings": {
            "permutations": config.permutations,
            "random_seed": config.random_seed,
        },
        "outputs": {
            "dependence_metrics": str(table_path),
            "plots": [],
        },
    }
    manifest_path = config.tables_dir / "phase4_dependence_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"[INFO] Dependence metrics written to {table_path}")
    print(f"[INFO] Manifest written to {manifest_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute dCor/HSIC dependence metrics for Phase 4.")
    parser.add_argument(
        "--cache-path",
        default="analysis/phase4/inputs/phase4_data_cache.pkl",
        help="Path to cached matrices pickle.",
    )
    parser.add_argument(
        "--metadata-path",
        default="analysis/phase4/inputs/phase4_matrix_manifest.json",
        help="Path to matrix metadata JSON.",
    )
    parser.add_argument(
        "--input-manifest",
        default="analysis/phase4/inputs/phase4_input_manifest.json",
        help="Phase 4 input manifest for default settings.",
    )
    parser.add_argument(
        "--tables-dir",
        default="analysis/phase4/tables",
        help="Directory to write dependence tables.",
    )
    parser.add_argument(
        "--plots-dir",
        default="analysis/phase4/plots",
        help="Directory to write optional dependence plots.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=None,
        help="Override permutation count (defaults to manifest value).",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="Override permutation seed (defaults to manifest value).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_manifest = json.loads(Path(args.input_manifest).read_text())
    permutations = args.permutations or input_manifest.get("default_permutations", 10000)
    random_seed = args.random_seed or input_manifest.get("random_seed", 42)

    config = DependenceConfig(
        cache_path=Path(args.cache_path).resolve(),
        metadata_path=Path(args.metadata_path).resolve(),
        input_manifest=Path(args.input_manifest).resolve(),
        tables_dir=Path(args.tables_dir).resolve(),
        plots_dir=Path(args.plots_dir).resolve(),
        permutations=permutations,
        random_seed=random_seed,
    )
    run_dependence(config)


if __name__ == "__main__":
    main()
