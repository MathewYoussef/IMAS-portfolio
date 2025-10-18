#!/usr/bin/env python3
"""
Phase 4 congruence analyses: PCA, Procrustes, RV, and PLS.

Consumes the cached Phase 4 matrices (phase4_data_cache.pkl) produced by
load_phase4_data.py and generates per-window/per-kind/per-modality
congruence metrics and plots.
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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import spatial
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA

PLOT_COLORS = {
    "Sigma": "#1f77b4",
    "Delta": "#ff7f0e",
    "12Oclock": "#2ca02c",
    "6Oclock": "#d62728",
}

MAX_COMPONENTS = 2


def load_json(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text())


@dataclass
class CongruenceConfig:
    cache_path: Path
    metadata_path: Path
    input_manifest: Path
    plots_dir: Path
    tables_dir: Path
    permutations: int
    random_seed: int


def load_cache(cache_path: Path) -> Dict[str, object]:
    with cache_path.open("rb") as handle:
        return pickle.load(handle)


def ensure_dirs(config: CongruenceConfig) -> None:
    config.plots_dir.mkdir(parents=True, exist_ok=True)
    config.tables_dir.mkdir(parents=True, exist_ok=True)


def compute_pca(matrix: np.ndarray, n_components: int = 2) -> Tuple[PCA, np.ndarray, np.ndarray]:
    components = min(n_components, min(matrix.shape))
    pca = PCA(n_components=components)
    scores = pca.fit_transform(matrix)
    loadings = pca.components_.T
    return pca, scores, loadings


def compute_procrustes(x: np.ndarray, y: np.ndarray, permutations: int, rng: random.Random) -> Dict[str, float]:
    mtx1, mtx2, disparity = spatial.procrustes(x, y)
    exceedances = 0
    for _ in range(permutations):
        permuted = rng.sample(list(range(mtx2.shape[0])), mtx2.shape[0])
        _, _, disp_perm = spatial.procrustes(mtx1, mtx2[permuted])
        if disp_perm <= disparity:
            exceedances += 1
    p_value = (exceedances + 1) / (permutations + 1)
    return {
        "disparity": float(disparity),
        "permutation_exceedances": exceedances,
        "permutation_p": float(p_value),
    }


def compute_rv(x: np.ndarray, y: np.ndarray, permutations: int, rng: random.Random) -> Dict[str, float]:
    x_centered = x - x.mean(axis=0, keepdims=True)
    y_centered = y - y.mean(axis=0, keepdims=True)
    numerator = np.trace(x_centered @ x_centered.T @ y_centered @ y_centered.T)
    denominator = math.sqrt(np.trace((x_centered @ x_centered.T) ** 2) * np.trace((y_centered @ y_centered.T) ** 2))
    rv_value = numerator / denominator if denominator != 0 else 0.0

    exceedances = 0
    for _ in range(permutations):
        permuted = rng.sample(list(range(y.shape[0])), y.shape[0])
        y_perm = y_centered[permuted]
        numerator_perm = np.trace(x_centered @ x_centered.T @ y_perm @ y_perm.T)
        denominator_perm = math.sqrt(np.trace((x_centered @ x_centered.T) ** 2) * np.trace((y_perm @ y_perm.T) ** 2))
        rv_perm = numerator_perm / denominator_perm if denominator_perm != 0 else 0.0
        if rv_perm >= rv_value:
            exceedances += 1
    p_value = (exceedances + 1) / (permutations + 1)
    return {
        "rv": float(rv_value),
        "permutation_exceedances": exceedances,
        "permutation_p": float(p_value),
    }


def compute_pls(x: np.ndarray, y: np.ndarray) -> Dict[str, object]:
    components = min(MAX_COMPONENTS, min(x.shape[1], y.shape[1], x.shape[0] - 1))
    if components < 1:
        return {"components": 0, "x_variance": None, "y_variance": None, "covariance": None, "warning": "Insufficient components"}
    pls = PLSRegression(n_components=components, scale=False)
    try:
        pls.fit(x, y)
        x_scores = pls.x_scores_
        y_scores = pls.y_scores_
        covariance = float(np.mean(pls.x_scores_ * pls.y_scores_))
        x_var = pls.x_weights_.var(axis=0).tolist()
        y_var = pls.y_weights_.var(axis=0).tolist()
        return {
            "components": components,
            "x_scores": x_scores.tolist(),
            "y_scores": y_scores.tolist(),
            "covariance": covariance,
            "x_variance": x_var,
            "y_variance": y_var,
            "warning": None,
        }
    except Exception as exc:
        return {"components": components, "warning": str(exc), "covariance": None, "x_variance": None, "y_variance": None}


def plot_pca_scores(
    scores_x: np.ndarray,
    scores_y: np.ndarray,
    dose_metadata: List[Dict[str, float]],
    window: str,
    kind: str,
    modality: str,
    plots_dir: Path,
    variance_x: List[float],
    variance_y: List[float],
) -> None:
    def components_available(scores: np.ndarray) -> int:
        return scores.shape[1]

    def extract_axes(scores: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if scores.shape[1] >= 2:
            return scores[:, 0], scores[:, 1]
        x = scores[:, 0]
        y = np.zeros_like(x)
        return x, y

    def variance_label(variance: List[float], index: int) -> float:
        return variance[index] if index < len(variance) else 0.0

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, scores, title, variance in [
        (axes[0], scores_x, "Reflectance", variance_x),
        (axes[1], scores_y, f"{modality.upper()}", variance_y),
    ]:
        x_vals, y_vals = extract_axes(scores)
        ax.scatter(x_vals, y_vals, color=PLOT_COLORS.get(kind, "gray"))
        for idx, meta in enumerate(dose_metadata):
            ax.annotate(meta["dose_id"], (x_vals[idx], y_vals[idx]))
        ax.set_title(
            f"{title} PCA Scores\n(PC1 {variance_label(variance, 0):.2f}, PC2 {variance_label(variance, 1):.2f})"
        )
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.grid(alpha=0.3)
    fig.suptitle(f"{kind} – {modality.upper()} ({window})")
    fig.tight_layout()
    outfile = plots_dir / f"pca_scores_{window}_{kind}_{modality}.png"
    fig.savefig(outfile, dpi=300)
    plt.close(fig)


def run_congruence(config: CongruenceConfig) -> None:
    ensure_dirs(config)
    cache = load_cache(config.cache_path)
    rng = random.Random(config.random_seed)

    congruence_records: List[Dict[str, object]] = []
    pca_score_tables: List[pd.DataFrame] = []
    pca_loading_tables: List[pd.DataFrame] = []

    for window, payload in cache["windows"].items():
        dose_metadata = payload["dose_metadata"]
        reflectance = payload["reflectance"]
        concentration = payload["concentration"]

        for kind, refl_info in reflectance.items():
            refl_matrix = np.array(refl_info["matrix"])
            refl_columns = refl_info["columns"]

            for modality, conc_info in concentration.items():
                conc_matrix = np.array(conc_info["matrix"])
                conc_columns = conc_info["columns"]

                pca_ref, scores_ref, loadings_ref = compute_pca(refl_matrix, MAX_COMPONENTS)
                pca_conc, scores_conc, loadings_conc = compute_pca(conc_matrix, MAX_COMPONENTS)

                components_used = min(pca_ref.n_components_, pca_conc.n_components_)
                scores_ref_used = scores_ref[:, :components_used]
                scores_conc_used = scores_conc[:, :components_used]
                variance_ref_used = pca_ref.explained_variance_ratio_[:components_used].tolist()
                variance_conc_used = pca_conc.explained_variance_ratio_[:components_used].tolist()

                shape_warning = None
                if scores_ref_used.shape != scores_conc_used.shape:
                    min_rows = min(scores_ref_used.shape[0], scores_conc_used.shape[0])
                    min_cols = min(scores_ref_used.shape[1], scores_conc_used.shape[1])
                    shape_warning = (
                        f"Adjusted shapes from {scores_ref_used.shape} and {scores_conc_used.shape} "
                        f"to {(min_rows, min_cols)} for congruence calculations."
                    )
                    scores_ref_used = scores_ref_used[:min_rows, :min_cols]
                    scores_conc_used = scores_conc_used[:min_rows, :min_cols]
                    variance_ref_used = variance_ref_used[:min_cols]
                    variance_conc_used = variance_conc_used[:min_cols]
                    components_used = min_cols

                if components_used == 0:
                    congruence_records.append(
                        {
                            "window": window,
                            "kind": kind,
                            "modality": modality,
                            "warning": "No usable PCA components after alignment.",
                        }
                    )
                    continue

                plot_pca_scores(
                    scores_ref_used,
                    scores_conc_used,
                    dose_metadata,
                    window,
                    kind,
                    modality,
                    config.plots_dir,
                    variance_x=variance_ref_used,
                    variance_y=variance_conc_used,
                )

                procrustes_stats = compute_procrustes(scores_ref_used, scores_conc_used, config.permutations, rng)
                rv_stats = compute_rv(scores_ref_used, scores_conc_used, config.permutations, rng)
                pls_stats = compute_pls(refl_matrix, conc_matrix)

                congruence_records.append(
                    {
                        "window": window,
                        "kind": kind,
                        "modality": modality,
                        "pca_components_reflectance": pca_ref.n_components_,
                        "pca_components_concentration": pca_conc.n_components_,
                        "pca_components_used": components_used,
                        "pca_variance_reflectance": pca_ref.explained_variance_ratio_.tolist(),
                        "pca_variance_concentration": pca_conc.explained_variance_ratio_.tolist(),
                        "pca_variance_reflectance_used": variance_ref_used,
                        "pca_variance_concentration_used": variance_conc_used,
                        "procrustes_disparity": procrustes_stats["disparity"],
                        "procrustes_permutation_exceedances": procrustes_stats["permutation_exceedances"],
                        "procrustes_permutation_p": procrustes_stats["permutation_p"],
                        "rv_value": rv_stats["rv"],
                        "rv_permutation_exceedances": rv_stats["permutation_exceedances"],
                        "rv_permutation_p": rv_stats["permutation_p"],
                        "pls_components": pls_stats["components"],
                        "pls_covariance": pls_stats["covariance"],
                        "pls_warning": pls_stats["warning"],
                        "permutations": config.permutations,
                        "random_seed": config.random_seed,
                        "zero_variance_reflectance": refl_info["scaling"]["zero_variance_columns"],
                        "zero_variance_concentration": conc_info["scaling"]["zero_variance_columns"],
                        "shape_warning": shape_warning,
                    }
                )

                scores_ref_df = pd.DataFrame(scores_ref_used, columns=[f"PC{i+1}" for i in range(components_used)])
                scores_ref_df["matrix"] = "reflectance"
                scores_ref_df["kind"] = kind
                scores_ref_df["modality"] = modality
                scores_ref_df["window"] = window
                scores_ref_df["dose_id"] = [meta["dose_id"] for meta in dose_metadata]
                pca_score_tables.append(scores_ref_df)

                scores_conc_df = pd.DataFrame(scores_conc_used, columns=[f"PC{i+1}" for i in range(components_used)])
                scores_conc_df["matrix"] = modality
                scores_conc_df["kind"] = kind
                scores_conc_df["modality"] = modality
                scores_conc_df["window"] = window
                scores_conc_df["dose_id"] = [meta["dose_id"] for meta in dose_metadata]
                pca_score_tables.append(scores_conc_df)

                loadings_ref_df = pd.DataFrame(loadings_ref[:, :components_used], columns=[f"PC{i+1}" for i in range(components_used)])
                loadings_ref_df["feature"] = refl_columns
                loadings_ref_df["matrix"] = "reflectance"
                loadings_ref_df["kind"] = kind
                loadings_ref_df["modality"] = modality
                loadings_ref_df["window"] = window
                pca_loading_tables.append(loadings_ref_df)

                loadings_conc_df = pd.DataFrame(loadings_conc[:, :components_used], columns=[f"PC{i+1}" for i in range(components_used)])
                loadings_conc_df["feature"] = conc_columns
                loadings_conc_df["matrix"] = modality
                loadings_conc_df["kind"] = kind
                loadings_conc_df["modality"] = modality
                loadings_conc_df["window"] = window
                pca_loading_tables.append(loadings_conc_df)

    congruence_df = pd.DataFrame(congruence_records)
    congruence_path = config.tables_dir / "congruence_summary.csv"
    congruence_df.to_csv(congruence_path, index=False)

    scores_path = config.tables_dir / "pca_scores.csv"
    pd.concat(pca_score_tables, ignore_index=True).to_csv(scores_path, index=False)

    loadings_path = config.tables_dir / "pca_loadings.csv"
    pd.concat(pca_loading_tables, ignore_index=True).to_csv(loadings_path, index=False)

    manifest = {
        "cache_path": str(config.cache_path),
        "metadata_path": str(config.metadata_path),
        "plots_dir": str(config.plots_dir),
        "tables_dir": str(config.tables_dir),
        "settings": {
            "permutations": config.permutations,
            "random_seed": config.random_seed,
            "max_components": MAX_COMPONENTS,
        },
        "outputs": {
            "congruence_summary": str(congruence_path),
            "pca_scores": str(scores_path),
            "pca_loadings": str(loadings_path),
            "plots": sorted(p.name for p in config.plots_dir.glob("*.png")),
        },
    }
    manifest_path = config.tables_dir / "phase4_congruence_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"[INFO] Congruence summary written to {congruence_path}")
    print(f"[INFO] Manifest written to {manifest_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run PCA/Procrustes/RV/PLS congruence analyses for Phase 4.")
    parser.add_argument(
        "--cache-path",
        default="analysis/phase4/inputs/phase4_data_cache.pkl",
        help="Path to cached matrices (pickle).",
    )
    parser.add_argument(
        "--metadata-path",
        default="analysis/phase4/inputs/phase4_matrix_manifest.json",
        help="Path to matrix metadata JSON.",
    )
    parser.add_argument(
        "--input-manifest",
        default="analysis/phase4/inputs/phase4_input_manifest.json",
        help="Phase 4 input manifest (for permutation/seed defaults).",
    )
    parser.add_argument(
        "--plots-dir",
        default="analysis/phase4/plots",
        help="Directory for PCA score plots.",
    )
    parser.add_argument(
        "--tables-dir",
        default="analysis/phase4/tables",
        help="Directory for congruence tables.",
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
    input_manifest = load_json(Path(args.input_manifest).resolve())
    permutations = args.permutations or input_manifest.get("default_permutations", 10000)
    random_seed = args.random_seed or input_manifest.get("random_seed", 42)

    config = CongruenceConfig(
        cache_path=Path(args.cache_path).resolve(),
        metadata_path=Path(args.metadata_path).resolve(),
        input_manifest=Path(args.input_manifest).resolve(),
        plots_dir=Path(args.plots_dir).resolve(),
        tables_dir=Path(args.tables_dir).resolve(),
        permutations=permutations,
        random_seed=random_seed,
    )
    run_congruence(config)


if __name__ == "__main__":
    main()
