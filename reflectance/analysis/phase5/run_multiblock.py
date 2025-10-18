#!/usr/bin/env python3
"""
Phase 5 multiblock fusion (stacked PLS approximation).

This script treats reflectance blocks (per-angle and optional composites) as X
and concentration modalities as Y, fits a PLSRegression model for each
window/modality, and records joint scores/loadings along with diagnostics.

NOTE: Python lacks a robust O2PLS implementation; this script uses a stacked
PLS approach (per block) as an approximation. Interpret results accordingly.
"""

from __future__ import annotations

import argparse
import json
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression


@dataclass
class FusionConfig:
    block_cache: Path
    input_manifest: Path
    tables_dir: Path
    plots_dir: Path
    include_exploratory: bool
    joint_components: int | None
    random_seed: int
    target_modalities: List[str] | None


def ensure_dirs(config: FusionConfig) -> None:
    config.tables_dir.mkdir(parents=True, exist_ok=True)
    config.plots_dir.mkdir(parents=True, exist_ok=True)


def load_cache(path: Path) -> Dict[str, object]:
    with path.open("rb") as handle:
        return pickle.load(handle)


def matrix_rank(matrix: np.ndarray) -> int:
    if matrix.size == 0:
        return 0
    return int(np.linalg.matrix_rank(matrix))


def stacked_matrix(blocks: Dict[str, Dict[str, object]], block_names: List[str]) -> Tuple[np.ndarray, List[str], Dict[str, Tuple[int, int]]]:
    matrices = []
    column_labels = []
    spans = {}
    start = 0
    for name in block_names:
        if name not in blocks:
            continue
        info = blocks[name]
        matrix = np.array(info["matrix"])
        matrices.append(matrix)
        cols = info["columns"]
        column_labels.extend([f"{name}:{col}" for col in cols])
        end = start + len(cols)
        spans[name] = (start, end)
        start = end
    if not matrices:
        return np.empty((0, 0)), [], {}
    stacked = np.hstack(matrices)
    return stacked, column_labels, spans


def adjust_components(requested: int, x_rank: int, y_rank: int) -> Tuple[int, str | None]:
    max_components = max(0, min(requested, x_rank, y_rank))
    warning = None
    if max_components < requested:
        warning = f"Joint components adjusted from {requested} to {max_components} due to rank limits ({x_rank}, {y_rank})."
    return max_components, warning


def plot_scores(window: str, modality: str, scores_x: np.ndarray, scores_y: np.ndarray, plots_dir: Path) -> str | None:
    if scores_x.shape[1] == 0 or scores_y.shape[1] == 0:
        return None
    pc1_x = scores_x[:, 0]
    pc1_y = scores_y[:, 0]
    pc2_x = scores_x[:, 1] if scores_x.shape[1] > 1 else np.zeros_like(pc1_x)
    pc2_y = scores_y[:, 1] if scores_y.shape[1] > 1 else np.zeros_like(pc1_y)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(pc1_x, pc1_y, color="#1f77b4", label="Component 1")
    ax.scatter(pc2_x, pc2_y, color="#ff7f0e", label="Component 2")
    ax.set_xlabel("Reflectance scores")
    ax.set_ylabel(f"{modality.upper()} scores")
    ax.set_title(f"PLS Scores ({window}, {modality})")
    ax.grid(alpha=0.3)
    ax.legend()
    outfile = plots_dir / f"multiblock_scores_{window}_{modality}.png"
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)
    return outfile.name


def run_fusion(config: FusionConfig) -> None:
    ensure_dirs(config)
    cache = load_cache(config.block_cache)
    manifest = json.loads(Path(config.input_manifest).read_text())
    default_components = manifest["default_components"]

    records_summary: List[Dict[str, object]] = []
    score_tables: List[pd.DataFrame] = []
    loading_tables: List[pd.DataFrame] = []
    plot_files: List[str] = []

    for window, entry in cache["blocks"].items():
        reflectance_blocks = entry["reflectance"]
        concentration_blocks = entry["concentration"]

        primary_blocks = manifest["reflectance_blocks"]["primary"]
        exploratory_blocks = manifest["reflectance_blocks"]["exploratory"] if config.include_exploratory else []
        selected_blocks = [b for b in primary_blocks if b in reflectance_blocks] + [b for b in exploratory_blocks if b in reflectance_blocks]

        if not selected_blocks:
            records_summary.append(
                {
                    "window": window,
                    "modality": None,
                    "blocks_used": [],
                    "warning": "No usable reflectance blocks.",
                }
            )
            continue

        x_matrix, x_labels, block_spans = stacked_matrix(reflectance_blocks, selected_blocks)
        x_rank = matrix_rank(x_matrix)

        modalities = config.target_modalities or list(concentration_blocks.keys())
        for modality in modalities:
            if modality not in concentration_blocks:
                continue

            y_info = concentration_blocks[modality]
            y_matrix = np.array(y_info["matrix"])
            y_labels = y_info["columns"]
            y_rank = matrix_rank(y_matrix)

            if x_matrix.size == 0 or y_matrix.size == 0:
                records_summary.append(
                    {
                        "window": window,
                        "modality": modality,
                        "blocks_used": selected_blocks,
                        "warning": "Empty X or Y matrix after filtering.",
                    }
                )
                continue

            requested_joint = config.joint_components or default_components.get("joint", 1)
            used_components, warning = adjust_components(requested_joint, x_rank, y_rank)
            if used_components == 0:
                records_summary.append(
                    {
                        "window": window,
                        "modality": modality,
                        "blocks_used": selected_blocks,
                        "warning": "No components retained (rank deficiency)."
                    }
                )
                continue

            pls = PLSRegression(n_components=used_components, scale=False)
            pls.fit(x_matrix, y_matrix)
            x_scores = pls.x_scores_
            y_scores = pls.y_scores_
            x_weights = pls.x_weights_
            y_weights = pls.y_weights_
            y_loadings = pls.y_loadings_

            plot_name = plot_scores(window, modality, x_scores, y_scores, config.plots_dir)
            if plot_name:
                plot_files.append(plot_name)

            record = {
                "window": window,
                "modality": modality,
                "blocks_used": selected_blocks,
                "x_rank": x_rank,
                "y_rank": y_rank,
                "requested_joint_components": requested_joint,
                "used_joint_components": used_components,
                "exploratory_blocks_included": config.include_exploratory,
                "warning": warning,
            }
            records_summary.append(record)

            scores_df = pd.DataFrame(x_scores, columns=[f"comp_{i+1}" for i in range(used_components)])
            scores_df["matrix"] = "reflectance"
            scores_df["window"] = window
            scores_df["modality"] = modality
            score_tables.append(scores_df)

            scores_y_df = pd.DataFrame(y_scores, columns=[f"comp_{i+1}" for i in range(used_components)])
            scores_y_df["matrix"] = modality
            scores_y_df["window"] = window
            scores_y_df["modality"] = modality
            score_tables.append(scores_y_df)

            loading_df = pd.DataFrame(x_weights, columns=[f"comp_{i+1}" for i in range(used_components)])
            loading_df["feature"] = x_labels[:x_weights.shape[0]]
            loading_df["matrix"] = "reflectance"
            loading_df["window"] = window
            loading_df["modality"] = modality
            loading_tables.append(loading_df)

            y_loading_df = pd.DataFrame(y_loadings, columns=[f"comp_{i+1}" for i in range(used_components)])
            y_loading_df["feature"] = y_labels[:y_loadings.shape[0]]
            y_loading_df["matrix"] = modality
            y_loading_df["window"] = window
            y_loading_df["modality"] = modality
            loading_tables.append(y_loading_df)

    summary_df = pd.DataFrame(records_summary)
    summary_path = config.tables_dir / "multiblock_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    if score_tables:
        scores_path = config.tables_dir / "multiblock_scores.csv"
        pd.concat(score_tables, ignore_index=True).to_csv(scores_path, index=False)
    else:
        scores_path = None

    if loading_tables:
        loadings_path = config.tables_dir / "multiblock_loadings.csv"
        pd.concat(loading_tables, ignore_index=True).to_csv(loadings_path, index=False)
    else:
        loadings_path = None

    manifest = {
        "block_cache": str(config.block_cache),
        "input_manifest": str(config.input_manifest),
        "tables_dir": str(config.tables_dir),
        "plots_dir": str(config.plots_dir),
        "settings": {
            "include_exploratory": config.include_exploratory,
            "joint_components": config.joint_components,
            "random_seed": config.random_seed,
        },
        "outputs": {
            "summary": str(summary_path),
            "scores": str(scores_path) if scores_path else None,
            "loadings": str(loadings_path) if loadings_path else None,
            "plots": plot_files,
        },
        "notes": "PLS approximation used; true O2PLS not available. Component counts capped by block ranks.",
    }
    manifest_path = config.tables_dir / "phase5_multiblock_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"[INFO] Multiblock summary written to {summary_path}")
    print(f"[INFO] Manifest written to {manifest_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run stacked PLS multiblock fusion for Phase 5.")
    parser.add_argument(
        "--block-cache",
        default="analysis/phase5/inputs/phase5_block_cache.pkl",
        help="Path to Phase 5 block cache.",
    )
    parser.add_argument(
        "--input-manifest",
        default="analysis/phase5/inputs/phase5_input_manifest.json",
        help="Phase 5 manifest with settings.",
    )
    parser.add_argument(
        "--tables-dir",
        default="analysis/phase5/tables",
        help="Directory to write table outputs.",
    )
    parser.add_argument(
        "--plots-dir",
        default="analysis/phase5/plots",
        help="Directory to write plots.",
    )
    parser.add_argument(
        "--include-exploratory",
        action="store_true",
        help="Include exploratory reflectance blocks (Sigma, Delta).",
    )
    parser.add_argument(
        "--joint-components",
        type=int,
        default=None,
        help="Override joint component count (defaults to manifest).",
    )
    parser.add_argument(
        "--modalities",
        nargs="+",
        default=None,
        help="Limit to specific concentration modalities (e.g., chrom dad latent).",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="Override seed (defaults to manifest).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = json.loads(Path(args.input_manifest).read_text())

    config = FusionConfig(
        block_cache=Path(args.block_cache).resolve(),
        input_manifest=Path(args.input_manifest).resolve(),
        tables_dir=Path(args.tables_dir).resolve(),
        plots_dir=Path(args.plots_dir).resolve(),
        include_exploratory=args.include_exploratory,
        joint_components=args.joint_components,
        random_seed=args.random_seed or manifest.get("random_seed", 42),
        target_modalities=args.modalities,
    )
    run_fusion(config)


if __name__ == "__main__":
    main()
