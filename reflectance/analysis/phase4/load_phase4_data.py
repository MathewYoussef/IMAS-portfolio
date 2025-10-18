#!/usr/bin/env python3
"""
Phase 4 data loader.

Read the Phase 3 join outputs, assemble z-scored reflectance and concentration
matrices per wavelength window and reflectance kind, and cache the results for
downstream congruence analyses (PCA/Procrustes/PLS/dCor/HSIC).
"""

from __future__ import annotations

import argparse
import json
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

REFLECTANCE_FEATURES = [
    "mean_reflectance",
    "median_reflectance",
    "min_reflectance",
    "max_reflectance",
    "range_reflectance",
    "area_reflectance",
    "mean_std",
    "median_std",
    "mean_mad",
    "median_mad",
]

CONCENTRATION_COLUMNS = {
    "chrom": ["chrom_total_mean", "chrom_oxidized_mean", "chrom_reduced_mean"],
    "dad": ["dad_total_mean", "dad_oxidized_mean", "dad_reduced_mean"],
    "latent": ["latent_total_mean", "latent_oxidized_mean", "latent_reduced_mean"],
}

REFLECTANCE_KINDS = ["Sigma", "Delta", "12Oclock", "6Oclock"]
DOSE_ORDER = ["dose_1", "dose_2", "dose_3", "dose_4", "dose_5", "dose_6"]


@dataclass
class LoaderConfig:
    phase4_manifest: Path
    cache_path: Path
    metadata_path: Path


def load_json(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text())


def order_by_dose(df: pd.DataFrame) -> pd.DataFrame:
    order_map = {label: idx for idx, label in enumerate(DOSE_ORDER)}
    return df.sort_values("dose_id", key=lambda s: s.map(order_map))


def zscore_matrix(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[int]]:
    means = matrix.mean(axis=0)
    stds = matrix.std(axis=0, ddof=0)
    zero_mask = np.isclose(stds, 0.0) | np.isnan(stds)
    stds_adj = stds.copy()
    stds_adj[zero_mask] = 1.0
    matrix_z = (matrix - means) / stds_adj
    matrix_z[:, zero_mask] = 0.0
    zero_indices = np.where(zero_mask)[0].tolist()
    return matrix_z, means, stds, zero_indices


def extract_reflectance_matrices(df: pd.DataFrame) -> Dict[str, Dict[str, object]]:
    matrices: Dict[str, Dict[str, object]] = {}
    for kind in REFLECTANCE_KINDS:
        subset = df[df["kind"] == kind].copy()
        if subset.empty:
            continue
        subset = order_by_dose(subset)
        feature_cols = [col for col in REFLECTANCE_FEATURES if col in subset.columns]
        matrix = subset[feature_cols].to_numpy(dtype=float)
        matrix_z, means, stds, zero_cols = zscore_matrix(matrix)
        matrices[kind] = {
            "matrix": matrix_z,
            "columns": feature_cols,
            "scaling": {
                "mean": means.tolist(),
                "std": stds.tolist(),
                "zero_variance_columns": [feature_cols[idx] for idx in zero_cols],
            },
        }
    return matrices


def extract_concentration_matrices(df: pd.DataFrame, modality: str) -> Dict[str, Dict[str, object]]:
    matrices: Dict[str, Dict[str, object]] = {}
    columns = [col for col in CONCENTRATION_COLUMNS[modality] if col in df.columns]
    if not columns:
        return matrices
    subset = order_by_dose(df[["dose_id"] + columns].drop_duplicates())
    matrix = subset[columns].to_numpy(dtype=float)
    matrix_z, means, stds, zero_cols = zscore_matrix(matrix)
    matrices[modality] = {
        "matrix": matrix_z,
        "columns": columns,
        "scaling": {
            "mean": means.tolist(),
            "std": stds.tolist(),
            "zero_variance_columns": [columns[idx] for idx in zero_cols],
        },
    }
    return matrices


def extract_dose_metadata(df: pd.DataFrame) -> List[Dict[str, float]]:
    subset = order_by_dose(df[["dose_id", "uva_mw_cm2", "uvb_mw_cm2"]].drop_duplicates())
    records = []
    for _, row in subset.iterrows():
        records.append(
            {
                "dose_id": row["dose_id"],
                "uva_mw_cm2": float(row["uva_mw_cm2"]),
                "uvb_mw_cm2": float(row["uvb_mw_cm2"]),
            }
        )
    return records


def build_loader(config: LoaderConfig) -> None:
    phase4_manifest = load_json(config.phase4_manifest)
    join_manifest = load_json(Path(phase4_manifest["phase3_join_manifest"]).resolve())

    windows_data: Dict[str, Dict[str, object]] = {}
    outputs = join_manifest.get("outputs", [])

    for entry in outputs:
        window = entry["window"]
        window_label = f"{int(window['lower_nm'])}_{int(window['upper_nm'])}nm"

        chrom_path = Path(entry["chrom_join"]).resolve()
        dad_path = Path(entry["dad_join"]).resolve()
        latent_path = Path(entry["latent_join"]).resolve()

        if not (chrom_path.exists() and dad_path.exists() and latent_path.exists()):
            print(f"[WARN] Missing join files for window {window_label}; skipping.")
            continue

        chrom_df = order_by_dose(pd.read_csv(chrom_path))
        dad_df = order_by_dose(pd.read_csv(dad_path))
        latent_df = order_by_dose(pd.read_csv(latent_path))

        reflectance_matrices = extract_reflectance_matrices(chrom_df)

        concentration_matrices: Dict[str, Dict[str, object]] = {}
        concentration_matrices.update(extract_concentration_matrices(chrom_df, "chrom"))
        concentration_matrices.update(extract_concentration_matrices(dad_df, "dad"))
        concentration_matrices.update(extract_concentration_matrices(latent_df, "latent"))

        dose_metadata = extract_dose_metadata(chrom_df)

        windows_data[window_label] = {
            "reflectance": reflectance_matrices,
            "concentration": concentration_matrices,
            "dose_metadata": dose_metadata,
        }

    cache_payload = {
        "windows": windows_data,
        "settings": {
            "permutations_default": phase4_manifest.get("default_permutations"),
            "random_seed_default": phase4_manifest.get("random_seed"),
        },
        "sources": {
            "phase4_manifest": str(config.phase4_manifest),
            "phase3_join_manifest": phase4_manifest["phase3_join_manifest"],
            "phase3_diagnostics_manifest": phase4_manifest["diagnostics_manifest"],
        },
    }

    with config.cache_path.open("wb") as handle:
        pickle.dump(cache_payload, handle)

    metadata = {
        "cache_path": str(config.cache_path),
        "settings": cache_payload["settings"],
        "windows": {},
    }

    for window_label, payload in windows_data.items():
        window_meta = {
            "dose_order": [record["dose_id"] for record in payload["dose_metadata"]],
            "dose_metadata": payload["dose_metadata"],
            "reflectance": {},
            "concentration": {},
        }
        for kind, info in payload["reflectance"].items():
            window_meta["reflectance"][kind] = {
                "columns": info["columns"],
                "shape": [len(payload["dose_metadata"]), len(info["columns"])],
                "zero_variance_columns": info["scaling"]["zero_variance_columns"],
            }
        for modality, info in payload["concentration"].items():
            window_meta["concentration"][modality] = {
                "columns": info["columns"],
                "shape": [len(payload["dose_metadata"]), len(info["columns"])],
                "zero_variance_columns": info["scaling"]["zero_variance_columns"],
            }
        metadata["windows"][window_label] = window_meta

    config.metadata_path.write_text(json.dumps(metadata, indent=2))

    print(f"[INFO] Cached Phase 4 matrices to {config.cache_path}")
    print(f"[INFO] Metadata written to {config.metadata_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assemble z-scored matrices for Phase 4 congruence analyses.")
    parser.add_argument(
        "--manifest",
        default="analysis/phase4/inputs/phase4_input_manifest.json",
        help="Phase 4 input manifest referencing Phase 3 outputs.",
    )
    parser.add_argument(
        "--cache-path",
        default="analysis/phase4/inputs/phase4_data_cache.pkl",
        help="Path to write the cached matrices pickle.",
    )
    parser.add_argument(
        "--metadata-path",
        default="analysis/phase4/inputs/phase4_matrix_manifest.json",
        help="Path to write a JSON summary of matrices and scaling.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = LoaderConfig(
        phase4_manifest=Path(args.manifest).resolve(),
        cache_path=Path(args.cache_path).resolve(),
        metadata_path=Path(args.metadata_path).resolve(),
    )
    build_loader(config)


if __name__ == "__main__":
    main()
