#!/usr/bin/env python3
"""
Phase 5 block loader.

Assemble multiblock structures for reflectance (per-angle, composites optional)
and concentration modalities based on the Phase 4 z-scored cache.
"""

from __future__ import annotations

import argparse
import json
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


@dataclass
class BlockConfig:
    input_manifest: Path
    phase4_cache: Path
    phase4_matrix_manifest: Path
    block_cache: Path
    block_manifest: Path


def load_json(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text())


def load_cache(path: Path) -> Dict[str, object]:
    with path.open("rb") as handle:
        return pickle.load(handle)


def filter_zero_variance(info: Dict[str, object]) -> Tuple[np.ndarray, List[str], List[str]]:
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


def assemble_blocks(config: BlockConfig) -> None:
    manifest = load_json(config.input_manifest)
    cache = load_cache(config.phase4_cache)
    matrix_manifest = load_json(config.phase4_matrix_manifest)

    primary_blocks = manifest["reflectance_blocks"]["primary"]
    exploratory_blocks = manifest["reflectance_blocks"]["exploratory"]

    block_payload: Dict[str, Dict[str, object]] = {}

    for window, payload in cache["windows"].items():
        dose_metadata = payload["dose_metadata"]
        reflectance = payload["reflectance"]
        concentration = payload["concentration"]

        window_entry: Dict[str, object] = {
            "dose_metadata": dose_metadata,
            "reflectance": {},
            "concentration": {},
            "notes": {
                "primary_blocks": primary_blocks,
                "exploratory_blocks": exploratory_blocks,
            },
        }

        for block_name in primary_blocks + exploratory_blocks:
            if block_name not in reflectance:
                continue
            matrix, kept_cols, dropped_cols = filter_zero_variance(reflectance[block_name])
            window_entry["reflectance"][block_name] = {
                "matrix": matrix,
                "columns": kept_cols,
                "dropped_columns": dropped_cols,
                "shape": matrix.shape,
            }

        for modality, info in concentration.items():
            matrix, kept_cols, dropped_cols = filter_zero_variance(info)
            window_entry["concentration"][modality] = {
                "matrix": matrix,
                "columns": kept_cols,
                "dropped_columns": dropped_cols,
                "shape": matrix.shape,
            }

        block_payload[window] = window_entry

    payload = {
        "blocks": block_payload,
        "settings": {
            "default_components": manifest.get("default_components", {}),
            "default_permutations": manifest.get("default_permutations"),
            "random_seed": manifest.get("random_seed"),
        },
        "sources": {
            "phase4_cache": str(config.phase4_cache),
            "phase4_matrix_manifest": str(config.phase4_matrix_manifest),
            "phase4_congruence_manifest": manifest.get("phase4_congruence_manifest"),
            "phase4_dependence_manifest": manifest.get("phase4_dependence_manifest"),
        },
    }

    with config.block_cache.open("wb") as handle:
        pickle.dump(payload, handle)

    summary = {
        "block_cache": str(config.block_cache),
        "settings": payload["settings"],
        "windows": {},
    }

    for window, entry in block_payload.items():
        summary["windows"][window] = {
            "dose_order": [meta["dose_id"] for meta in entry["dose_metadata"]],
            "reflectance_blocks": {
                name: {
                    "columns": info["columns"],
                    "dropped_columns": info["dropped_columns"],
                    "shape": info["shape"],
                }
                for name, info in entry["reflectance"].items()
            },
            "concentration_blocks": {
                name: {
                    "columns": info["columns"],
                    "dropped_columns": info["dropped_columns"],
                    "shape": info["shape"],
                }
                for name, info in entry["concentration"].items()
            },
        }

    config.block_manifest.write_text(json.dumps(summary, indent=2))
    print(f"[INFO] Block cache written to {config.block_cache}")
    print(f"[INFO] Block manifest written to {config.block_manifest}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assemble Phase 5 multiblock structures from Phase 4 cache.")
    parser.add_argument(
        "--manifest",
        default="analysis/phase5/inputs/phase5_input_manifest.json",
        help="Phase 5 manifest referencing Phase 4 cache.",
    )
    parser.add_argument(
        "--block-cache",
        default="analysis/phase5/inputs/phase5_block_cache.pkl",
        help="Output pickle for assembled blocks.",
    )
    parser.add_argument(
        "--block-manifest",
        default="analysis/phase5/inputs/phase5_block_manifest.json",
        help="Output JSON summarising block metadata.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = load_json(Path(args.manifest).resolve())
    config = BlockConfig(
        input_manifest=Path(args.manifest).resolve(),
        phase4_cache=Path(manifest["phase4_cache"]).resolve(),
        phase4_matrix_manifest=Path(manifest["phase4_matrix_manifest"]).resolve(),
        block_cache=Path(args.block_cache).resolve(),
        block_manifest=Path(args.block_manifest).resolve(),
    )
    assemble_blocks(config)


if __name__ == "__main__":
    main()
