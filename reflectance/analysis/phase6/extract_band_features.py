#!/usr/bin/env python3
"""
Extract band-level reflectance features from the canonical dose summaries.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


@dataclass
class BandWindow:
    label: str
    lower_nm: float
    upper_nm: float


@dataclass
class BandConfig:
    manifest_path: Path
    reflectance_stats: Path
    reflectance_composites: Path
    wavelength_grid: Path
    band_windows: List[BandWindow]
    output_path: Path


def load_band_windows(entries: List[Dict[str, float]]) -> List[BandWindow]:
    return [BandWindow(**entry) for entry in entries]


def load_wavelength_grid(path: Path) -> np.ndarray:
    return np.load(path)


def extract_mean_columns(df: pd.DataFrame) -> List[str]:
    return sorted(col for col in df.columns if col.startswith("mean_"))


def compute_band_metrics(values: np.ndarray, wavelengths: np.ndarray) -> Dict[str, float]:
    metrics = {
        "band_mean": float(np.mean(values)),
        "band_median": float(np.median(values)),
        "band_min": float(np.min(values)),
        "band_max": float(np.max(values)),
        "band_range": float(np.max(values) - np.min(values)),
        "band_area": float(np.trapz(values, wavelengths)),
    }
    return metrics


def process_dataframe(
    df: pd.DataFrame,
    kind_column: str,
    kind_mapping: Dict[str, str],
    wavelengths: np.ndarray,
    band_windows: List[BandWindow],
    source: str,
) -> List[Dict[str, object]]:
    mean_columns = extract_mean_columns(df)
    mean_values = df[mean_columns].to_numpy(dtype=float)
    records: List[Dict[str, object]] = []

    for index, row in df.iterrows():
        kind_raw = row[kind_column]
        kind = kind_mapping.get(kind_raw, kind_raw)
        dose_id = row["dose_id"]
        uva = row["uva_mw_cm2"]
        uvb = row["uvb_mw_cm2"]
        values = mean_values[index]

        for band in band_windows:
            mask = (wavelengths >= band.lower_nm) & (wavelengths <= band.upper_nm)
            if not mask.any():
                continue
            band_values = values[mask]
            band_wavelengths = wavelengths[mask]
            metrics = compute_band_metrics(band_values, band_wavelengths)
            record = {
                "dose_id": dose_id,
                "kind": kind,
                "band_label": band.label,
                "band_lower_nm": band.lower_nm,
                "band_upper_nm": band.upper_nm,
                "uva_mw_cm2": uva,
                "uvb_mw_cm2": uvb,
                "source": source,
            }
            record.update(metrics)
            records.append(record)
    return records


def extract_band_features(config: BandConfig) -> None:
    manifest = json.loads(config.manifest_path.read_text())
    wavelengths = load_wavelength_grid(config.wavelength_grid)
    band_windows = load_band_windows(manifest["band_windows"])

    stats_df = pd.read_csv(config.reflectance_stats)
    composites_df = pd.read_csv(config.reflectance_composites)

    angle_records = process_dataframe(
        stats_df,
        kind_column="angle",
        kind_mapping={"12Oclock": "12Oclock", "6Oclock": "6Oclock"},
        wavelengths=wavelengths,
        band_windows=band_windows,
        source="angle",
    )

    composite_records = process_dataframe(
        composites_df,
        kind_column="composite",
        kind_mapping={"Sigma": "Sigma", "Delta": "Delta"},
        wavelengths=wavelengths,
        band_windows=band_windows,
        source="composite",
    )

    records = angle_records + composite_records
    df = pd.DataFrame(records)
    config.output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(config.output_path, index=False)
    print(f"[INFO] Band features written to {config.output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract reflectance band features for Phase 6 baseline models.")
    parser.add_argument(
        "--manifest",
        default="analysis/phase6/inputs/phase6_input_manifest.json",
        help="Phase 6 input manifest.",
    )
    parser.add_argument(
        "--output",
        default="analysis/phase6/tables/band_features.csv",
        help="Output CSV path for band features.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = json.loads(Path(args.manifest).read_text())
    config = BandConfig(
        manifest_path=Path(args.manifest).resolve(),
        reflectance_stats=Path(manifest["canonical_reflectance_stats"]).resolve(),
        reflectance_composites=Path(manifest["canonical_reflectance_composites"]).resolve(),
        wavelength_grid=Path(manifest["wavelength_grid"]).resolve(),
        band_windows=load_band_windows(manifest["band_windows"]),
        output_path=Path(args.output).resolve(),
    )
    extract_band_features(config)


if __name__ == "__main__":
    main()
