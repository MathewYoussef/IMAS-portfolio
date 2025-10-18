#!/usr/bin/env python3
"""
Extract advanced spectral descriptors for Phase 9.

The script consumes canonical reflectance spectra, computes peak/AUC/
continuum-removed/FWHM/curvature metrics for each configured wavelength
window, generates quick QA samples, and records provenance in a feature
manifest so downstream phases can ingest the enriched features.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

DEFAULT_MANIFEST = Path("analysis/phase9/inputs/phase9_input_manifest.json")


@dataclass
class BandWindow:
    label: str
    lower_nm: float
    upper_nm: float
    fixed_peak_nm: float | None = None
    notes: str | None = None


def load_manifest(path: Path) -> Dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract advanced spectral features for Phase 9.")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="Path to the Phase 9 manifest.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/phase9/tables/advanced_band_features.csv"),
        help="Destination CSV for advanced features.",
    )
    parser.add_argument(
        "--feature-manifest",
        type=Path,
        default=Path("analysis/phase9/tables/phase9_feature_manifest.json"),
        help="Manifest describing generated feature tables and parameters.",
    )
    return parser.parse_args()


def extract_mean_columns(df: pd.DataFrame) -> List[str]:
    mean_columns = [col for col in df.columns if col.startswith("mean_")]
    mean_columns.sort(key=lambda name: int(name.split("_")[1]))
    return mean_columns


def resolve_window_mask(
    wavelengths: np.ndarray,
    window: BandWindow,
) -> Tuple[np.ndarray, float, float, float, float]:
    lower = window.lower_nm
    upper = window.upper_nm
    if lower == 0.0 and upper == 0.0:
        lower = float(wavelengths[0])
        upper = float(wavelengths[-1])
    mask = (wavelengths >= lower) & (wavelengths <= upper)
    if not np.any(mask):
        return mask, lower, upper, float("nan"), float("nan")
    selected = wavelengths[mask]
    return mask, lower, upper, float(selected[0]), float(selected[-1])


def compute_band_metrics(
    band_values: np.ndarray,
    band_wavelengths: np.ndarray,
    fixed_peak_nm: float | None = None,
) -> Dict[str, float]:
    if band_values.size == 0:
        return {
            "peak_reflectance": float("nan"),
            "peak_wavelength_nm": float("nan"),
            "auc": float("nan"),
            "continuum_depth": float("nan"),
            "fwhm_nm": float("nan"),
            "max_curvature": float("nan"),
        }

    series = band_values.astype(float)
    spacing = float(np.mean(np.diff(band_wavelengths))) if band_wavelengths.size > 1 else float("nan")

    peak_index = int(np.argmax(series))
    peak_reflectance = float(series[peak_index])
    peak_wavelength = float(band_wavelengths[peak_index])
    auc = float(np.trapz(series, band_wavelengths))

    continuum = np.linspace(series[0], series[-1], len(series))
    continuum_removed = continuum - series
    depth = float(np.max(continuum_removed))
    if depth < 0:
        depth = 0.0

    if depth > 0 and len(continuum_removed) > 1:
        half_depth = depth / 2
        mask = continuum_removed >= half_depth
        if np.any(mask):
            masked_wavelengths = band_wavelengths[mask]
            if masked_wavelengths.size > 1:
                fwhm = float(masked_wavelengths[-1] - masked_wavelengths[0])
            else:
                fwhm = 0.0
        else:
            fwhm = 0.0
    else:
        fwhm = 0.0

    if band_wavelengths.size >= 3 and not np.isnan(spacing) and spacing > 0:
        second_derivative = np.gradient(np.gradient(series, spacing), spacing)
        max_curvature = float(np.max(np.abs(second_derivative)))
    else:
        max_curvature = float("nan")

    metrics = {
        "peak_reflectance": peak_reflectance,
        "peak_wavelength_nm": peak_wavelength,
        "auc": auc,
        "continuum_depth": depth,
        "fwhm_nm": fwhm,
        "max_curvature": max_curvature,
    }
    if fixed_peak_nm is not None and band_wavelengths.size > 0:
        idx = int(np.argmin(np.abs(band_wavelengths - fixed_peak_nm)))
        metrics["fixed_peak_wavelength_nm"] = float(band_wavelengths[idx])
        metrics["fixed_peak_reflectance"] = float(series[idx])
        metrics["fixed_continuum_depth"] = float(continuum_removed[idx])
    else:
        metrics["fixed_peak_wavelength_nm"] = float("nan")
        metrics["fixed_peak_reflectance"] = float("nan")
        metrics["fixed_continuum_depth"] = float("nan")
    return metrics


def dose_sort_key(series: pd.Series) -> pd.Series:
    extracted = series.str.extract(r"(\d+)").astype(float)
    return extracted[0].fillna(np.inf)


def process_dataframe(
    df: pd.DataFrame,
    kind_column: str,
    kind_mapping: Dict[str, str],
    wavelengths: np.ndarray,
    mean_columns: Sequence[str],
    band_windows: Sequence[BandWindow],
    source: str,
) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    mean_matrix = df[mean_columns].to_numpy(dtype=float)

    for idx, row in df.iterrows():
        values = mean_matrix[idx]
        kind_raw = row[kind_column]
        kind = kind_mapping.get(kind_raw, kind_raw)
        for window in band_windows:
            mask, lower, upper, actual_lower, actual_upper = resolve_window_mask(wavelengths, window)
            if not np.any(mask):
                continue

            band_values = values[mask]
            band_wavelengths = wavelengths[mask]
            metrics = compute_band_metrics(band_values, band_wavelengths, window.fixed_peak_nm)

            record: Dict[str, Any] = {
                "dose_id": row["dose_id"],
                "kind": kind,
                "source": source,
                "band_label": window.label,
                "band_lower_nm": lower,
                "band_upper_nm": upper,
                "actual_band_lower_nm": actual_lower,
                "actual_band_upper_nm": actual_upper,
                "uva_mw_cm2": float(row.get("uva_mw_cm2", np.nan)),
                "uvb_mw_cm2": float(row.get("uvb_mw_cm2", np.nan)),
                "trim_fraction": float(row.get("trim_fraction", np.nan)),
            }
            for sample_field in ("samples_total", "samples_used", "samples_12", "samples_6"):
                if sample_field in row:
                    record[sample_field] = float(row[sample_field])
            record.update(metrics)
            records.append(record)
    return records


def generate_qa_samples(
    stats_df: pd.DataFrame,
    wavelengths: np.ndarray,
    mean_columns: Sequence[str],
    band_windows: Sequence[BandWindow],
    qa_checks: Dict[str, Any],
) -> List[str]:
    saved: List[str] = []
    profile_path = qa_checks.get("profile_samples")
    continuum_path = qa_checks.get("continuum_samples")

    if not profile_path and not continuum_path:
        return saved

    primary_window = next((window for window in band_windows if window.label != "peak_band_auto"), band_windows[0])
    mask, _, _, _, _ = resolve_window_mask(wavelengths, primary_window)
    if not np.any(mask):
        return saved

    subset = stats_df.loc[stats_df["angle"].isin(["12Oclock", "12oClock", "12OClock"])]
    if subset.empty:
        subset = stats_df
    if subset.empty:
        return saved

    subset = subset.assign(_dose_order=dose_sort_key(subset["dose_id"])).sort_values("_dose_order", kind="mergesort")
    subset = subset.drop(columns="_dose_order")
    band_wavelengths = wavelengths[mask]

    if profile_path:
        records: List[Dict[str, Any]] = []
        for _, row in subset.iterrows():
            values = row[mean_columns].to_numpy(dtype=float)[mask]
            for wavelength, reflectance in zip(band_wavelengths, values.astype(float)):
                records.append(
                    {
                        "dose_id": row["dose_id"],
                        "angle": row["angle"],
                        "band_label": primary_window.label,
                        "wavelength_nm": float(wavelength),
                        "reflectance": float(reflectance),
                    }
                )
        profile_df = pd.DataFrame(records)
        path_profiles = Path(profile_path)
        path_profiles.parent.mkdir(parents=True, exist_ok=True)
        profile_df.to_csv(path_profiles, index=False)
        saved.append(str(path_profiles))

    if continuum_path:
        representative = subset.iloc[min(len(subset) - 1, 2)]
        values = representative[mean_columns].to_numpy(dtype=float)[mask]
        continuum = np.linspace(values[0], values[-1], len(values))
        continuum_removed = continuum - values
        cont_df = pd.DataFrame(
            {
                "dose_id": representative["dose_id"],
                "angle": representative["angle"],
                "band_label": primary_window.label,
                "wavelength_nm": band_wavelengths.astype(float),
                "continuum_minus_reflectance": continuum_removed.astype(float),
            }
        )
        path_cont = Path(continuum_path)
        path_cont.parent.mkdir(parents=True, exist_ok=True)
        cont_df.to_csv(path_cont, index=False)
        saved.append(str(path_cont))

    return saved


def write_feature_manifest(
    manifest_path: Path,
    source_manifest: Path,
    features_df: pd.DataFrame,
    processing_defaults: Dict[str, Any],
    qa_artifacts: Sequence[str],
) -> None:
    numeric_columns = [
        "peak_reflectance",
        "peak_wavelength_nm",
        "auc",
        "continuum_depth",
        "fwhm_nm",
        "max_curvature",
        "fixed_peak_reflectance",
        "fixed_peak_wavelength_nm",
        "fixed_continuum_depth",
    ]
    summary_stats: Dict[str, Dict[str, float]] = {}
    for column in numeric_columns:
        if column in features_df:
            values = pd.to_numeric(features_df[column], errors="coerce")
            summary_stats[column] = {
                "min": float(values.min()) if not values.isna().all() else float("nan"),
                "max": float(values.max()) if not values.isna().all() else float("nan"),
                "median": float(values.median()) if not values.isna().all() else float("nan"),
            }

    manifest_payload = {
        "source_manifest": str(source_manifest),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "record_count": int(len(features_df)),
        "bands": sorted(features_df["band_label"].dropna().unique().tolist()),
        "kinds": sorted(features_df["kind"].dropna().unique().tolist()),
        "processing": processing_defaults,
        "qa_artifacts": list(qa_artifacts),
        "metrics_summary": summary_stats,
        "status": "complete",
        "notes": (
            "Advanced spectral descriptors computed with endpoint continuum removal applied to denoised spectra "
            "(no additional smoothing). Delta composites are difference spectra, so peak/area metrics may be negative."
        ),
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest_payload, indent=2))


def main() -> None:
    args = parse_args()
    manifest = load_manifest(args.manifest)

    spectra_stats = Path(manifest["canonical_spectra"]["reflectance_stats"])
    composites_path = Path(manifest["canonical_spectra"]["reflectance_composites"])
    wavelength_grid = Path(manifest["canonical_spectra"]["wavelength_grid"])

    for path in (spectra_stats, composites_path, wavelength_grid):
        if not path.exists():
            raise FileNotFoundError(f"Required input missing: {path}")

    stats_df = pd.read_csv(spectra_stats)
    composites_df = pd.read_csv(composites_path)
    wavelengths = np.load(wavelength_grid)
    mean_columns = extract_mean_columns(stats_df)

    band_windows = [BandWindow(**entry) for entry in manifest.get("feature_windows", [])]
    if not band_windows:
        raise ValueError("No feature windows defined in manifest.")

    processing_defaults = manifest.get("processing_defaults", {})
    continuum_method = processing_defaults.get("continuum_method", "linear_hull")
    derivative_spacing = processing_defaults.get("derivative_spacing_nm", 0.5)

    angle_records = process_dataframe(
        stats_df,
        kind_column="angle",
        kind_mapping={"12Oclock": "12Oclock", "6Oclock": "6Oclock"},
        wavelengths=wavelengths,
        mean_columns=mean_columns,
        band_windows=band_windows,
        source="angle",
    )

    composite_mean_columns = extract_mean_columns(composites_df)
    composite_records = process_dataframe(
        composites_df,
        kind_column="composite",
        kind_mapping={"Sigma": "Sigma", "Delta": "Delta"},
        wavelengths=wavelengths,
        mean_columns=composite_mean_columns,
        band_windows=band_windows,
        source="composite",
    )

    records = angle_records + composite_records
    features_df = pd.DataFrame(records)
    if features_df.empty:
        raise RuntimeError("No features were computed; check band window configuration.")

    features_df = features_df.sort_values(["band_label", "kind", "dose_id"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    features_df.to_csv(args.output, index=False)

    qa_checks = manifest.get("qa_checks", {})
    qa_artifacts = generate_qa_samples(
        stats_df=stats_df,
        wavelengths=wavelengths,
        mean_columns=mean_columns,
        band_windows=band_windows,
        qa_checks=qa_checks,
    )

    write_feature_manifest(
        manifest_path=args.feature_manifest,
        source_manifest=args.manifest,
        features_df=features_df,
        processing_defaults={
            "smoothing": "none",
            "continuum_method": continuum_method,
            "derivative_spacing_nm": derivative_spacing,
        },
        qa_artifacts=qa_artifacts,
    )

    print(f"[INFO] Advanced features written to {args.output}")
    if qa_artifacts:
        print(f"[INFO] QA samples written: {', '.join(qa_artifacts)}")
    print(f"[INFO] Feature manifest written to {args.feature_manifest}")


if __name__ == "__main__":
    main()
