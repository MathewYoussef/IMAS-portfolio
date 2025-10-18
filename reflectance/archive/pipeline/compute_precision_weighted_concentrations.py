#!/usr/bin/env python3
"""Compute precision-weighted latent scytonemin concentrations per sample and treatment.

The workflow:
1. Fit a Deming regression (λ = 1) between chromatogram and DAD concentrations for each state.
2. Calibrate DAD readings into the chromatogram scale using the fitted intercept/slope.
3. Estimate measurement variance from the residuals and fuse the two observations via inverse-variance weighting.
4. Report per-sample latent concentrations with standard errors and treatment-level trimmed summaries.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

from dose_metadata import attach_dose_metadata, iter_dose_records


SCYTONEMIN_STATES = ("total", "oxidized", "reduced")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate precision-weighted latent concentrations from chromatogram and DAD inputs."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("Combined_Scytonemin_Concentrations.csv"),
        help="CSV containing chromatogram and DAD concentrations per sample.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("aggregated_reflectance"),
        help="Directory where precision-weighted outputs will be written.",
    )
    parser.add_argument(
        "--trim-fraction",
        type=float,
        default=0.2,
        help="Trim fraction used for treatment-level aggregation.",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=1e-6,
        help="Floor added to variance estimates to avoid division by zero.",
    )
    return parser.parse_args()


def deming_fit(x: np.ndarray, y: np.ndarray, lambda_ratio: float = 1.0) -> Tuple[float, float]:
    """Closed-form Deming regression with error-variance ratio lambda."""
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    x_c = x - x_mean
    y_c = y - y_mean

    s_xx = np.mean(x_c**2)
    s_yy = np.mean(y_c**2)
    s_xy = np.mean(x_c * y_c)

    if np.isclose(s_xy, 0.0):
        raise ValueError("Covariance is zero; cannot fit Deming regression.")

    term = s_yy - lambda_ratio * s_xx
    slope = (term + np.sqrt(term**2 + 4 * lambda_ratio * s_xy**2)) / (2 * s_xy)
    intercept = y_mean - slope * x_mean
    return float(intercept), float(slope)


def trimmed_mean(values: np.ndarray, trim_fraction: float) -> Tuple[float, float, int]:
    """Return trimmed mean and std along axis 0 with symmetric trimming."""
    values = np.asarray(values, dtype=float)
    n = values.size
    if n == 0:
        return np.nan, np.nan, 0
    if trim_fraction < 0 or trim_fraction >= 0.5:
        raise ValueError("trim_fraction must be in [0, 0.5).")
    if trim_fraction == 0 or n <= 2:
        trimmed = np.sort(values)
    else:
        cut = int(np.floor(trim_fraction * n))
        if cut * 2 >= n:
            cut = max(0, (n - 1) // 2)
        sorted_vals = np.sort(values)
        trimmed = sorted_vals[cut : n - cut]
        if trimmed.size == 0:
            trimmed = sorted_vals
    mean = float(trimmed.mean())
    std = float(trimmed.std(ddof=1)) if trimmed.size > 1 else np.nan
    return mean, std, trimmed.size


def compute_latent_concentrations(
    df: pd.DataFrame, epsilon: float
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-sample latent concentrations and model summary stats."""
    base = df[["sample_id"]].copy()
    base["treatment"] = df["sample_id"].str[0].map(lambda t: f"treatment_{t}")
    base["dose_id"] = df["sample_id"].str[0].astype(int).map(lambda x: f"dose_{7 - x}")
    base["uva_mw_cm2"] = base["dose_id"].map(lambda label: attach_dose_metadata(label).uva_mw_cm2)
    base["uvb_mw_cm2"] = base["dose_id"].map(lambda label: attach_dose_metadata(label).uvb_mw_cm2)

    summaries: Dict[str, Dict[str, float]] = {}

    results = base.copy()

    for state in SCYTONEMIN_STATES:
        chrom_col = f"chrom_{state}_mg_per_gDW"
        dad_col = f"dad_{state}_mg_per_gDW"

        chrom = df[chrom_col].to_numpy(dtype=float)
        dad = df[dad_col].to_numpy(dtype=float)
        mask = np.isfinite(chrom) & np.isfinite(dad)

        if mask.sum() < 3:
            raise ValueError(f"Not enough valid pairs to fit Deming regression for {state}.")

        chrom_valid = chrom[mask]
        dad_valid = dad[mask]

        intercept, slope = deming_fit(chrom_valid, dad_valid, lambda_ratio=1.0)
        dad_calibrated = (dad - intercept) / slope

        residual_array = chrom - dad_calibrated
        finite_residuals = residual_array[np.isfinite(residual_array)]
        resid_var = float(np.var(finite_residuals, ddof=1)) if finite_residuals.size > 1 else 0.0
        var_each = max(resid_var / 2.0, epsilon)

        chrom_var = np.full(chrom.shape, var_each)
        dad_var = np.full(dad.shape, var_each)

        fused = (chrom / chrom_var + dad_calibrated / dad_var) / (1.0 / chrom_var + 1.0 / dad_var)
        fused_se = np.sqrt(1.0 / (1.0 / chrom_var + 1.0 / dad_var))

        results[f"{state}_chrom"] = chrom
        results[f"{state}_dad"] = dad
        results[f"{state}_dad_calibrated"] = dad_calibrated
        results[f"{state}_latent"] = fused
        results[f"{state}_latent_se"] = fused_se
        results[f"{state}_residual"] = residual_array

        summaries[state] = {
            "intercept": intercept,
            "slope": slope,
            "residual_variance": resid_var,
            "assumed_measurement_variance": var_each,
        }

    summary_df = pd.DataFrame.from_dict(summaries, orient="index").reset_index().rename(columns={"index": "state"})
    return results, summary_df


def aggregate_by_treatment(
    sample_df: pd.DataFrame, trim_fraction: float
) -> pd.DataFrame:
    records = []
    for (treatment, dose_id), group in sample_df.groupby(["treatment", "dose_id"]):
        dose_meta = attach_dose_metadata(dose_id)
        record = {
            "treatment": treatment,
            "dose_id": dose_id,
            "trim_fraction": trim_fraction,
            "sample_count": len(group),
            "uva_mw_cm2": dose_meta.uva_mw_cm2,
            "uvb_mw_cm2": dose_meta.uvb_mw_cm2,
        }
        for state in SCYTONEMIN_STATES:
            latent_col = f"{state}_latent"
            se_col = f"{state}_latent_se"
            mean, std, used = trimmed_mean(group[latent_col].dropna().to_numpy(), trim_fraction)
            record[f"{state}_latent_mean_trimmed"] = mean
            record[f"{state}_latent_std_trimmed"] = std
            record[f"{state}_latent_n_used"] = used
            record[f"{state}_latent_se_median"] = float(np.median(group[se_col].dropna().to_numpy())) if group[se_col].notna().any() else np.nan
        records.append(record)
    return pd.DataFrame(records).sort_values("treatment")


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    if "sample_id" not in df.columns:
        raise ValueError("Input CSV must contain 'sample_id'.")

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    sample_latent_df, summary_df = compute_latent_concentrations(df, args.epsilon)
    treatment_latent_df = aggregate_by_treatment(sample_latent_df, args.trim_fraction)

    sample_path = outdir / "precision_weighted_concentrations.csv"
    treatment_path = outdir / "precision_weighted_concentrations_treatment.csv"
    summary_path = outdir / "precision_weighted_concentration_model_stats.csv"
    meta_path = outdir / "precision_weighted_concentration_metadata.json"

    sample_latent_df.to_csv(sample_path, index=False)
    treatment_latent_df.to_csv(treatment_path, index=False)
    summary_df.to_csv(summary_path, index=False)

    meta = {
        "states": SCYTONEMIN_STATES,
        "trim_fraction": args.trim_fraction,
        "epsilon": args.epsilon,
        "files": {
            "sample": str(sample_path),
            "treatment": str(treatment_path),
            "model_stats": str(summary_path),
        },
        "uva_uvb_by_dose": {
            record.label: {
                "uva_mw_cm2": record.uva_mw_cm2,
                "uvb_mw_cm2": record.uvb_mw_cm2,
            }
            for record in iter_dose_records()
        },
    }
    meta_path.write_text(json.dumps(meta, indent=2))

    print(f"Wrote {sample_path} with {len(sample_latent_df)} samples.")
    print(f"Wrote {treatment_path} with {len(treatment_latent_df)} treatments.")
    print(f"Model summary saved to {summary_path}.")


if __name__ == "__main__":
    main()
