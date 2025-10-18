#!/usr/bin/env python3
"""
Analyze chromatogram vs DAD dry-weight-normalized scytonemin concentrations.

For each scytonemin state (total, oxidized, reduced) the script:
  * Computes paired statistics between chromatogram and DAD estimates.
  * Runs a paired t-test and simple linear regression to check for scale/offset differences.
  * Generates scatter plots with identity lines and sample-wise line plots.
  * Stores a CSV summary table with key metrics for downstream reporting.
"""

import argparse
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats


def _configure_matplotlib_cache() -> None:
    """Ensure matplotlib can write its cache within the repository."""
    if "MPLCONFIGDIR" in os.environ:
        return

    script_dir = Path(__file__).resolve().parent
    candidates = [
        script_dir / ".matplotlib-cache",
        script_dir / ".cache",
        script_dir,
    ]

    for candidate in candidates:
        try:
            candidate.mkdir(exist_ok=True)
        except Exception:
            # Directory may already exist or be read-only; ignore and try next.
            pass
        if candidate.exists() and os.access(candidate, os.W_OK):
            os.environ["MPLCONFIGDIR"] = str(candidate)
            return

    # Fall back to current working directory (may still be read-only, but it's the best effort).
    cwd = Path.cwd()
    if os.access(cwd, os.W_OK):
        os.environ["MPLCONFIGDIR"] = str(cwd)


_configure_matplotlib_cache()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


STATES = ("total", "oxidized", "reduced")


@dataclass
class ComparisonStats:
    state: str
    sample_count: int
    chrom_mean: float
    dad_mean: float
    mean_difference: float
    std_difference: float
    mean_ratio: float
    std_ratio: float
    t_statistic: float
    t_pvalue: float
    slope: float
    intercept: float
    r_value: float
    regression_pvalue: float


def load_data(csv_path: Path) -> pd.DataFrame:
    """Load the combined scytonemin dataset."""
    df = pd.read_csv(csv_path)
    required_cols = ["sample_id"]
    for state in STATES:
        required_cols.extend(
            [
                f"chrom_{state}_mg_per_gDW",
                f"dad_{state}_mg_per_gDW",
            ]
        )

    missing = sorted(set(required_cols) - set(df.columns))
    if missing:
        missing_str = ", ".join(missing)
        raise ValueError(f"Missing expected columns: {missing_str}")

    return df


def compute_state_statistics(df: pd.DataFrame, state: str) -> ComparisonStats:
    """Calculate paired statistics between chromatogram and DAD measures."""
    chrom_col = f"chrom_{state}_mg_per_gDW"
    dad_col = f"dad_{state}_mg_per_gDW"

    subset = df[["sample_id", chrom_col, dad_col]].dropna()
    chrom_values = subset[chrom_col].to_numpy()
    dad_values = subset[dad_col].to_numpy()

    if chrom_values.size == 0:
        raise ValueError(f"No overlapping data for state '{state}'.")

    differences = dad_values - chrom_values
    ratios = []
    for chrom_val, dad_val in zip(chrom_values, dad_values):
        if chrom_val != 0:
            ratios.append(dad_val / chrom_val)

    mean_ratio = float(np.mean(ratios)) if ratios else np.nan
    std_ratio = float(np.std(ratios, ddof=1)) if len(ratios) > 1 else np.nan

    if chrom_values.size > 1:
        t_res = stats.ttest_rel(dad_values, chrom_values, nan_policy="omit")
        t_stat = float(t_res.statistic)
        t_pvalue = float(t_res.pvalue)
        slope, intercept, r_value, p_value, _ = stats.linregress(
            chrom_values, dad_values
        )
    else:
        t_stat = np.nan
        t_pvalue = np.nan
        slope = np.nan
        intercept = np.nan
        r_value = np.nan
        p_value = np.nan

    return ComparisonStats(
        state=state,
        sample_count=int(chrom_values.size),
        chrom_mean=float(np.mean(chrom_values)),
        dad_mean=float(np.mean(dad_values)),
        mean_difference=float(np.mean(differences)),
        std_difference=float(np.std(differences, ddof=1))
        if differences.size > 1
        else np.nan,
        mean_ratio=mean_ratio,
        std_ratio=std_ratio,
        t_statistic=t_stat,
        t_pvalue=t_pvalue,
        slope=float(slope) if not np.isnan(slope) else np.nan,
        intercept=float(intercept) if not np.isnan(intercept) else np.nan,
        r_value=float(r_value) if not np.isnan(r_value) else np.nan,
        regression_pvalue=float(p_value) if not np.isnan(p_value) else np.nan,
    )


def plot_scatter(df: pd.DataFrame, state: str, outdir: Path) -> Path:
    """Create a scatter plot comparing chromatogram vs DAD concentrations."""
    chrom_col = f"chrom_{state}_mg_per_gDW"
    dad_col = f"dad_{state}_mg_per_gDW"

    subset = df[["sample_id", chrom_col, dad_col]].dropna()
    if subset.empty:
        raise ValueError(f"No data to plot for state '{state}'.")

    chrom_values = subset[chrom_col].to_numpy()
    dad_values = subset[dad_col].to_numpy()

    min_axis = min(np.min(chrom_values), np.min(dad_values))
    max_axis = max(np.max(chrom_values), np.max(dad_values))
    padding = 0.05 * (max_axis - min_axis) if max_axis != min_axis else 0.1

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(chrom_values, dad_values, s=40, alpha=0.75, label="Samples")
    ax.plot(
        [min_axis - padding, max_axis + padding],
        [min_axis - padding, max_axis + padding],
        linestyle="--",
        color="black",
        label="Identity",
    )
    ax.set_xlabel("Chromatogram concentration (mg gDW⁻¹)")
    ax.set_ylabel("DAD concentration (mg gDW⁻¹)")
    ax.set_title(f"{state.capitalize()} scytonemin: Chrom vs DAD")
    ax.legend()
    ax.grid(True, alpha=0.3)

    output_path = outdir / f"{state}_chrom_vs_dad_scatter.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def plot_sample_traces(df: pd.DataFrame, state: str, outdir: Path) -> Path:
    """Plot chromatogram and DAD concentrations across samples to compare patterns."""
    chrom_col = f"chrom_{state}_mg_per_gDW"
    dad_col = f"dad_{state}_mg_per_gDW"

    subset = df[["sample_id", chrom_col, dad_col]].dropna()
    if subset.empty:
        raise ValueError(f"No data to plot for state '{state}'.")

    subset = subset.sort_values("sample_id")
    sample_ids = subset["sample_id"].tolist()
    chrom_values = subset[chrom_col].to_numpy()
    dad_values = subset[dad_col].to_numpy()

    x = np.arange(len(sample_ids))

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(x, chrom_values, marker="o", label="Chromatogram")
    ax.plot(x, dad_values, marker="s", label="DAD")
    ax.set_xticks(x)
    ax.set_xticklabels(sample_ids, rotation=45, ha="right")
    ax.set_ylabel("Concentration (mg gDW⁻¹)")
    ax.set_title(f"{state.capitalize()} scytonemin pattern across samples")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)

    output_path = outdir / f"{state}_chrom_vs_dad_by_sample.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def format_stats_line(stats_obj: ComparisonStats) -> str:
    """Return a concise, readable summary line for CLI output."""
    return (
        f"{stats_obj.state.capitalize():8s} | n={stats_obj.sample_count:3d} | "
        f"mean_diff={stats_obj.mean_difference:6.3f} | "
        f"mean_ratio={stats_obj.mean_ratio:6.3f} | "
        f"paired_t p={stats_obj.t_pvalue:7.3g} | "
        f"slope={stats_obj.slope:6.3f} | intercept={stats_obj.intercept:6.3f} | "
        f"reg p={stats_obj.regression_pvalue:7.3g}"
    )


def run_analysis(csv_path: Path, outdir: Path) -> Dict[str, Dict[str, Path]]:
    """Execute the comparison workflow."""
    df = load_data(csv_path)

    outdir.mkdir(parents=True, exist_ok=True)
    plot_paths: Dict[str, Dict[str, Path]] = {}
    stats_rows: List[ComparisonStats] = []

    print(f"Loaded {len(df)} rows from {csv_path}")
    print("State-by-state comparison metrics:")

    for state in STATES:
        state_stats = compute_state_statistics(df, state)
        stats_rows.append(state_stats)
        print("  " + format_stats_line(state_stats))

        scatter_path = plot_scatter(df, state, outdir)
        trace_path = plot_sample_traces(df, state, outdir)
        plot_paths[state] = {
            "scatter": scatter_path,
            "sample_traces": trace_path,
        }

    summary_df = pd.DataFrame([asdict(row) for row in stats_rows])
    summary_path = outdir / "chrom_vs_dad_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSummary table written to {summary_path}")

    return plot_paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare chromatogram and DAD dry-weight corrected scytonemin concentrations "
            "for total, oxidized, and reduced fractions."
        )
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("Combined_Scytonemin_Concentrations.csv"),
        help="Path to the combined concentrations CSV file.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("chrom_vs_dad_outputs"),
        help="Directory to store generated plots and summary files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plot_paths = run_analysis(args.csv, args.outdir)

    print("\nGenerated plots:")
    for state, paths in plot_paths.items():
        for label, path in paths.items():
            print(f"  {state.capitalize():8s} {label:14s} -> {path}")


if __name__ == "__main__":
    main()
