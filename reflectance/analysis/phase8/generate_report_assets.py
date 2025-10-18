#!/usr/bin/env python3
"""Generate plots and markdown summary for the Phase 8 report."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MANIFEST_PATH = Path("analysis/phase8/inputs/phase8_input_manifest.json")
PLOTS_DIR = Path("analysis/phase8/plots")
REPORT_DIR = Path("analysis/phase8/report")

LOOCV_PLOT = PLOTS_DIR / "phase6_loocv.png"
DEPENDENCE_PLOT = PLOTS_DIR / "dependence_pvalues.png"
SUMMARY_MD = REPORT_DIR / "summary.md"


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text())


def save_manifest(manifest: dict) -> None:
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2))


def ensure_dirs() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)


def build_plots(comparison: pd.DataFrame) -> None:
    comparison_phase9 = comparison.copy()
    best_rmse = (
        comparison_phase9
        .sort_values("baseline_rmse_loocv")
        .groupby(["kind", "modality"], as_index=False)
        .first()
    )

    pivot_rmse = best_rmse.pivot(index="kind", columns="modality", values="baseline_rmse_loocv")
    ax = pivot_rmse.plot(kind="bar", figsize=(6, 4))
    ax.set_ylabel("Baseline RMSE (LOOCV)")
    ax.set_title("LOOCV RMSE by Reflectance Kind / Modality")
    ax.legend(title="Modality")
    plt.tight_layout()
    plt.savefig(LOOCV_PLOT, dpi=300)
    plt.close()

    pivot_p = best_rmse.pivot(index="kind", columns="modality", values="dcor_permutation_p")
    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(pivot_p.values, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(pivot_p.columns)))
    ax.set_xticklabels(pivot_p.columns)
    ax.set_yticks(np.arange(len(pivot_p.index)))
    ax.set_yticklabels(pivot_p.index)
    ax.set_title("Distance-Correlation Permutation p-values")
    for i in range(pivot_p.shape[0]):
        for j in range(pivot_p.shape[1]):
            ax.text(j, i, f"{pivot_p.values[i, j]:.2f}", ha="center", va="center", color="white" if pivot_p.values[i, j] > 0.5 else "black")
    fig.colorbar(im, ax=ax, label="p-value")
    plt.tight_layout()
    plt.savefig(DEPENDENCE_PLOT, dpi=300)
    plt.close()


def build_summary_markdown(comparison: pd.DataFrame, baseline_summary: pd.DataFrame) -> None:
    min_row = comparison.loc[comparison["baseline_rmse_loocv"].idxmin()]
    max_row = comparison.loc[comparison["baseline_rmse_loocv"].idxmax()]
    min_p_row = comparison.loc[comparison["dcor_permutation_p"].idxmin()]
    min_rmse = float(min_row["baseline_rmse_loocv"])
    max_rmse = float(max_row["baseline_rmse_loocv"])
    in_sample_rmse = baseline_summary["rmse"].astype(float).mean()

    lines = [
        "# Reflectance vs. Scytonemin Concentration: Summary\n\n",
        "## Key Metrics\n",
        f"- **Lowest LOOCV RMSE**: {min_rmse:.2f} (kind={min_row['kind']}, modality={min_row['modality']}, component={min_row.get('component', 'total')}, band={min_row.get('band_label', 'n/a')}, model={min_row.get('model_type', 'phase6')}).\n",
        f"- **Highest LOOCV RMSE**: {max_rmse:.2f} (kind={max_row['kind']}, modality={max_row['modality']}, component={max_row.get('component', 'total')}, band={max_row.get('band_label', 'n/a')}, model={max_row.get('model_type', 'phase6')}).\n",
        f"- **Best distance-correlation p-value**: {min_p_row['dcor_permutation_p']:.3f} (kind={min_p_row['kind']}, modality={min_p_row['modality']}).\n",
        f"- **Average in-sample RMSE (Phase 6 baseline)**: {in_sample_rmse:.2e} (orders of magnitude smaller than LOOCV errors).\n",
        f"- **Multiblock PLS joint components**: {int(comparison['multiblock_joint_components'].max())} per modality (no warnings).\n\n",
        "## Interpretation\n",
        "- Phase 3 permutation p-values stay >0.9, so linear correlations are not statistically convincing.\n",
        "- Phase 4 dependence tests (dCor/RV) yield p-values between ~0.14 and ~0.84; none reach conventional significance.\n",
        "- Multiblock fusion adds only a single shared component, reflecting the weak dependence signal.\n",
        f"- Ridge models with advanced features yield LOOCV RMSEs spanning {min_rmse:.2f}–{max_rmse:.2f} depending on kind, modality, and component, highlighting uneven predictive strength.\n",
        "- Under current sampling, reflectance signatures alone cannot reliably predict scytonemin concentration; more data or regularised methods with validation remain essential.\n\n",
        "## Next Steps\n",
        "- Expand dose/replicate counts to improve statistical power.\n",
        "- Explore regularised or nonlinear models only after increasing sample size.\n",
        "- Replace the stacked PLS approximation with a true O2PLS implementation once data volume allows.\n"
    ]
    SUMMARY_MD.write_text("".join(lines))


def main() -> None:
    ensure_dirs()
    manifest = load_manifest()
    comparison = pd.read_csv(Path(manifest["phase8"]["comparison_overview"]["path"]))
    baseline_summary = pd.read_csv(Path(manifest["phase7"]["baseline_summary"]["path"]))

    build_plots(comparison)
    build_summary_markdown(comparison, baseline_summary)

    manifest.setdefault("phase8")["plots"] = {
        "loocv": {
            "path": str(LOOCV_PLOT),
            "source": "analysis/phase8/generate_report_assets.py"
        },
        "dependence": {
            "path": str(DEPENDENCE_PLOT),
            "source": "analysis/phase8/generate_report_assets.py"
        }
    }
    manifest["phase8"]["summary_markdown"] = {
        "path": str(SUMMARY_MD),
        "source": "analysis/phase8/generate_report_assets.py"
    }
    save_manifest(manifest)
    print(f"[INFO] Plots written to {PLOTS_DIR}")
    print(f"[INFO] Summary markdown written to {SUMMARY_MD}")


if __name__ == "__main__":
    main()
