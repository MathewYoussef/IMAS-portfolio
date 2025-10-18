# Reflectance Workspace

This directory mirrors the full analysis stack used to build the Phase 9
reflectance equations and the canonical concentration cross-checks. It is
intended to travel inside the IMAS portfolio repo as a self-contained workspace;
no files from the top-level IMAS tooling are required to reproduce the results
here.

## Layout

| Path | Description |
| --- | --- |
| `analysis/phase3` – `analysis/phase9` | Historical analysis phases (features, dependence tests, modelling, reports). Each phase ships with its own `README.md` and manifests describing the inputs it expects. |
| `analysis/phase9/equations/` | Curated ridge equations, documentation, and verification tooling. See `equation_reference.md` for every coefficient table and audit trail. |
| `canonical_dataset/` | Canonical dose-level tables used across all scripts. |
| `archive/` | Legacy copies of earlier equation drops retained for provenance. |
| `tests/verify_canonical_dataset.py` | Regression test that compares regenerated artefacts against the immutable archive bundle. |
| Script entry points (e.g. `aggregate_reflectance.py`, `build_canonical_dataset.py`) | CLI utilities referenced throughout the manifests. |

## Quickstart

```bash
cd reflectance
python3 -m venv .venv && source .venv/bin/activate
pip install pandas numpy

# Rebuild advanced features and ridge baselines
python analysis/phase9/scripts/extract_advanced_features.py
python analysis/phase9/scripts/run_advanced_baseline.py

# Update downstream summaries
python analysis/phase7/build_summary.py
python analysis/phase8/build_report_tables.py
python tests/verify_canonical_dataset.py
```

The commands above regenerate the Phase 9 advanced feature matrix, rerun the
ridge baselines, and rebuild the Phase 7/8 artefacts that surface the headline
metrics (Σ oxidized LOOCV ≈ 0.07; DAD oxidized LOOCV ≈ 0.17).

## Equation Verification

The `analysis/phase9/equations/` folder includes:

- `model_inventory.json` – curated ridge coefficients with reporting vs QA tags.
- `equation_reference.md` – human-readable tables with metrics, feature subsets,
  validation pointers, and links back to Phase 8/9 reports.
- `evaluate_equations.py` – recomputes in-sample predictions and checks them
  against `advanced_baseline_predictions.csv` (see
  `reports/evaluation_summary.csv`).
- `feature_subset_sweep.py` – LOOCV ablation runner for single/two-feature
  subsets (`outputs/feature_subset_sweep_summary.csv`).
- `validation_harness.py` – regenerates per-dose residuals for each curated
  equation (`outputs/validation_summary.csv` & `validation_folds.csv`).
- `PROGRESS.md` – step-by-step audit trail covering every artefact.

## Reports & Documentation

- `REPORT_PHASE9.md` captures the Phase 9 baseline study (feature engineering,
  modelling, interpretation).
- Earlier phases keep their own READMEs under `analysis/phaseX/README.md`.
- `Introdcutin_to_reflectance_vs_concentration.md` and
  `reflectance_part_2_robust_mean_reflectance_data.md` summarise the broader
  reflectance vs concentration narrative.

## Validation Summary

- Σ oxidized ridge (chrom, Σ, 406 nm) remains the best performer with LOOCV
  RMSE ≈ 0.0748 and max held-out residual ≈ 0.074.
- DAD oxidized Σ follows at LOOCV ≈ 0.166; totals stay noisier (chrom 12 o’clock
  LOOCV ≈ 0.303, DAD Σ ≈ 0.554).
- Feature subset sweeps show Σ oxidized is stable with `fwhm_nm` alone, while
  totals require width/peak combinations to stay within ≈ 0.04 RMSE of the full
  model.
- Dependence tests (Phase 4) remain non-significant; ridge equations are useful
  for QA but do not overturn the “n = 6” statistical limitation.

All provenance manifests live alongside their generating scripts so auditors can
rebuild any table or plot deterministically.
