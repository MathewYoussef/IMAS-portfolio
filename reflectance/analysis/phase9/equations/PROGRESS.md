# Phase 9 Equation Build Log

## 2025-10-18 – Inventory Setup
- Reviewed `analysis/phase9/tables/advanced_baseline_summary.csv` and filtered to ridge rows.
- Captured headline combinations (best LOOCV per target):
  - `chrom_oxidized · Sigma · peak_band_auto` → LOOCV ≈ 0.0748, R² ≈ 0.959.
  - `dad_oxidized · Sigma · peak_band_auto` → LOOCV ≈ 0.1661, R² ≈ 0.956.
  - `chrom_total · 12Oclock · 320_480nm` → LOOCV ≈ 0.3030, R² ≈ 0.436.
  - `dad_total · Sigma · peak_band_auto` → LOOCV ≈ 0.5537, R² ≈ 0.764.
  - `chrom_reduced · Delta · 320_480nm` → LOOCV ≈ 0.4401, R² ≈ 0.798 (flagged QA).
  - `dad_reduced · 6Oclock · peak_band_auto` → LOOCV ≈ 0.3844, R² ≈ 0.831 (flagged QA).
  - `latent_total · Sigma · peak_band_auto` → LOOCV ≈ 0.2896, R² ≈ 0.751.
- Decision: document oxidized + totals as reporting equations; keep reduced entries flagged QA-only unless downstream modelling requests them.
- Materialised `analysis/phase9/equations/model_inventory.json` with seven curated entries and a `reporting_role` flag.

## 2025-10-18 – Equation Reference
- Authored `analysis/phase9/equations/equation_reference.md` with per-model tables (feature order, coefficients, intercept, metrics, reporting role).
- Manually recreated `analysis/phase9/equations/model_inventory.json` after scripted writes were blocked to preserve reporting vs QA tagging.
- Cross-checked oxidized entries map to rows 80–83 of `analysis/phase9/tables/advanced_baseline_summary.csv`; totals/reduced sections reference their respective bands.

## 2025-10-18 – Sanity-Check Script
- Added `analysis/phase9/equations/evaluate_equations.py` to recompute predictions from curated coefficients.
- Script reads the feature matrix, canonical targets, and baseline predictions to validate in-sample values within the configured tolerance.
- Initial run deferred in this read-only sandbox; run locally via `python analysis/phase9/equations/evaluate_equations.py`.

## 2025-10-18 – Evaluation Run
- Executed `python analysis/phase9/equations/evaluate_equations.py --role reporting --output analysis/phase9/equations/reports/evaluation_summary.csv`.
- Console output and CSV agree with catalogue metrics; maximum delta between recomputed and stored predictions was < 2e-16 for all reporting models.
- Stored artefact in `analysis/phase9/equations/reports/evaluation_summary.csv` for auditors.

## 2025-10-18 – Feature Subset Sweep
- Added `analysis/phase9/equations/feature_subset_sweep.py` to exercise predefined single/two-feature ablations (ridge α=0.1) against the reporting inventory.
- Ran `python analysis/phase9/equations/feature_subset_sweep.py --role reporting --include-full` and captured outputs under `analysis/phase9/equations/outputs/`.
  - Summary: `feature_subset_sweep_summary.csv`
  - Per-fold diagnostics: `feature_subset_sweep_folds.csv`
- Key observation: Σ oxidized ridge is effectively unchanged with `fwhm_nm` alone (ΔRMSE ≈ -3e-05), while chrom_total@12Oclock benefits from width-driven subsets (RMSE ≈ 0.259 with `fwhm_nm` only).

## 2025-10-18 – Evaluation Replay
- From the project root, ran `python analysis/phase9/equations/evaluate_equations.py --inventory analysis/phase9/equations/model_inventory.json --features analysis/phase9/tables/advanced_band_features.csv --targets canonical_dataset/dose_level_canonical_summary.csv --predictions analysis/phase9/tables/advanced_baseline_predictions.csv --role reporting --output analysis/phase9/equations/reports/evaluation_summary.csv`.
- Output confirms the five reporting equations reproduce in-sample metrics exactly (max prediction delta ≈ 1e-16; all within tolerance).
- Persisted the summary to `analysis/phase9/equations/reports/evaluation_summary.csv` for auditors.

## 2025-10-18 – Validation Harness (Task 5)
- Added `analysis/phase9/equations/validation_harness.py` to replay held-out folds using the curated coefficients.
- Ran `python analysis/phase9/equations/validation_harness.py --role reporting --summary-output analysis/phase9/equations/outputs/validation_summary.csv --fold-output analysis/phase9/equations/outputs/validation_folds.csv`.
- Results confirm reporting models reproduce catalogue RMSE/MAE exactly (max residual ≈ 0.34 for the 12 o’clock total); per-dose residuals archived for auditors.

## 2025-10-18 – Documentation Refresh (Task 6)
- Folded subset and validation findings into `analysis/phase9/equations/equation_reference.md`, linking to `feature_subset_sweep_*` and `validation_*` outputs plus the Phase 8/Phase 9 reports that cite headline metrics.
- Listed verification artefacts under a new “Verification Outputs” section for quick auditor access.
