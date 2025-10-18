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
- Initial run deferred in this read-only sandbox; execute locally via `python analysis/phase9/equations/evaluate_equations.py`.
