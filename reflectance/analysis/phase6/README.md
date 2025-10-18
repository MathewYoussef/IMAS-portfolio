# Phase 6 Workspace

Phase 6 builds a minimal, interpretable baseline linking simple reflectance band features to concentration targets.

## Structure

- `inputs/phase6_input_manifest.json` references:
  - Phase 4 cache/metadata (reflectance matrices).
  - Phase 4 dependence/congruence manifests (context).
  - Phase 5 multiblock manifest (for provenance).
  - Band definitions (default 360–410 nm), model defaults, and RNG seed.
- `plots/` and `tables/` store feature tables, baseline model outputs, and diagnostics.

## Workflow (planned)

1. Extract band features (`analysis/phase6/extract_band_features.py`):
   - Load the Phase 4 cache; compute band-level summaries (mean, median, min/max, band area, band depth if defined) for each reflectance kind (12Oclock, 6Oclock, Σ, Δ).
   - Write `tables/band_features_{window}.csv` with dose/kind metadata and band features.
2. Run baseline models (`analysis/phase6/run_baseline_models.py`):
   - Fit simple regressions (e.g., linear) between band features and concentration targets (chrom_total, dad_total, latent_total).
   - Support leave-one-dose-out or simple cross-validation; log RMSE/MAE, slopes/intercepts, and residual diagnostics.
   - Current implementation uses linear regression with leave-one-dose-out predictions. Outputs:
     - `tables/baseline_summary.csv`
     - `tables/baseline_predictions.csv`
     - `tables/phase6_baseline_manifest.json`
     - `plots/baseline_fit_{target_kind}.png`
3. Maintain a manifest (`tables/phase6_baseline_manifest.json`) documenting inputs, band windows, model settings, permutation/validation choices, and outputs.

## Notes

- Canonical data remain untouched; Phase 6 consumes cached reflectance matrices from Phase 4.
- Feature extraction should note the exact wavelength subset and any scaling/normalisation applied.
- Small sample size (n = 6 doses) warrants cautious interpretation; document cross-validation or permutation choices.
- If exploratory windows (beyond 360–410 nm) are added, extend the manifest and note them as such.
