# Phase 9 Equation Documentation & Validation Plan

## Objectives
1. Document the ridge-based concentration equations (totals, oxidized, reduced) derived from Phase 9 features.
2. Provide tooling to re-evaluate those equations, explore smaller feature subsets, and verify predictions against known concentrations.

## Tasks
1. **Scope & Inventory**
   - Confirm which Phase 9 ridge models to enshrine (default: the rows already cited in Phase 8/Phase 9 reports—Σ oxidized peak_band_auto plus any totals/reduced needed downstream).
   - Decide upfront whether reduced components (LOOCV > 0.6) stay QA-only or get documented alongside totals/oxidized.
   - Decide if QA variants (linear/LASSO, alternate bands) should be included; if so, tag each entry as `reporting` or `qa` up front.
   - Finalize folder layout for equation artefacts (`analysis/phase9/equations/`) and mirror the existing manifest structure.
2. **Equation Reference Doc**
   - Auto-generate a Markdown summary capturing, per curated model, with table snippets and clickable references (pointing to source rows such as `analysis/phase9/tables/advanced_baseline_summary.csv:80`):
     - Target/modality/kind/band/component
     - Feature columns, coefficients, intercept
     - Fit metrics (RMSE, MAE, R², LOOCV RMSE/MAE)
     - Preprocessing assumptions
     - Brief prose on fixed-peak descriptors (e.g., 406 nm handling) and any fixed_peak_* terms
3. **Sanity-Check Script**
   - Implement script that loads `advanced_band_features.csv` + canonical concentrations, recomputes predictions using the exact feature ordering stored in `feature_columns`, and flags any drift beyond tight tolerances.
   - Provide an optional flag to rebuild features directly from `canonical_dataset/dose_reflectance_stats.csv` for full replay; verify manifest hashes before executing.
   - Persist coefficients/intercepts in a JSON manifest; expose CLI options to focus on specific targets/components and reuse the existing LOOCV splits.
4. **Feature-Subset Exploration (Optional Extension)**
   - Extend/check script to run LOOCV sweeps over a predefined shortlist of ablations for the headline Σ oxidized equation (e.g., AUC-only, FWHM-only, AUC+FWHM) with a reproducible seed; capture per-fold errors.
   - Persist results under `analysis/phase9/equations/outputs/`; note in the doc that n = 6 makes exhaustive subsets noisy and require a meaningful LOOCV margin before recommending slimmed equations.
5. **Validation Harness (Optional Extension)**
   - Replay held-out folds explicitly (with optional bootstrap noise around existing spectra) rather than drawing random unlabelled samples.
   - Emit per-dose errors plus aggregate RMSE/MAE with a reproducible seed; avoid validating on synthetic spectra far from the manifold.
6. **Documentation Updates**
   - Append findings (best equations, tolerances, minimal-feature results) to the Markdown doc; link scripts/outputs and cross-link to `REPORT_PHASE9.md` and the Phase 8 summary.
   - Add a “Verification Outputs” section summarising generated CSV/plots; note which models serve as QA vs. production alongside the small-sample caveat.
7. **Integration Hooks**
   - Update Phase 9 README/manifest (and, if needed, `analysis/phase9/manifest.json`) accordingly.
   - Add a feature flag so Phase 7/8 builders can optionally ingest equation outputs without disrupting existing artefacts.

## Open Questions
- Which tolerance thresholds define an “acceptable” prediction (absolute/relative)?
- Should linear/LASSO equations be included alongside ridge for comparison?
- Do we promote minimal-feature fits into Phase 7/8, or treat them as exploratory QA only?
