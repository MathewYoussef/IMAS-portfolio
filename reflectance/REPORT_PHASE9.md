# Phase 9 Reflectance Baseline Report

## 1. Motivation & Design
- **Why Phase 9?** Earlier baselines (Phase 6) relied on simple band averages and showed weak predictive power (LOOCV RMSE > 1 across most combinations). We hypothesised that richer spectral descriptors—peak intensity, band area, continuum depth, FWHM, curvature—might capture concentration trends more faithfully.
- **Data foundations** (`analysis/phase9/inputs/phase9_input_manifest.json`): canonical trimmed spectra (per angle and Σ/Δ composites), wavelength grid, crosswalks, and dose-level concentration tables, plus Phase 6 outputs for comparison.
- **Guiding constraints:** maintain provenance through manifests, keep QA artefacts for auditors, and ensure downstream Phase 7/8 reporting consumes only vetted results.

## 2. Feature Engineering (“How”)
- **Script** `analysis/phase9/scripts/extract_advanced_features.py`
  - **How:** slice wavelength windows (320–480 nm, 360–410 nm, full-grid auto peak) and compute peak, AUC, continuum depth (linear hull), FWHM, and second-derivative curvature for each dose × kind combination.
  - **Why:** capture shape/peak information absent from mean-only descriptors; operate directly on the trimmed spectra (no Savitzky–Golay smoothing), and when a window supplies `fixed_peak_nm` (e.g. 360–410 nm locked at **406 nm**) also report `fixed_peak_reflectance`/`fixed_continuum_depth` so statistics stay aligned across doses.
  - **Outputs:** `analysis/phase9/tables/advanced_band_features.csv` (72 rows) + QA CSVs (`example_band_profiles.csv`, `example_continuum_removed.csv`). Manifest notes document negative Δ behaviour, “no extra smoothing,” and summarises the fixed-peak metrics alongside the traditional descriptors.

## 3. Baseline Modelling (“How & Why”)
- **Script** `analysis/phase9/scripts/run_advanced_baseline.py`
  - **How:** fit linear, ridge, and optional LASSO regressions per modality/kind/band using the enriched feature matrix (fixed-peak columns filled where defined), recording leave-one-dose-out (LOOCV) diagnostics plus coefficient traces.
  - **Why:** ridge remains the reporting model—regularisation mitigates n = 6 overfitting while preserving interpretability—whereas linear and LASSO runs are retained for QA if they are re-enabled in the manifest.
  - **Outputs:** `analysis/phase9/tables/advanced_baseline_summary.csv`, `…/_predictions.csv`, `…/_manifest.json`, and plots under `analysis/phase9/plots/`. The manifest now lists totals, oxidized, and reduced components for every modality/kind/band combination.

### Key Ridge Metrics (LOOCV RMSE)
| Modality | Component | Kind      | Band (Fixed Peak)        | RMSE | Interpretation |
|----------|-----------|-----------|--------------------------|------|----------------|
| chrom    | oxidized  | Sigma     | peak_band_auto (406 nm)  | **0.07** | Best-in-class fit; Σ oxidized ridge lines up with the 406/384 nm rebound pattern across doses. |
| dad      | oxidized  | Sigma     | peak_band_auto (406 nm)  | 0.17 | DAD’s oxidized component mirrors the same fixed peak, albeit with higher residuals. |
| chrom    | total     | 12Oclock  | 360–410 nm (fixed 406 nm)| 0.35 | Angle-level totals carry middling residuals once the peak is locked. |
| chrom    | oxidized  | Δ         | 360–410 nm (fixed 406 nm)| 0.30 | Δ oxidized improves but remains noisier than Σ/angle counterparts. |
| dad      | total     | 6Oclock   | 320–480 nm               | **1.90** | Worst case; 6 o’clock DAD totals still fail to respect the dose ordering. |

*(Full statistics: `analysis/phase9/tables/advanced_baseline_summary.csv` — ridge rows.)*

#### Feature diagnostics
- **Descriptors in play:** each ridge run ingests peak height, peak wavelength, trapezoidal AUC, continuum depth, FWHM, and second-derivative curvature; fixed-peak reflectance/continuum depth are active on windows that lock the wavelength (e.g., 360–410 nm at 406 nm).
- **Σ oxidized (peak band @ 406 nm):** the dominant weights fall on AUC (≈ 1.7 × 10⁻³) and FWHM (≈ −3.3 × 10⁻²), with curvature contributing a smaller stabilising term (≈ −3.9 × 10⁻⁴). This combination delivers LOOCV RMSE ≈ 0.07 with in-sample R² ≈ 0.96 (see the `chrom_oxidized` / `Sigma` / `peak_band_auto` row in `analysis/phase9/tables/advanced_baseline_summary.csv`).
- **Angle totals (12 o’clock 360–410 nm):** fixed-peak reflectance and continuum depth terms appear explicitly. For the oxidized component, the ridge fit leans on peak wavelength shift (≈ 6.4 × 10⁻²) and FWHM (≈ −3.3 × 10⁻²), while the fixed-depth term (≈ −1.4 × 10⁻⁴) helps enforce alignment across doses (`analysis/phase9/tables/advanced_baseline_summary.csv:10`).
- **Reduced forms:** broad-window AUC and FWHM terms dominate, but the curvature coefficients stay near-zero, underscoring the weak link between reduced concentration and reflectance shape (see the Δ reduced rows in `analysis/phase9/tables/advanced_baseline_summary.csv` around the 320–480 nm window, where rmse_loocv sits near 0.44 and R² remains positive).
- **DAD oxidized:** mirrors the Σ oxidized pattern—AUC (≈ 1.2 × 10⁻²) and FWHM (≈ −7.8 × 10⁻²) anchor the fit, curvature (≈ −1.2 × 10⁻³) fine-tunes residuals—yielding LOOCV RMSE ≈ 0.17 with in-sample R² ≈ 0.96 (see the `dad_oxidized` / `Sigma` / `peak_band_auto` row in `analysis/phase9/tables/advanced_baseline_summary.csv`).
- **Second-derivative curvature:** contributes across oxidized combos but carries orders-of-magnitude smaller weights than FWHM/continuum depth, so it acts as a refinement rather than the leading driver.
- **Observed fixed-peak reflectance:** the regenerated tables show 12 o’clock reflectance@406 nm tracing 0.174 → 0.122 → 0.139 → 0.153 → 0.220 → 0.217 (`analysis/phase9/tables/advanced_band_features.csv:26-31`), while the canonical trimmed stats at 384 nm follow 0.164 → 0.111 → 0.126 → 0.140 → 0.210 → 0.201 (`../canonical_dataset/dose_reflectance_stats.csv`). Both sequences bottom out at dose 2, peak at dose 5, and ease slightly at dose 6.

### Dependence Tests (Phase 3/4 recaps)
- Best dCor permutation p-value: **0.138** (`analysis/phase8/tables/comparison_overview.csv`). Pearson/Spearman/Kendall permutation p-values remain ≥ 0.90.
- Interpretation: richer features improve regression accuracy for select combinations but do **not** yield statistically significant global dependence; reflectance alone remains insufficient to recover concentration ordering at n = 6.

## 4. Downstream Integration (“How results propagate”)
- **Phase 7** `analysis/phase7/build_summary.py`: merges Phase 6 baselines and Phase 9 ridge records, tagging rows with `feature_set`, `model_type`, and the new component field so oxidized/reduced metrics survive aggregation.
- **Phase 8** `analysis/phase8/build_report_tables.py`: filters to `feature_set == "phase9"` so `comparison_overview.csv`, plots, and summary markdown reflect ridge runs. The refreshed tables now surface totals vs. oxidized vs. reduced results alongside band metadata.
- **Validation commands:**  
  ```
  python analysis/phase9/scripts/extract_advanced_features.py
  python analysis/phase9/scripts/run_advanced_baseline.py
  python analysis/phase7/build_summary.py
  python analysis/phase8/build_report_tables.py
  python analysis/phase8/generate_report_assets.py
  ```

## 5. Interpretation & Learnings (“What the stats mean”)
1. **Peak-driven improvements** – Σ composites benefit most; ridge lowers LOOCV RMSE to ~0.07–0.17 for the oxidized components (vs >1 previously). Angle totals sit in the 0.30–0.35 range once 406 nm is locked, while Δ retains larger residuals, signalling intrinsic variance/noise in the difference spectrum.
2. **Fixed-wavelength checks (406 nm / 384 nm)** – Locking the band sample keeps comparisons wavelength-aligned; at 12 o’clock the 406 nm trace dips at dose 2, climbs through dose 5, then softens slightly at dose 6 (the 384 nm trace mirrors this), matching the ridge-driven improvements.
3. **Regularisation matters** – Linear fits remain illusory (perfect in-sample, poor LOOCV); LASSO is now wired for QA but, like linear, offers little benefit on six samples. Ridge strikes a practical balance and is therefore the only model surfaced in the report.
4. **Statistical caution** – Despite better fits, permutation tests remain non-significant. The n=6 regime still limits inference; improvements are exploratory rather than definitive.
5. **Other bands worth watching:** broad 320–480 nm windows still provide the most stable handle on reduced concentrations, while the 360–410 nm fixed window remains the clearest hook for angle-level totals. No alternative automated peak surpassed the 406 nm anchor in this run; however, raw spectra confirm the 384 nm reflectance follows the same U-shaped pattern (see `../canonical_dataset/dose_reflectance_stats.csv`).

## 6. Conclusions & Next Steps
- **Conclusion:** Advanced spectral descriptors plus ridge regularisation materially improve predictions for Σ composites and selected angle bands but do not overcome the fundamental sample-size and Δ-spectrum limitations. Reflectance alone still cannot robustly infer concentrations across all modalities.
- **Recommended next steps:**
  1. Highlight Σ vs. Δ behaviour in the Phase 8 narrative (why peaks help Σ but not Δ).
  2. Retain QA artefacts (linear/LASSO metrics, plots) for auditors or future models.
  3. Prioritise new dose/replicate measurements before investing in more complex models; additional data is the clearest path to tighter LOOCV and meaningful statistical tests.
  4. Optional: explore non-linear kernels or multi-block methods once n increases, leveraging the Phase 9 feature tables already in place.

The Phase 9 pipeline, manifests, and reporting assets are now aligned and ready for final dissemination or publication. Ridge-based summaries represent the current best-performing baseline given the available data.
