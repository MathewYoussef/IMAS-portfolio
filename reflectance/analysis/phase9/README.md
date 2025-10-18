# Phase 9 Workspace

Phase 9 extends the reflectance feature pipeline to derive peak, continuum-removed, and curvature descriptors directly from the canonical spectra before routing results into modelling and reporting.

## Structure
- `inputs/phase9_input_manifest.json` records the canonical spectra, crosswalk CSVs, Phase 6 baselines, and planned output locations for advanced feature tables and diagnostics.
- `scripts/` hosts reusable tooling (`extract_advanced_features.py`, `run_advanced_baseline.py`) that work strictly from manifest entries.
- `tables/` holds derived feature matrices, QA samples, and validation metrics; `plots/` remains available for optional visual diagnostics.

## Planned Workflow
1. **Assemble inputs** – Load trimmed-mean spectra (`dose_reflectance_stats.csv`, composites) and align replicate metadata via the sample/concentration crosswalks noted in the manifest.
2. **Compute advanced band metrics** – For each configured window, calculate peak height & wavelength, trapezoidal area, continuum-removed band depth/FWHM, and finite-difference slope/curvature statistics; persist results with clear column naming.
3. **Quality checks** – Export representative band curves (CSV) versus dose, validate FWHM/peak indices against wavelength bounds, and log any smoothing/continuum settings in the manifest notes.
4. **Model integration** – Feed the enriched features into `run_advanced_baseline.py` (linear/ridge/LASSO), mirror the LOOCV diagnostics from Phase 6, and export summaries for Phase 7/8 ingestion.

## Notes
- Stick to additive append-only behaviour: regenerate outputs under new filenames rather than editing canonical tables.
- **Savitzky–Golay smoothing is permanently disabled; Phase 9 operates on the denoised canonical spectra exactly as recorded.**
- **Fixed-band metrics:** bands may specify a `fixed_peak_nm` (e.g. 406 nm for the 360–410 nm band) so `fixed_peak_reflectance`/`fixed_continuum_depth` track the *same* wavelength for every dose and angle.
- Δ composites are difference spectra, so peak height and AUC values can be negative; confirm sign conventions before modelling.
- Continuum removal currently uses a straight line between band endpoints (“linear hull”); note this approximation if future work adopts a richer convex-hull approach.
- Auto-peak windows report their configured bounds as `0.0`; rely on `actual_band_lower_nm`/`actual_band_upper_nm` in downstream scripts. When a window specifies `fixed_peak_nm` (e.g., 360–410 nm locked at 406 nm), the extractor records `fixed_peak_reflectance`/`fixed_continuum_depth` so dose-to-dose comparisons stay aligned.
- Advanced baselines now cover chrom/dad totals as well as oxidized and reduced components; ridge fits feed the consolidated report by default.
- Feature QA: peak height, AUC, continuum depth, FWHM, and max curvature all feed the ridge diagnostics, with oxidized Σ components showing the strongest alignment (LOOCV ≈ 0.07) and reduced forms remaining weakest despite the broader 320–480 nm window.
- Linear and LASSO variants are treated as QA runs. To generate them, add the desired entries to `inputs/phase9_input_manifest.json` under `modeling.models` before rerunning the baseline script.
- The sample size constraint (six doses) still limits inference; document uncertainty alongside any apparent improvements.
- Update upstream manifests (Phase 7/8) only after Phase 9 artefacts pass QA so the consolidated report remains reproducible.
