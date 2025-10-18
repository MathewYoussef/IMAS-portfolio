# Phase 3 Workspace

This directory hosts artefacts generated while executing the Phase 3 diagnostics.

## Loader workflow

- Run `python analysis/phase3/load_phase3_data.py` to ingest all canonical tables required for Phase 3.  
- The script enforces canonical dose ordering, checks the expected Σ/Δ composites and angle metadata, and writes a summary manifest to `analysis/phase3/phase3_data_manifest.json`.  
- Use the manifest to drive notebooks/noted scripts so every downstream join shares the same wavelength window (`mean_000` corresponds to 300 nm with 0.5 nm spacing by default).

## Derived feature tables

- Run `python analysis/phase3/derive_reflectance_features.py --window 320:480` (default) to generate per-angle and Σ/Δ summaries for the 320–480 nm slice of the 300–600 nm grid.  
- Outputs land in this directory:
  - `phase3_reflectance_features_angles_320_480nm.csv`
  - `phase3_reflectance_features_composites_320_480nm.csv`
  - `phase3_reflectance_feature_manifest.json` documents window bounds (320–480 nm, 0.5 nm step, 321 points) and feature columns (`mean_reflectance`, `median_reflectance`, `min_reflectance`, `max_reflectance`, `range_reflectance`, `area_reflectance`, `mean_std`, `median_std`, `mean_mad`, `median_mad`).  
- Add additional windows by passing multiple `--window lower:upper` flags; the script appends each run to the manifest and writes uniquely suffixed filenames (e.g., `_350_420nm`). Record the bounds/index ranges here so downstream joins can trace the calculations.

## Concentration joins

- Run `python analysis/phase3/assemble_targets.py` to merge the reflectance feature tables with canonical chromatogram (`dose_summary.csv`), DAD (`dose_dad_concentrations.csv`), and latent (`precision_weighted_concentrations_treatment.csv`) summaries.  
- The script reads the feature manifest so every wavelength window is processed automatically, stacks per-angle and Σ/Δ metrics into a single `kind` column, and writes:
  - `phase3_join_reflectance_chrom_{window}.csv`
  - `phase3_join_reflectance_dad_{window}.csv`
  - `phase3_join_reflectance_latent_{window}.csv`
  - `phase3_concentration_target_manifest.json` capturing source files, window bounds, and output paths.  
- Dose metadata (`dose_id`, `uva_mw_cm2`, `uvb_mw_cm2`), sample counts, and uncertainty fields from each modality are preserved; column renames (e.g., `chrom_total_mg_per_gDW_trimmed_mean` → `chrom_total_mean`) are encoded in `analysis/phase3/assemble_targets.py`.

## Diagnostics (Phase 3.4)

- Run `python analysis/phase3/run_diagnostics.py --permutations 10000 --random-seed 42` to generate scatter plots, Bland–Altman plots, dispersion comparisons, and correlation/regression tables for every join listed in `phase3_concentration_target_manifest.json`.  
- Outputs live under `analysis/phase3/diagnostics/`:
  - Plots: `diagnostics/plots/*.png` (one scatter + Bland–Altman per modality × metric × window). Σ/Δ and per-angle kinds share axes with distinct markers/colours.  
  - Tables (CSV):  
    - `correlations.csv` — includes permutation-based p-values for Pearson/Spearman/Kendall, keyed by `window`, `modality`, `metric`.  
    - `regression_results.csv` — OLS slope/intercept/r/p-value/SE per pairing.  
    - `bland_altman_summary.csv` — bias and ±1.96 σ limits of agreement.  
    - `dispersion_comparisons.csv` — long-form table listing reflectance vs concentration dispersion metrics (SD vs SD/MAD; latent comparisons note `target_metric=None` where MAD/SE is unavailable).  
  - `phase3_diagnostics_manifest.json` records plot filenames and diagnostic settings (permutations, random seed).  
- Small-n caveat: n = 6 doses; interpretation relies on permutation tests (default 10 000 shuffles) and leave-one-dose-out logic embedded in regression tables. Update the seed/count in the manifest if you rerun with different parameters.

## Notes

- `reflectance_sample_mapping.csv` is intentionally left untouched because reflectance and concentration samples cannot be joined one-to-one (destructive sampling). Analysts should work strictly at the dose level until the mapping is refreshed.  
- Keep derived tables, plots, and notebooks for each Phase 3 sub-step inside this folder so the canonical dataset remains immutable between rebuilds.  
- Record any wavelength-window overrides or alternative feature definitions alongside the artefacts so later phases can trace parameter choices.
- Phase 4 scripts consume `analysis/phase3/phase3_concentration_target_manifest.json` (and diagnostics outputs) via `analysis/phase4/inputs/phase4_input_manifest.json`; update Phase 3 manifests first whenever new windows/features are generated so downstream phases stay in sync.
