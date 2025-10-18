# Phase 4 Workspace

Phase 4 examines dose-level congruence between reflectance aggregates and concentration modalities using the Phase 3 join outputs.

## Structure

- `inputs/phase4_input_manifest.json` references:
  - `analysis/phase3/phase3_concentration_target_manifest.json` (reflectance↔concentration joins).
  - `analysis/phase3/diagnostics/phase3_diagnostics_manifest.json` (Phase 3.4 plots/tables for context).  
  Extend `default_windows` when additional reflectance windows are generated in Phase 3.
- `plots/` and `tables/` will hold PCA/PLS score plots, Procrustes/RV results, and dependence metrics (dCor/HSIC).

## Workflow

1. Run `python analysis/phase4/load_phase4_data.py` to ingest the Phase 3 join tables, extract matrices per reflectance `kind` (Σ, Δ, 12Oclock, 6Oclock), and z-score them while recording scaling parameters (per column) and dose ordering.  
   - Cache: `inputs/phase4_data_cache.pkl` stores z-scored matrices plus scaling metadata.  
   - Summary: `inputs/phase4_matrix_manifest.json` lists windows, feature columns, zero-variance columns, dose order, and default permutation/seed settings.
2. Run `python analysis/phase4/run_congruence.py --permutations 10000 --random-seed 42` (defaults pulled from the manifest) to compute PCA, Procrustes, RV, and two-component PLS statistics for every window/kind/modality combination.  
   - Outputs:  
     - Plots (`plots/pca_scores_{window}_{kind}_{modality}.png`) showing reflectance vs concentration PCA scores annotated by dose.  
     - Tables: `tables/congruence_summary.csv`, `tables/pca_scores.csv`, `tables/pca_loadings.csv`, plus `tables/phase4_congruence_manifest.json` capturing settings (components, permutations, seed) and generated artefacts.  
   - The congruence summary records variance explained, permutation exceedances, and any warnings (e.g., rank issues or shape adjustments).
3. Run `python analysis/phase4/run_dependence.py --permutations 10000 --random-seed 42` to compute distance correlation and HSIC for each window/kind/modality pair.  
   - Outputs: `tables/dependence_metrics.csv` (tidy table with statistics, permutation exceedances, kernel bandwidths, zero-variance drops) and `tables/phase4_dependence_manifest.json`.  
   - HSIC uses an RBF kernel with the median heuristic in z-scored space; the chosen bandwidths are logged per comparison.
4. Write outputs to `plots/` and `tables/`, and maintain a manifest (e.g., `phase4_congruence_manifest.json`) documenting inputs, parameters (components, permutations, seeds), and generated artefacts.

## Notes

- Canonical data remain untouched; Phase 4 scripts consume only Phase 3 artefacts.  
- Document all seeds, permutation counts, and bandwidth choices in manifests/README for reproducibility.  
- Small sample size (n = 6 doses) warrants cautious interpretation; emphasise this in reports and manifests.  
- When standardising matrices, note that scaling is performed per window and `kind`, and store the parameters so later phases can reuse or invert the transformation if needed.
