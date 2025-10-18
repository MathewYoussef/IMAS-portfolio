# Phase 5 Workspace

Phase 5 explores multiblock fusion between reflectance aggregates and concentration modalities.

## Structure

- `inputs/phase5_input_manifest.json` references:
  - Phase 4 cache (`analysis/phase4/inputs/phase4_data_cache.pkl`) and metadata (`phase4_matrix_manifest.json`).
  - Phase 4 congruence and dependence manifests for provenance.
  - Primary reflectance blocks (`12Oclock`, `6Oclock`) and optional exploratory blocks (`Sigma`, `Delta`).
  - Default joint/specific component counts and permutation/seed settings.
- `plots/` and `tables/` will collect multiblock outputs (scores/loadings, variance summaries, diagnostics).

## Workflow (planned)

1. Run `python analysis/phase5/load_phase5_data.py` to assemble multiblock structures:
   - Pull z-scored matrices for 12Oclock and 6Oclock (primary blocks) from the Phase 4 cache.
   - Optionally assemble exploratory blocks (e.g., Σ, Δ) when enabled.
   - Record block-level metadata (features, zero-variance drops, dose ordering).
   - Cache: `inputs/phase5_block_cache.pkl` (numpy arrays per block).  
   - Summary: `inputs/phase5_block_manifest.json` listing blocks, shapes, active/dropped columns, and default component/permutation settings.
2. Implement fusion scripts (e.g., `analysis/phase5/run_multiblock.py`) to run O2PLS/MB-PLS or equivalent across the selected blocks and concentration modalities. Document package choice (e.g., pyOPLS/pypls) and component settings in manifests.
   - Current implementation uses a stacked PLS approximation (`python analysis/phase5/run_multiblock.py --include-exploratory` to add Σ/Δ).  
   - Outputs: `tables/multiblock_summary.csv`, `tables/multiblock_scores.csv`, `tables/multiblock_loadings.csv`, `tables/phase5_multiblock_manifest.json`, and plots (`plots/multiblock_scores_{window}_{modality}.png`).  
   - Component counts are capped by block ranks, and the manifest records any automatic downgrades or warnings.
3. Emit tables (variance explained, joint/specific scores) and plots (block score scatter, loading heatmaps, variance contributions) under `tables/` and `plots/`, and track outputs via a manifest (e.g., `phase5_multiblock_manifest.json`).
4. Keep reproducibility controls (component counts, permutations, seed) in the manifest so reruns don’t require code edits.

## Notes

- Primary analyses keep 12Oclock and 6Oclock as separate reflectance blocks; any combined-angle or Σ/Δ analyses should be flagged as exploratory in manifests/outputs.
- Canonical data remain untouched; Phase 5 scripts rely solely on cached matrices from prior phases.
- Document any zero-variance feature drops, permutation iterations, or cross-validation settings so downstream consumers understand robustness levels.
- Record package versions (and installation quirks) for O2PLS/MB-PLS in the README/manifest if third-party libraries are used.
