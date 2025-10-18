# Phase 7 Comparison

This folder summarises metrics across Phases 3–6 to support a consolidated narrative.

## Contents

- `consolidated_report.ipynb` – notebook scaffold (now with a loader cell) referencing diagnostics, congruence, multiblock, and baseline outputs listed in `inputs.json`.
- `inputs.json` – manifest mapping each staged file to its original Phase 3–6 source for provenance.
- Local snapshots of key tables (`phase3_*`, `phase4_*`, `phase5_multiblock_summary.csv`, `phase6_baseline_loocv.csv`, `phase6_baseline_summary.csv`) copied from earlier phases for quick comparison.
- (Optional) future scripts/notebooks for generating combined tables/plots.
- Generated summary tables:
  - `modality_summary.csv` – Phase 3 correlation & regression metrics (total means per modality).
  - `kind_modality_summary.csv` – Phase 4 dCor/Procrustes/RV, Phase 5 multiblock settings, and Phase 6 baseline LOOCV metrics (kind × modality).

## Next Steps

- Use the notebook loader cell to bring the staged tables into data frames and build unified comparison plots/tables.
- Document key findings (e.g., permutation p-values, variance explained, LOOCV errors) for reporting.
