# Downstream Reflectance Handoff — Hybrid Bundle

This folder merges every unique artefact from `downstream_reflectance_1`, `downstream_reflectance_2`, and `downstream_reflectance_3` so the reflectance team can work from a single, comprehensive package.

## Structure
- `data/` – Canonical trimmed-mean summaries, sequential deltas, UVA/UVB trend stats, cross-assay calibration tables, the replicate matrix, and the aggregated reflectance table.
- `docs/` – Narrative results, summary notes, plan, Deming alignment guidance, figure captions, the UVA/UVB dose metadata, and the `narative-fixes.md` errata sheet.
- `figures/` – Union of all downstream-ready PNGs (dose trajectories, peak comparison, UV regime context, sequential deltas, Deming concordance, quadratic diagnostics/curvature, reduced vs oxidized, replicate overlays, etc.).

## Suggested Use
1. Use `data/dose_level_summary.csv` and `data/dose_summary.csv` as the canonical trimmed-mean tables for concentration ↔ reflectance comparisons.
2. Align reflectance regressions to UVA/UVB axes with `data/dose_trend_stats.csv` and `data/dose_pattern_summary.csv`.
3. Mirror stepwise behaviour against `data/dose_pattern_sequential_deltas.csv` and the companion graphics in `figures/`.
4. Calibrate reflectance-derived concentrations with `data/chrom_dad_alignment.csv` and the plots `figures/cross_assay_concordance.png` and `figures/figure5_assay_concordance.png`.
5. Keep narrative and methodological context handy via `docs/results-narrative.md`, `docs/summary.md`, `docs/plan.md`, and the UVA/UVB metadata in `docs/dose_metadata.md`.

