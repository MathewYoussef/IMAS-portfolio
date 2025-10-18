# Downstream Reflectance Handoff (v3)

This bundle copies the trimmed-mean dose results, supporting documentation, and the key figures the reflectance spectroscopy team needs to align their spectra-to-dose work with the chromatogram/DAD concentration analysis.

## Folder Map
- `data/` — Canonical dose summaries, sequential deltas, UVA/UVB trend coefficients, cross-assay calibration tables, and the replicate-level concentration matrix.
- `docs/` — Narrative descriptions, interpretation notes, planning context, and the UVA/UVB dose metadata table.
- `figures/` — PNG plots illustrating dose trajectories, peak timing, sequential deltas, UV regime context, quadratic curvature, and cross-assay agreement.

## Integration Checklist
1. Aggregate reflectance spectra by dose (20 % trimmed mean + 2 000-draw bootstrap CI) and overlay against `data/dose_level_summary.csv`.
2. Confirm that reflectance trajectories peak at dose₄ (totals/oxidized) and dose₅ (reduced) and soften after the UVB drop (dose₅→₆); refer to `figures/dose_trajectories.png` and `figures/peak_comparison.png`.
3. Compare reflectance-derived step deltas and UVB-weighted trends with `data/dose_pattern_sequential_deltas.csv` and `data/dose_trend_stats.csv`.
4. Use `data/chrom_dad_alignment.csv` and `figures/cross_assay_concordance.png` to calibrate reflectance outputs onto the chromatogram or DAD scale as needed.
5. Carry forward the baseline-subtraction caveat in `docs/results-narrative.md` when interpreting negative oxidized control values.

