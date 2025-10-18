# Downstream Reflectance Package

This folder aggregates the dose–response summaries, narratives, and figures needed to compare chromatogram/DAD concentrations against reflectance spectra.

## Data (`data/`)
- `Combined_Scytonemin_Concentrations.csv` — replicate-level concentrations (Chrom & DAD).
- `dose_level_summary.csv` — 20% trimmed means, SDs, and 95% bootstrap CIs (canonical table).
- `aggregated_reflectance_dose_summary.csv` — duplicate trimmed summary provided for convenience.
- `dose_summary.csv` — the same trimmed means exported for the reflectance pipeline.
- `dose_pattern_summary.csv` — UVA + UVA² coefficients and Kendall τ trend checks.
- `dose_trend_stats.csv` — linear UVA/UVB regressions with slopes, CIs, and correlations.
- `dose_pattern_sequential_deltas.csv` — bootstrap dose-to-dose change estimates with 95% CIs.
- `chrom_dad_alignment.csv` — cross-assay Pearson r plus Deming slopes/intercepts.
- `dose_metadata.md` — UVA/UVB intensities for doses 1–6.

## Narratives (`narrative/`)
- `results-narrative.md` — polished storyline with peaks, CIs, Δ vs control, and trend interpretation.
- `dose_relationship_results.md` — methods/results reference tying tables to analyses.
- `summary.md` — executive summary of trimmed-mean findings.
- `plan.md` — original analysis plan (context for integration).
- `narative-fixes.md` — tracked revisions and quality gates.

## Figures (`figures/`)
Key PNGs used across the manuscript and handoff decks:
- Dose trajectories, sequential deltas, peak comparison, UV context, cross-assay concordance, reduced vs oxidized, quadratic diagnostics, UVA trends, replicate vs trimmed, and supporting journey figures with captions.

## Suggested Use
1. Align reflectance dose IDs using `data/dose_metadata.md`.
2. Aggregate reflectance means (prefer trimmed/robust) and overlay against `data/dose_level_summary.csv`.
3. Check whether reflectance peaks (dose₄ totals, dose₅ reduced) and step changes match the three resolved transitions in `data/dose_pattern_sequential_deltas.csv`.
4. If magnitudes differ, reference `data/chrom_dad_alignment.csv` for calibration offsets before fitting concentration↔reflectance relationships.
