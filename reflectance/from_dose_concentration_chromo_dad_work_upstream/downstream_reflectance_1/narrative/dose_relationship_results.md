# Dose–Concentration Pattern Summary

## Framing
Inside `mean_concentration_cannon_4/` we analyse 20 % trimmed means (middle three replicates) because replicate-level tests (n = 30) were dominated by biological scatter.

## Data Pipeline & Reminders
- Replicates: `Combined_Scytonemin_Concentrations.csv`.
- Trimmed means: `dose_level_summary.csv` (20 % trimmed mean, trimmed SD, 95 % percentile bootstrap CI with 2 000 resamples).
- Trend diagnostics: `dose_trend_stats.csv`, `dose_pattern_summary.csv`.
- Sequential deltas: `dose_pattern_sequential_deltas.csv`.
- Assay alignment: `chrom_dad_alignment.csv`.
- Component definitions: chromatogram/DAD “total” values come from direct quantification; oxidized/reduced pools are spectral components and are not additive, so reduced estimates can exceed totals.
- Baseline artefact: DAD oxidized controls are slightly negative after baseline subtraction—retain them but explain the adjustment.
- Canonical data: reference `dose_level_summary.csv` elsewhere rather than maintaining duplicate trimmed-mean tables.

## Mean-Level Regressions (UVA axis)
Weighted least squares (weights = 1/trimmed-SD²) yield (`dose_trend_stats.csv`):
- Chrom total slope: 0.091 mg·gDW⁻¹·(mW·cm⁻²)⁻¹ (SE 0.075, R² = 0.27).
- Chrom reduced slope: 0.390 (SE 0.085, R² = 0.84).
- DAD total slope: 0.134 (SE 0.132, R² = 0.20).
- DAD reduced slope: 0.472 (SE 0.122, R² = 0.79).

Pearson correlations with UVA span ~0.64–0.75 for totals and ~0.69–0.75 for reduced pools (lower bound ≈ 0.69 for DAD reduced). Slopes remain positive, but the final dose falls, so characterise the behaviour as “mostly increasing with a late dip.”

## Quadratic Fits (UVA + UVA²)
Coefficients from `dose_pattern_summary.csv`:

| metric | β₁ (± SE) | p(β₁) | β₂ (± SE) | p(β₂) | Kendall τ (replicates) | Kendall p |
|:--|--:|--:|--:|--:|--:|--:|
| Chrom Total | 0.632 ± 0.132 | 0.017 | −0.149 ± 0.035 | 0.024 | 0.280 | 0.040 |
| Chrom Oxidized | 0.171 ± 0.208 | 0.470 | −0.039 ± 0.061 | 0.567 | 0.132 | 0.336 |
| Chrom Reduced | 0.653 ± 0.222 | 0.060 | −0.127 ± 0.100 | 0.294 | 0.260 | 0.057 |
| DAD Total | 0.927 ± 0.383 | 0.094 | −0.231 ± 0.108 | 0.122 | 0.255 | 0.062 |
| DAD Oxidized | 0.354 ± 0.535 | 0.556 | −0.073 ± 0.155 | 0.668 | 0.201 | 0.142 |
| DAD Reduced | 0.770 ± 0.238 | 0.048 | −0.145 ± 0.103 | 0.254 | 0.285 | 0.037 |

Interpretation:
- β₁ captures the slope near UVA = 0. Negative β₂ terms encode the softening between dose₄ and dose₅ while UVB is still rising; the UVB decline occurs between dose₅ and dose₆.
- Replicate-level Kendall τ values around 0.26–0.29 (totals/reduced) flag modest ordering even before aggregation.

## Sequential Deltas
Bootstrap deltas from `dose_pattern_sequential_deltas.csv` (95 % CIs):
- Chrom total: +0.40, +0.08, +0.24, −0.18, +0.02.
- Chrom reduced: +0.60 [0.23, 1.30], −0.16, +0.34, +0.27, −0.38.
- DAD total: +0.89, +0.09, +0.35, −0.40, +0.15.
- DAD reduced: +0.89 [0.39, 1.76], −0.23, +0.32, +0.07, −0.18.

Only the initial jump (dose₁→₂) has a CI excluding zero; later intervals overlap zero. Sign sequences: totals `+ + + − +`, reduced `+ − + + −`.

## Interpretation
- Chromatogram totals peak at 0.986 mg·gDW⁻¹ (dose₄), slide to 0.807 mg·gDW⁻¹ (dose₅), and close at 0.831 mg·gDW⁻¹.
- Chromatogram reduced reaches 1.790 mg·gDW⁻¹ (dose₅) before declining to 1.411 mg·gDW⁻¹.
- DAD totals rise to 1.844 mg·gDW⁻¹ (dose₄), drop to 1.440 mg·gDW⁻¹, and rebound to 1.585 mg·gDW⁻¹.
- DAD reduced peaks at 2.093 mg·gDW⁻¹ (dose₅) and falls to 1.913 mg·gDW⁻¹.
- Summarise the response as “mostly increasing with a late dip,” explicitly noting the non-additive component definitions and the negative oxidized baseline.

## Follow-Up
- Reference the delta signs, regression slopes, and Deming summaries when aligning with reflectance analyses or manuscript figures.
- Document how negative oxidized controls are treated (retain, clip, or re-baseline).
- Consider leave-one-out tests on the trimmed means or monotonic trend assessments (e.g., Page’s test) if additional confirmation is required.
