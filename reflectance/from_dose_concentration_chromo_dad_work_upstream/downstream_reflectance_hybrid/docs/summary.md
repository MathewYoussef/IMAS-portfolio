# Dose-Level Robust Mean Findings

## Objective
We revisited the dose experiment with 20 % trimmed means to see whether collapsing replicates exposes a clearer trend than replicate-level models (n = 30) managed to detect.

## Data & Robust Aggregation
- Replicates per dose are in `Combined_Scytonemin_Concentrations.csv`.
- `dose_level_summary.csv` supplies trimmed means, trimmed SDs, and 95 % percentile bootstrap CIs derived from 2 000 resamples.
- “Total” reflects direct chromatogram/DAD quantification; oxidized and reduced pools are independent components that need not add up to the total.
- DAD oxidized estimates turn slightly negative at the control because of baseline subtraction—retain them but explain the artefact when presenting.

## UVA–UVB Structure
- UVA steps from 0 to 3.185 mW·cm⁻²; UVB climbs from 0 to 0.768 mW·cm⁻² and only drops to 0.707 mW·cm⁻² on the final dose (`aggregated_reflectance/dose_metadata.md:2-7`).
- The main softening in the means takes place between dose₄ and dose₅ while UVB is still increasing; the UVB decline coincides with the final step.
- Replicate IDs run in reverse dose order; consult the metadata when merging.

## Dose-Level Behaviour
- Chromatogram totals: 0.274 → 0.670 → 0.747 → 0.986 → 0.807 → 0.831 mg·gDW⁻¹.
- Chromatogram reduced reaches 1.790 mg·gDW⁻¹ at dose₅ before dropping to 1.411 mg·gDW⁻¹ at dose₆.
- DAD totals peak at 1.844 mg·gDW⁻¹ (dose₄), dip to 1.440 mg·gDW⁻¹, and recover to 1.585 mg·gDW⁻¹.
- DAD reduced tops out at 2.093 mg·gDW⁻¹ (dose₅) and falls to 1.913 mg·gDW⁻¹.
- Weighted linear slopes versus UVA remain positive (0.091, 0.390, 0.472 mg·gDW⁻¹·(mW·cm⁻²)⁻¹ for chrom total/reduced and DAD reduced respectively) yet the late decline means the trend is “mostly increasing with a late dip.”
- Pearson r with UVA ranges from about 0.64–0.75 for totals and 0.69–0.75 for reduced pools (`dose_trend_stats.csv:2-13`).

## Sequential Deltas
`dose_pattern_sequential_deltas.csv` captures the + + + − + signature for totals and + − + + − for reduced pools; only the first transition has a 95 % CI that excludes zero.

## Chromatogram ↔ DAD Alignment
Cross-assay agreement remains strong (Pearson r = 0.985/0.957/0.948; Deming slopes 1.88/2.03/1.09 for total/oxidized/reduced), so downstream models can treat either assay as interchangeable at the mean level.

## Interpretation & Next Steps
Robust means highlight a strong rise through dose₄, peaks for reduced pools at dose₅, and a late decline once UVB steps down. Describe the behaviour as “mostly increasing with a late dip,” note the non-additive component definitions, and remember the n = 6 sample size keeps results exploratory.

Follow-up:
1. Document handling of the negative oxidized baseline.
2. Reference `dose_level_summary.csv` as the canonical trimmed-mean table instead of duplicating CSVs.
3. Ensure figures, captions, and manuscript text echo the mid-course dip and corrected correlation ranges (lower bound ≈ 0.69 for DAD reduced).
