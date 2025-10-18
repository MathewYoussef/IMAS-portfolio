# Exploring the Consequence of Leveraging a Robust Mean for Concentration Data

Exploring the consequence of leveraging a robust mean for concentration data, collapsing upon replicates, and thereby excluding innate biological variance has led to the following exploration. UVA increased stepwise from 0.000 to 3.185 mW·cm⁻², while UVB rose from 0.000 to 0.768 mW·cm⁻² and then eased to 0.707 mW·cm⁻² at dose₆. Across both Chrom and DAD readouts, concentrations climb through the mid-range and soften at the top end; reduced fractions carry most of the signal, with oxidized pools smaller and less stable. Using 20 % trimmed means to aggregate n = 5 replicates per dose and percentile bootstrap 95 % CIs based on 2,000 resamples, we re-examined how concentration responds to the UVA–UVB dose schedule.

## Trimmed Mean Profiles

Across both chromatogram ("Chrom") and diode array detector ("DAD") readouts, trimmed means reveal a common shape: concentrations climb through the mid-range, then soften at the top end. The reduced fraction carries most of the signal; oxidized components are smaller and less stable.

- Chrom Total peaked at 0.986 mg·gDW⁻¹ at dose₄ (95 % CI 0.668–1.197).
- Chrom Reduced peaked at 1.790 mg·gDW⁻¹ at dose₅ (95 % CI 0.823–2.917).
- DAD Total peaked at 1.844 mg·gDW⁻¹ at dose₄ (95 % CI 1.423–2.212).
- DAD Reduced peaked at 2.093 mg·gDW⁻¹ at dose₅ (95 % CI 1.056–3.574). Totals partially rebound at dose₆, whereas reduced pools continue downward as UVB steps down (0.768→0.707 mW·cm⁻²) while UVA continues upward.
- Chrom Oxidized peaked at 0.611 mg·gDW⁻¹ at dose₄ (95 % CI 0.062–1.081).
- DAD Oxidized peaked at 0.951 mg·gDW⁻¹ at dose₄ (95 % CI 0.191–1.656).

## Changes Relative to Control

Relative to the control (dose₁), effect sizes were: Chrom Total 0.274→0.986 (Δ = 0.712 mg·gDW⁻¹; +260 %), Chrom Reduced 0.740→1.790 (Δ = 1.050; +142 %), DAD Total 0.520→1.844 (Δ = 1.324; +254 %), DAD Reduced 1.050→2.093 (Δ = 1.043; +99 %), DAD Oxidized −0.114→0.951 (Δ = 1.065). The negative control mean (−0.114 mg·gDW⁻¹; 95 % CI −0.532 to 0.757) is consistent with a baseline-subtraction artifact—values may cross zero due to baseline correction—so we retain it and explain the sign.

## Peak-Dose Uncertainty (95% Bootstrap CIs)

Peak-dose uncertainty illustrates the late-dose softening:

- Chrom Total at dose₄: 0.986 mg·gDW⁻¹ (95 % CI 0.668–1.197).
- Chrom Reduced at dose₅: 1.790 mg·gDW⁻¹ (95 % CI 0.823–2.917).
- DAD Total at dose₄: 1.844 mg·gDW⁻¹ (95 % CI 1.423–2.212).
- DAD Reduced at dose₅: 2.093 mg·gDW⁻¹ (95 % CI 1.056–3.574).
- Chrom Oxidized at dose₄: 0.611 mg·gDW⁻¹ (95 % CI 0.062–1.081).
- DAD Oxidized at dose₄: 0.951 mg·gDW⁻¹ (95 % CI 0.191–1.656).

## Stepwise Changes

Bootstrap differences between successive doses show the signature step pattern:

- Totals (Chrom and DAD): + + + − +.
- Reduced (Chrom and DAD): + − + + −.

Of the sequential dose-to-dose differences, three transitions are resolved at the 95 % level: Chrom Reduced dose₁→₂ (Δ = 0.601 mg·gDW⁻¹; 95 % CI 0.235–1.296), DAD Reduced dose₁→₂ (Δ = 0.891 mg·gDW⁻¹; 95 % CI 0.391–1.760), and DAD Reduced dose₃→₄ (Δ = 0.319 mg·gDW⁻¹; 95 % CI 0.001–1.133).

All other step CIs straddle zero, which is expected at this sample size and does not contradict the smoother trend analyses below.

## Trend Across Dose (Rank-Based and Parametric)

Using weighted least squares fits with UVA and UVA² (reporting the linear term β₁) and confirming with Kendall’s τ, we observe:

**Chrom Total**: linear slope 0.632 ± 0.132 (p = 0.017); quadratic term −0.149 ± 0.035 (p = 0.024); Kendall’s τ = 0.280 (p = 0.040).  
**DAD Reduced**: linear slope 0.770 ± 0.238 (p = 0.048); quadratic term −0.145 ± 0.103 (p = 0.254); Kendall’s τ = 0.285 (p = 0.037).  
**Chrom Reduced**: linear slope 0.653 ± 0.222 (p = 0.060); quadratic term −0.127 ± 0.100 (p = 0.294); Kendall’s τ = 0.260 (p = 0.057).  
**DAD Total**: linear slope 0.927 ± 0.383 (p = 0.094); quadratic term −0.231 ± 0.108 (p = 0.122); Kendall’s τ = 0.255 (p = 0.062).  
**Chrom Oxidized**: linear slope 0.171 ± 0.208 (p = 0.470); quadratic term −0.039 ± 0.061 (p = 0.567); Kendall’s τ = 0.132 (p = 0.336).  
**DAD Oxidized**: linear slope 0.354 ± 0.535 (p = 0.556); quadratic term −0.073 ± 0.155 (p = 0.668); Kendall’s τ = 0.201 (p = 0.142).

Interpretation: Chrom Total shows a clear positive trend with concave-down curvature (both terms p < 0.05), consistent with the mid-range rise followed by a late dip. DAD Reduced also trends upward (linear p ≈ 0.048; Kendall’s τ p ≈ 0.037). Chrom Reduced and DAD Total are directionally similar but narrowly miss 0.05 thresholds at n = 6 doses; oxidized components show no significant trend.

## Relation to UVA and UVB Axes

Relation to UVA and UVB axes. Regressions of trimmed means on the physical dose axes show a UVB-dominant association for the reduced pools, with UVA still contributory: Chrom Reduced vs UVB slope 1.218 (95 % CI 0.718–1.717; r = 0.900); vs UVA slope 0.390 (0.154–0.626; r = 0.748). DAD Reduced vs UVB slope 1.537 (0.697–2.378; r = 0.824); vs UVA slope 0.472 (0.133–0.812; r = 0.689). Totals are positively related to both axes (r ≈ 0.65–0.83) with wider slope CIs at this sample size.

## Cross-Assay Agreement (Chrom vs DAD)

Cross-assay agreement (Chrom↔DAD) is high: Pearson r = 0.985 (Totals), 0.957 (Oxidized), 0.948 (Reduced). Deming fits (error ratio = 1) yield slopes ≈ 1.88 (Total), 2.03 (Oxidized), and 1.09 (Reduced) with intercept ≈ 0.33 mg·gDW⁻¹ for the Reduced fraction, indicating near-proportionality for Reduced and calibration offsets for the other pools.

## Synthesis

The robust-mean analysis reveals a mostly increasing trajectory through dose₄ followed by a late softening. The reduced component is the main driver of the dose signal and correlates most strongly with UVB. Assays agree on shape across doses; absolute scaling differs in totals and oxidized pools, while the reduced pool is nearly one-to-one. Given n = 5 per dose, stepwise CIs are wide, so inference should lean on the trend tests and axis regressions rather than any single jump.

## Notes on Artifacts and Handling

The slight negative DAD oxidized mean at the control (−0.114 mg·gDW⁻¹; 95 % CI −0.532 to 0.757) arises from baseline subtraction and is retained (not truncated to zero) to preserve unbiased uncertainty.
