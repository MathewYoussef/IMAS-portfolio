# Dose-Level Robust Mean Analysis Plan

## Overview
We will explore dose-dependent behaviour after collapsing biological replicates with robust mean estimators. The goal is to contrast dose-level trends against the previously reported replicate-level null findings, and to prepare inputs for downstream reflectance modelling.

## Step 1 — Build Dose-Level Summaries
- Enumerate the six unique (UVA, UVB) dose pairs and highlight that UVB dips only on the final step (0.768 → 0.707 mW·cm⁻²).
- Compute 20 % trimmed means (middle three replicates) and trimmed SDs for chromatogram and DAD totals/oxidized/reduced.
- Generate 95 % percentile bootstrap CIs with 2 000 resamples; document the three-replicate sensitivity in the analysis notes.
- Maintain a single canonical trimmed-mean file (`dose_level_summary.csv`) and reference it elsewhere instead of duplicating CSVs.

## Step 2 — Quantify Trends Using Robust Means
- Treat the dose summary as *n* = 6, compute Pearson/Spearman/Kendall correlations against UVA/UVB, and run weighted least squares fits (weights = 1/trimmed-SD²).
- Report slopes, SEs, R², and 95 % CIs while explicitly characterising the trend as “mostly increasing with a mid/late dip.”

## Step 3 — Maintain Replicate-Level Context
- Keep the replicate-level ANOVA, permutation, and robust tests (n = 30) to show how biological variance obscures dose ordering.
- Summarise the contrast: replicate analyses stay null; trimmed means climb through dose₄/₅ before declining.
- Build visuals overlaying replicates with trimmed means ± CI to show the variance collapse and the late drop.

## Step 4 — Align Chromatogram and DAD at Dose Level
- Using the dose-level robust means, compute Pearson/Spearman correlations and run Deming regression between chromatogram and DAD totals/oxidized/reduced.
- Document the assumed error variance ratio for Deming (default 1:1 unless assay-specific variances are estimated) and provide slopes/intercepts with 95 % CIs to demonstrate assay concordance.

## Step 5 — Prepare Outputs for Reflectance Integration
- Point reflectance workflows to `dose_level_summary.csv` instead of shipping duplicate summaries. Reference the existing metadata table for the (UVA, UVB) mapping and final-step UVB drop.

## Step 6 — Narrative Updates
- Synchronise summaries, captions, and scripts so they:
  - Describe the trend as “mostly increasing with a mid/late dip” and correct the Pearson r lower bound (~0.69 for DAD reduced).
  - Explain the non-additive total vs. oxidized/reduced definitions and call out the negative oxidized baseline.
- Emphasise the dual perspective (replicate heterogeneity vs. dose ordering) and the exploratory nature of conclusions based on six dose means.
