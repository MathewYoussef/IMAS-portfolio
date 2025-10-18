# Reflectance Part 2: Robust Mean Reflectance Data

# Reflectance Part 2: Robust Mean Reflectance Data

## Fundamental framing
- **Objective**: Prove that reflectance-derived rankings reproduce the concentration ordering captured by the dose-level DAD (primary) and chromatogram (secondary) summaries; dose curves remain supporting context, not targets.
- **Question**: Can carefully preprocessed reflectance spectra quantify scytonemin concentrations when replicate-level responses are noisy but dose-level robust means show a consistent UVA/UVB trend?
- **Answer**: Yes—if we curate dose-level (UVA/UVB pair) concentration targets, quantify their dose trends, and validate reflectance→concentration models against those dose summaries while documenting replicate variance.
- **Principle**: Dose-level robust means provide the “pattern truth”; replicate distributions remain available to illustrate biological scatter.

## Phase 0 – Premise check (reaffirm decision to focus on DAD)
- Accept that stress biology yields U-shaped dose responses; this does not invalidate reflectance calibration.
- Scytonemin’s optical fingerprint is well constrained (UVA absorber near 370–384 nm), so reflectance carries a physical signal even when dose misbehaves.
- Decision: pursue reflectance→concentration modeling; treat dose purely as contextual information.

## Phase 1 – Upgrade the concentration targets (no new wet work)
- Reprocess DAD cubes (sample×time×λ) with PARAFAC2 to isolate scytonemin factors and accommodate retention shifts.
- Compare DAD totals vs oxidized/reduced components using method-comparison statistics (Deming/Passing–Bablok, Bland–Altman) to quantify assay offsets.
- Produce precision-weighted latent concentrations from DAD-derived factors (with associated standard errors) for use as the training target. Track dose-level (UVA/UVB) robust means, trimmed SD/MAD, and bootstrap CIs alongside replicate distributions to support dose-level modeling while preserving replicate variance narratives.
- Outcome: a defensible, uncertainty-aware concentration label derived exclusively from DAD spectra, summarised in `canonical_dataset/dose_summary.csv`.

## Phase 2 – Robust reflectance preprocessing (dose aggregates)
- Replicate-level: keep the 150→1 sample/angle trimmed mean (10 % trim), plus within-sample dispersion (trimmed SD/MAD) to quantify measurement noise.
- Dose-level: for each UVA/UVB dose, combine the five sample spectra per angle via a robust estimator (e.g., 20 % trimmed mean, Huber mean). Retain dose-level dispersion (trimmed SD/MAD + bootstrap CI) so reflectance variance aligns with concentration CIs.
- Angles remain separate; build composite curves Σ = (12 o’clock + 6 o’clock)/2 and Δ = 12 o’clock − 6 o’clock from the dose-level summaries to capture overall pigment load vs directional contrast. Store results in `canonical_dataset/dose_reflectance_composites.csv`.
- Outcome: dose-scale reflectance signatures (means + spreads) ready for direct comparison with dose-level concentration targets.

## Identifier policy
- Every aggregated table now carries physical UVA/UVB intensities (`uva_mw_cm2`, `uvb_mw_cm2`) alongside the categorical `dose_id` so dose alignment always occurs on explicit irradiance values.
- Regenerate existing CSV/JSON outputs after pulling these code changes; older artefacts lacking the UVA/UVB columns should be considered deprecated.

## Phase 3 – Dose-level pattern alignment (no DAD spectral features yet)
### Objectives
- Establish whether dose-level reflectance aggregates (per angle, Σ, Δ) reproduce the ranking and dispersion of canonical concentration targets (chromatogram, DAD, latent).
- Produce reproducible feature tables, joins, plots, and statistics that flow cleanly into Phase 4 multivariate work.

### Plan
1. **Data loading & validation**  
   - Build a single loader (script or notebook) that ingests the canonical tables (`canonical_dataset/dose_reflectance_stats.csv`, `…/dose_reflectance_composites.csv`, `…/dose_reflectance_features_angles.csv`, `…/dose_reflectance_features_composites.csv`, `…/dose_dad_concentrations.csv`, `…/dose_summary.csv`, `…/precision_weighted_concentrations_treatment.csv`, `…/dose_level_canonical_summary.csv`).  
   - Enforce canonical dose ordering (`dose_1` → `dose_6`) and assert that Σ/Δ composites have the expected angle metadata.  
   - Record the wavelength window(s) used for downstream summaries (default 320–480 nm; note index start/step).  
   - Flag the stale `reflectance_sample_mapping.csv` in the README so analysts do not attempt sample-level joins.
2. **Reflectance feature derivation**  
   - Reuse the existing dose-level feature tables; derive any additional wavelength-subset metrics (e.g., alternative band integrations, normalized indices) without regenerating the trimmed base tables.  
   - Store derived outputs in `analysis/phase3/` (or equivalent) so the canonical bundle remains frozen.  
   - Verify Σ/Δ vs angle metadata consistency before publishing any new tables.
3. **Concentration target assembly**  
   - Pull dose-level DAD stats (totals/oxidized/reduced, trimmed SD/MAD), chromatogram stats (matching fields from `dose_summary.csv`), and latent precision-weighted totals (plus SE).  
   - Harmonise column names/units to simplify joins (e.g., align `_trimmed_mean` suffixes).  
   - Note uncertainty formats (MAD vs SD vs SE) for later variance comparisons.
4. **Pairwise diagnostics**  
   - Join reflectance features with each concentration modality on `dose_id` / UVA / UVB.  
   - Generate scatter plots, linear fits, and Bland–Altman plots for Σ, Δ, and per-angle metrics against DAD, chromatogram, and latent targets.  
   - Annotate plots with dose labels, replicate counts, and trim fractions; highlight any outliers or crossings.
5. **Variance & uncertainty review**  
   - Compare reflectance SD/MAD against DAD and chromatogram dispersion; document how latent SEs are handled (convert to SD assuming √n or report separately).  
   - Produce bar/violin charts showing which doses exhibit elevated variance in each modality.
6. **Statistical tests**  
   - Compute Pearson/Spearman/Kendall correlations for each reflectance metric vs each concentration target; supplement with permutation or bootstrap p-values given n = 6.  
   - Fit simple OLS models (optionally Deming if variance estimates warrant it) and run leave-one-dose-out fits to gauge robustness.  
   - Reuse and extend the existing correlation machinery (`canonical_dataset/dose_reflectance_dad_correlations.json`) to keep outputs consistent.
7. **Comparative assessment (DAD vs chromatogram vs latent)**  
   - Summarise which concentration modality aligns most closely with reflectance metrics by comparing correlation strength, regression slope/intercept stability, and residual structure.  
   - Document doses where modalities diverge and how reflectance behaves relative to each.  
   - Capture qualitative notes for follow-on Phase 4 modelling choices.
8. **Reporting & archive**  
   - Save all derived tables, plots, and notebooks/scripts to `analysis/phase3/` with a README covering: data sources, loader command, wavelength windows, metric definitions, correlation/regression settings, and the stale sample-mapping caveat.  
   - List follow-up questions or open issues before advancing to Phase 4.

### Outcome
- A reproducible diagnostic bundle showing whether reflectance aggregates preserve the concentration ranking/variance patterns across chromatogram, DAD, and latent targets, ready to anchor Phase 4 multivariate analyses.

## Phase 4 – Modal pattern congruence (dose scale)
- Run PCA/PLS on dose-level reflectance aggregates (per angle and Σ/Δ) and on dose-level concentration summaries (chromatogram, DAD, latent). Compare ordinations via Procrustes and RV (permutation-tested).
- Include Σ/Δ representations when probing cross-modal structure with latent concentrations.
- Apply dependence measures (distance correlation, HSIC) between reflectance features (per angle and composites) and concentrations to detect nonlinear relationships, acknowledging the small sample size.
- Outcome: confirmation that doses occupy similar geometric arrangements across modalities at the group level.

## Phase 5 – Optional multiblock fusion
- Once dose-level congruence is established, consider O2PLS/MB-PLS on dose-aggregated reflectance blocks (12 o’clock, 6 o’clock, Σ, Δ) and dose-level concentration matrices to extract joint vs angle-specific signals.
- Angle blocks remain explicit, but fusion operates solely on dose summaries to respect the destructive sampling design.

## Phase 6 – Minimal, auditable band model
- Compute continuum-removed band depth at ~384 nm and band area over 360–410 nm using dose-level Σ spectra; regress these against dose-level concentrations (weighted by inverse CI²).
- Provides an interpretable baseline showing whether simple band metrics reproduce dose ordering before scaling up to multivariate models.

## Execution Checklist
1. Aggregate reflectance replicates → sample/angle trimmed means (with within-sample dispersion).  
   Result gives us: denoised, sample-level spectra that suppress outlier shots while preserving measurement noise estimates, leading to: reliable inputs for treatment-level summarisation.

2. Collapse sample-level spectra within each treatment (per angle) using a robust estimator and form Σ (angle average) and Δ (angle difference).  
   Result gives us: treatment-scale reflectance means and spreads for UV-facing, UV-opposing, and composite views, leading to: directly comparable spectral signatures for each treatment.

3. Summarise DAD concentrations by treatment (robust means/medians, trimmed SD/MAD).  
   Result gives us: treatment-level concentration benchmarks with uncertainty, leading to: a concentration reference aligned with the reflectance aggregation (`canonical_dataset/dose_dad_concentrations.csv`).

4. Build or confirm the dose crosswalk (archived at `archive/aggregated_reflectance/dose_alignment_crosswalk.csv`) and keep it under version control for audit.  
   Result gives us: explicit linkage between reflectance folders and canonical dose IDs, leading to: confidence that comparisons rely on consistent group definitions.  
   _Note: destructive sampling means sample-level IDs differ between modalities; a 30/30 sample-level join is impossible, so all downstream work operates purely at the treatment level._

5. Compare treatment-level means—scatter plots and correlations of reflectance Σ/Δ/angle band metrics vs DAD total/oxidised/reduced (e.g. using `canonical_dataset/dose_reflectance_dad_summary.csv`).  
   Result gives us: evidence of whether treatments rank the same across modalities, leading to: acceptance or refinement of the reflectance features.

6. Compare treatment-level spreads—bars/violins of reflectance dispersion vs DAD dispersion.  
   Result gives us: insight into whether treatments with high concentration variability also show high spectral variability, leading to: validation (or not) of reflectance sensitivity to within-treatment heterogeneity.

7. Optional leave-one-treatment-out regressions using aggregated metrics (no DAD spectra yet).  
   Result gives us: quantitative confirmation of directional consistency despite n=6, leading to: a go/no-go decision on advancing to multivariate modeling.

8. Rank and pattern concordance checks (future work).  
   Result gives us: Kendall/Spearman rank alignment, z-scored profiles, and pairwise agreement tallies across treatments, leading to: explicit evidence of whether reflectance “high vs low” patterns mirror DAD concentrations despite absolute offsets.

## Canonical dataset
- The harmonised reference bundle lives in `canonical_dataset/`. It contains the trimmed reflectance tables, DAD/latent concentration summaries, composite spectra, comparison outputs, and an aggregated `dose_level_canonical_summary.csv`, all keyed by `dose_id`, `uva_mw_cm2`, and `uvb_mw_cm2`.
- Provenance: `canonical_dataset/README.json` records the build command (`python build_canonical_dataset.py --output canonical_dataset --archive`) and the archived legacy inputs (`archive/aggregated_reflectance/`). Use that JSON plus the archived CSVs for any future audit or re-generation.
- Guardrail: run `python tests/verify_canonical_dataset.py` after rebuilding to diff the canonical artefacts against the archived source within a 1e-10 tolerance and catch accidental drift.
- Legacy artefacts under `archive/aggregated_reflectance/` remain available for audit, but all downstream analysis should depend on the canonical files to guarantee consistent identifiers and metadata.
## Success criteria & deliverables
- Reflectance calibration hits RMSEP and CCC targets under grouped nested CV, with no significant bias in Bland–Altman plots.
- Procrustes/RV show strong agreement (>0.7/0.6 with p < 0.01) between reflectance and DAD ordinations.
- dCor/HSIC confirm dependence between reflectance and concentration.
- Final outputs: robust-mean reflectance dataset, DAD-derived concentration table with uncertainties, model performance reports, pattern-congruence statistics, and documentation of preprocessing/intercept policies.

This roadmap keeps the focus squarely on DAD spectra as the authoritative reference while leveraging reflectance measurements to deliver a defensible, physics-aware quantification of scytonemin.
