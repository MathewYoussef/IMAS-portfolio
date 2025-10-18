# Narative Fixes

## A) Priority fixes (do these first)

### P0 — Contradictions, omissions, and misstatements

#### Resolve UVA-slope inconsistency (linear UVA vs. UVA + UVA²)

**Problem:** `dose_relationship_results.md:17–21` reports weighted UVA-only slopes (0.091–0.472) that contradict the UVA + UVA² fits (β₁ ≈ 0.63–0.93 for totals; 0.77 for DAD reduced).  
**Action:** Choose one of two paths and apply consistently:  
- **Preferred:** Replace the UVA-only paragraph with the UVA + UVA² coefficients and p-values (the canonical model used elsewhere).  
- **Alternate:** Keep the UVA-only slopes, but label them explicitly as a separate analysis (“axis-specific linear regression, WLS”) and move them to a short Methods-adjacent subsection to avoid co-mingling models in the Results.  
**Acceptance:** Nowhere in the Results are linear UVA numbers presented as if they were the UVA + UVA² β₁; model labels appear in the same sentence as the coefficients.

#### Restore the dose progression and the “reduced-fraction dominates” narrative

**Problem:** Results omit the explicit UVA/UVB sequence and the canonical statement about reduced pools carrying the signal.  
**Action:** Add a one-sentence dose geometry and a one-sentence takeaway immediately after the opening line about robust means.  
**Acceptance:** Both sentences appear before any peak/delta claims.

#### Insert relative-to-control effect sizes (Δ and % increase)

**Problem:** Reader cannot gauge magnitude (e.g., +260% Chrom Total; +99% DAD Reduced).  
**Action:** Add a compact paragraph listing control→peak values, Δ, and % for each metric.  
**Acceptance:** All five metrics include start value, peak value, absolute Δ, percent rise, with units.

#### Document peak-dose uncertainty (95% bootstrap CIs)

**Problem:** Peaks listed without their canonical CIs.  
**Action:** Append the 95% CI after each peak mean (Chrom Total dose₄ 0.668–1.197, etc.).  
**Acceptance:** Every peak mean in Results is immediately followed by its 95% CI and dose label.

#### Complete the stepwise story (include the third resolved jump)

**Problem:** DAD Reduced dose₃→₄ significant increase (CI 0.001–1.133) is missing.  
**Action:** Add this as the third resolved transition alongside the two early jumps.  
**Acceptance:** The three resolved steps are named verbatim (Chrom Red 1→2; DAD Red 1→2; DAD Red 3→4).

#### Put UVA & UVB axis regressions and cross-assay agreement back in

**Problem:** UVB-dominant association and Chrom↔DAD concordance are asserted elsewhere but missing here.  
**Action:** Add a short “Relation to UVA and UVB” paragraph quoting slopes, CIs, and r for Reduced; then a “Cross-assay agreement” paragraph with Pearson r and Deming slope/intercept (quote exact coefficients; don’t round into “≥ 0.95”).  
**Acceptance:** Reduced vs. UVB and vs. UVA slopes with CIs appear; Pearson r are 0.985/0.957/0.948, and Deming parameters are present including the ≈ 0.33 mg·gDW⁻¹ intercept for Reduced.

#### Fix the late-dose phrasing error

**Problem:** Text says “UVB steps back up” at dose₆; UVB actually decreases (0.768→0.707), while UVA increases.  
**Action:** Replace that clause with a sentence reflecting “UVA up, UVB down.”  
**Acceptance:** No sentence in Results states or implies that UVB increases at dose₆.

#### Keep the negative oxidized control note

**Problem:** Must retain and explain baseline-subtraction artifact.  
**Action:** Keep the sentence verbatim and—if space permits—add “values may cross zero due to baseline correction.”  
**Acceptance:** The negative control is explicitly acknowledged with its CI.

## B) Concrete text patches (ready-to-paste)

Apply these to `dose_relationship_results.md` (and mirror in `summary.md` where relevant).

### B1) Intro block (dose geometry + reduced-dominance)

Add immediately after the first sentence about robust means:  
`UVA increased stepwise from 0.000 to 3.185 mW·cm⁻², while UVB rose from 0.000 to 0.768 mW·cm⁻² and then eased to 0.707 mW·cm⁻² at dose₆. Across both Chrom and DAD readouts, concentrations climb through the mid-range and soften at the top end; reduced fractions carry most of the signal, with oxidized pools smaller and less stable.`

### B2) Correct the late-dose clause (Totals rebound vs Reduced decline)

Replace `Totals partially rebounded at dose₆, whereas reduced pools continue downward once the UVB steps back up.` with  
`Totals partially rebound at dose₆, whereas reduced pools continue downward as UVB steps down (0.768→0.707 mW·cm⁻²) while UVA continues upward.`

### B3) Effect sizes vs. control (Δ and %)

Insert under “Dose-level trimmed means”:  
`Relative to control (dose₁), effect sizes were: Chrom Total 0.274→0.986 (Δ = 0.712 mg·gDW⁻¹; +260%), Chrom Reduced 0.740→1.790 (Δ = 1.050; +142%), DAD Total 0.520→1.844 (Δ = 1.324; +254%), DAD Reduced 1.050→2.093 (Δ = 1.043; +99%), DAD Oxidized −0.114→0.951 (Δ = 1.065).`

### B4) Peak-dose CIs

Append to each peak sentence (example):  
`…Chrom Total peaked at 0.986 mg·gDW⁻¹ at dose₄ (95% CI 0.668–1.197).`  
Repeat for Chrom Reduced dose₅, DAD Total dose₄, DAD Reduced dose₅, Chrom Oxidized dose₄, DAD Oxidized dose₄.

### B5) Stepwise transitions (include the third)

Replace the “only first jumps” sentence with:  
`Of the sequential dose-to-dose differences, three transitions are resolved at the 95% level: Chrom Reduced dose₁→₂ (Δ = 0.601; CI 0.235–1.296), DAD Reduced dose₁→₂ (Δ = 0.891; CI 0.391–1.760), and DAD Reduced dose₃→₄ (Δ = 0.319; CI 0.001–1.133).`

### B6) Clarify model labeling for trend fits

Prefix the trend paragraph with:  
`Using weighted least squares fits with UVA and UVA² (reporting the linear term β₁) and confirming with Kendall’s τ, we observe…`  
This header cures the ambiguity that previously invited the UVA-only confusion.

### B7) Axis regressions (UVA & UVB; UVB-dominant for Reduced)

Insert a short paragraph after the trend section:  
`Relation to UVA and UVB axes. Regressions of trimmed means on the physical dose axes show a UVB-dominant association for the reduced pools, with UVA still contributory: Chrom Reduced vs UVB slope 1.218 (95% CI 0.718–1.717; r = 0.900); vs UVA slope 0.390 (0.154–0.626; r = 0.748). DAD Reduced vs UVB slope 1.537 (0.697–2.378; r = 0.824); vs UVA slope 0.472 (0.133–0.812; r = 0.689). Totals are positively related to both axes (r ≈ 0.65–0.83) with wider slope CIs at this sample size.`

### B8) Cross-assay agreement (exact, not ≥ thresholds)

Replace any “≥ 0.95” phrasing with:  
`Cross-assay agreement (Chrom↔DAD) is high: Pearson r = 0.985 (Totals), 0.957 (Oxidized), 0.948 (Reduced). Deming fits (error ratio = 1) yield slopes ≈ 1.88 (Total), 2.03 (Oxidized), and 1.09 (Reduced) with intercept ≈ 0.33 mg·gDW⁻¹ for the Reduced fraction, indicating near-proportionality for Reduced and calibration offsets for the other pools.`

### B9) Linear UVA-only paragraph (if you decide to keep it)

Move to Methods-adjacent box and relabel:  
`Alternate analysis (axis-specific linear UVA, WLS) — Slopes 0.091–0.472 mg·gDW⁻¹·(mW·cm⁻²)⁻¹ summarize the UVA-only association without curvature. These complement, but do not replace, the canonical UVA + UVA² trend fits reported in the main text.`

## C) Placement & sourcing map (what to edit where)

- `dose_relationship_results.md`  
  - Lines 17–21: Replace/relocate UVA-only slopes (B1, B6, B7, B9).  
  - Late-dose sentence: apply fix (B2).  
  - Add Δ/% block (B3).  
  - Add peak CIs inline (B4).  
  - Stepwise transitions (B5).  
  - Axis regressions and cross-assay agreement (B7, B8).  
  - Negative oxidized control note retained.
- `summary.md`  
  - Ensure the same Δ/% effects, peak CIs, three resolved steps, UVB-dominant paragraph, and exact r values appear.  
  - Replace any “only the first step is significant” with the three-step sentence (B5).  
  - Replace any “r ≥ 0.95” with exact coefficients (B8).
- `plan.md` / `captions.md`  
  - Caption for the stepwise plot mentions three resolved transitions.  
  - Caption for the trends plot states UVA up & UVB down at dose₆.  
  - If you include axis scatter plots, captions state UVB-dominant for Reduced.
- (If present) `aggregated_reflectance/dose_metadata.md`  
  - Verify dose₆ note: UVA rises; UVB dips (no “UVB steps up” language).

## D) Quality gates (acceptance tests you can literally tick off)

### Numbers & units

- Every peak mean is immediately followed by a dose index and 95% CI with units.  
- Each metric has a control→peak sentence with Δ and % change.  
- Late-dose sentence says UVA↑, UVB↓ at dose₆ (no “UVB up”).  
- Three resolved steps named exactly: Chrom Red 1→2; DAD Red 1→2; DAD Red 3→4.

### Model hygiene

- Where β₁ is quoted, the sentence names the model: “UVA + UVA² (WLS).”  
- If UVA-only slopes are present anywhere, they’re segregated and labeled “Alternate analysis (linear UVA, WLS).”

### Axis & concordance

- Reduced vs UVB and vs UVA slopes with CIs and r are present.  
- Cross-assay Pearson r are 0.985/0.957/0.948 (no “≥”).  
- Deming slopes and intercept ≈ 0.33 mg·gDW⁻¹ for Reduced are reported.

### Artifacts & disclaimers

- Negative DAD oxidized control is acknowledged with CI and baseline-subtraction rationale.

### Narrative coherence

- The “reduced carries the signal; oxidized smaller/unstable” sentence appears in the opening Results paragraph.  
- No internal contradictions remain between prose and the dose table.

## E) Optional polish (high-leverage, low-effort)

- **Figure callouts:** At the start of each subsection add “(Fig. X)” so readers see the relevant plot when the claim lands.  
- **Panel scaling sanity:** For the stepwise figure, keep consistent y-limits across metrics so the sign pattern + + + − + and + − + + − is visually comparable.  
- **Notation pass:** Use dose₁–dose₆ consistently; units as mg·gDW⁻¹, mW·cm⁻²; italicize p only, not CI.  
- **Method breadcrumb:** One clause early in Results — “20% trimmed means; 2,000-draw percentile bootstrap; n = 5 per dose” — keeps readers oriented without sending them to Methods.

## Why this works

You’re repairing truth-conditions (numbers and models), narrative scaffolding (dose geometry → peaks/uncertainty → steps → trends → axes → concordance), and reader calibration (effect sizes and CIs). That sequence matches how minds internalize causal structure: from what happened (dose), to how big (Δ/% + CIs), to how it changes (steps/trends), to what drives it (axes), to can we trust it (cross-assay). When you’ve applied the patches, read the Results out loud once, listening for the rhythm: geometry → magnitude → uncertainty → shape → driver → agreement → caveat. If it flows in that order and every quantitative claim is visibly tethered to a source number, you’ve turned the statistical journey into a map the reader can actually follow.

