# Introduction to Reflectance vs Concentration

## Observed Measurements
- Dry-weight-normalized data show that chromatogram and diode array totals do not equal the oxidized plus reduced fractions.
- The offset is consistent across samples, so the effect is systematic rather than random instrument noise.

## Why oxidized + reduced ≠ total
### Diode array spectra
- The DAD integrates all absorbing molecules eluting between 2.4 and 3.1 minutes across the monitored wavelengths.
- Any co-extracted pigments or metabolites that absorb within that window add to the apparent total signal even if they do not fall cleanly into the oxidized or reduced scytonemin reference channels.

### Chromatogram areas
- The chromatogram signal represents the area under the curve at 384 nm across the retention span that includes both the reduced (~1.10 min) and oxidized (~1.36 min) peaks.
- Molecules co-eluting with scytonemin and absorbing at 384 nm contribute to the area, so the integrated total reflects more than just scytonemin.

## Biological contributions
- Environmental samples contain diverse pigments, phenolics, and degradation products that co-elute with scytonemin.
- These compounds can shift the oxidized/reduced balance or introduce additional absorbance, creating totals that exceed the sum of the scytonemin-specific fractions.
- Sample matrices are not pure; treating them as if they matched the scytonemin control standard inherently overshoots concentrations compared with a pure reference.

## Analytical implications
- Totals are best interpreted as “scytonemin-equivalent” absorbance within the specified retention window, rather than a strict oxidized + reduced mass balance.
- Chromatogram and DAD concentrations are not interchangeable: paired tests show significant offsets (e.g., total mean difference ~0.70 mg gDW⁻¹) and DAD routinely exceeds chromatogram estimates.
- Despite the offset, both measurement paths track each sample in near lockstep (Pearson r ≈0.95–0.96 across total, oxidized, reduced fractions), so relative trends with respect to sample ID are preserved.
- Given the additional spectral detail available from the diode array detector, downstream modeling will prioritize the DAD dataset for feature extraction and linking reflectance to concentration.

## Canonical data reference
- The harmonised working tables (reflectance summaries, DAD statistics, latent concentrations, and comparison outputs) now live under `canonical_dataset/` and are keyed consistently by `dose_id`, `uva_mw_cm2`, and `uvb_mw_cm2`.
- Run `python tests/verify_canonical_dataset.py` after any regeneration to confirm the canonical bundle still matches the archived source within floating-point tolerance.
- Legacy aggregations that pre-date the harmonisation have been archived under `archive/aggregated_reflectance/` and should only be consulted for historical auditing; all new analysis should source inputs from the canonical bundle.
- For provenance, consult `canonical_dataset/README.json`. It captures the build command (`python build_canonical_dataset.py --output canonical_dataset --archive`) and reiterates where the archived source tables reside, so audits can be reproducibly traced back to the original aggregations.
