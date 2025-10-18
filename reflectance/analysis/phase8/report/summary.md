# Reflectance vs. Scytonemin Concentration: Summary

## Key Metrics
- **Lowest LOOCV RMSE**: 0.07 (kind=Sigma, modality=chrom, component=oxidized, band=peak_band_auto, model=ridge).
- **Highest LOOCV RMSE**: 1.90 (kind=6Oclock, modality=dad, component=total, band=320_480nm, model=ridge).
- **Best distance-correlation p-value**: 0.138 (kind=6Oclock, modality=chrom).
- **Average in-sample RMSE (Phase 6 baseline)**: 2.79e-12 (orders of magnitude smaller than LOOCV errors).
- **Multiblock PLS joint components**: 1 per modality (no warnings).

## Interpretation
- Phase 3 permutation p-values stay >0.9, so linear correlations are not statistically convincing.
- Phase 4 dependence tests (dCor/RV) yield p-values between ~0.14 and ~0.84; none reach conventional significance.
- Multiblock fusion adds only a single shared component, reflecting the weak dependence signal.
- Ridge models with advanced features yield LOOCV RMSEs spanning 0.07–1.90 depending on kind, modality, and component; oxidized Σ fits lean on peak height, AUC, FWHM, and curvature around the fixed 406 nm band, whereas reduced forms still respond weakly even with the broader 320–480 nm descriptors.
- Under current sampling, reflectance signatures alone cannot reliably predict scytonemin concentration; more data or regularised methods with validation remain essential.

## Next Steps
- Expand dose/replicate counts to improve statistical power.
- Explore regularised or nonlinear models only after increasing sample size.
- Replace the stacked PLS approximation with a true O2PLS implementation once data volume allows.
