# Phase 9 Equation Reference

This document summarises the curated ridge equations derived from Phase 9 features. 
Each entry is sourced from `analysis/phase9/tables/advanced_baseline_summary.csv` and flagged as reporting or QA-only.

## chrom_oxidized – Sigma / peak_band_auto (REPORTING)


| Field | Value |
|---|---|
| Modality | chrom |
| Component | oxidized |
| Estimator | ridge (ridge) |
| LOOCV RMSE | 0.074790 |
| LOOCV MAE | 0.052151 |
| In-sample RMSE | 0.036600 |
| In-sample MAE | 0.024526 |
| In-sample R² | 0.958856 |


**Feature Columns & Coefficients**


| Feature | Weight |
|---|---|
| peak_reflectance | 0.000063 |
| peak_wavelength_nm | -0.049194 |
| auc | 0.001685 |
| continuum_depth | -0.000258 |
| fwhm_nm | -0.032634 |
| max_curvature | -0.000385 |
| fixed_peak_reflectance | 0.000000 |
| fixed_continuum_depth | 0.000000 |
| Intercept | 33.670264 |

**Notes**
- Source row: `analysis/phase9/tables/advanced_baseline_summary.csv` (target=chrom_oxidized, kind=Sigma, band=peak_band_auto).
- Features computed from `analysis/phase9/tables/advanced_band_features.csv`; no smoothing applied, fixed peaks captured via `fixed_peak_*` columns where applicable.
- Concentration targets drawn from `canonical_dataset/dose_level_canonical_summary.csv`.


## dad_oxidized – Sigma / peak_band_auto (REPORTING)


| Field | Value |
|---|---|
| Modality | dad |
| Component | oxidized |
| Estimator | ridge (ridge) |
| LOOCV RMSE | 0.166129 |
| LOOCV MAE | 0.134155 |
| In-sample RMSE | 0.074774 |
| In-sample MAE | 0.055429 |
| In-sample R² | 0.956211 |


**Feature Columns & Coefficients**


| Feature | Weight |
|---|---|
| peak_reflectance | 0.001918 |
| peak_wavelength_nm | -0.055113 |
| auc | 0.011666 |
| continuum_depth | 0.000482 |
| fwhm_nm | -0.077675 |
| max_curvature | -0.001232 |
| fixed_peak_reflectance | 0.000000 |
| fixed_continuum_depth | 0.000000 |
| Intercept | 42.389713 |

**Notes**
- Source row: `analysis/phase9/tables/advanced_baseline_summary.csv` (target=dad_oxidized, kind=Sigma, band=peak_band_auto).
- Features computed from `analysis/phase9/tables/advanced_band_features.csv`; no smoothing applied, fixed peaks captured via `fixed_peak_*` columns where applicable.
- Concentration targets drawn from `canonical_dataset/dose_level_canonical_summary.csv`.


## chrom_total – 12Oclock / 320_480nm (REPORTING)


| Field | Value |
|---|---|
| Modality | chrom |
| Component | total |
| Estimator | ridge (ridge) |
| LOOCV RMSE | 0.302980 |
| LOOCV MAE | 0.221067 |
| In-sample RMSE | 0.165875 |
| In-sample MAE | 0.114924 |
| In-sample R² | 0.436234 |


**Feature Columns & Coefficients**


| Feature | Weight |
|---|---|
| peak_reflectance | 0.002673 |
| peak_wavelength_nm | 0.000000 |
| auc | -0.004233 |
| continuum_depth | 0.001069 |
| fwhm_nm | -0.020796 |
| max_curvature | -0.000028 |
| fixed_peak_reflectance | 0.000000 |
| fixed_continuum_depth | 0.000000 |
| Intercept | 2.001393 |

**Notes**
- Source row: `analysis/phase9/tables/advanced_baseline_summary.csv` (target=chrom_total, kind=12Oclock, band=320_480nm).
- Features computed from `analysis/phase9/tables/advanced_band_features.csv`; no smoothing applied, fixed peaks captured via `fixed_peak_*` columns where applicable.
- Concentration targets drawn from `canonical_dataset/dose_level_canonical_summary.csv`.


## dad_total – Sigma / peak_band_auto (REPORTING)


| Field | Value |
|---|---|
| Modality | dad |
| Component | total |
| Estimator | ridge (ridge) |
| LOOCV RMSE | 0.553727 |
| LOOCV MAE | 0.522883 |
| In-sample RMSE | 0.199569 |
| In-sample MAE | 0.198020 |
| In-sample R² | 0.764162 |


**Feature Columns & Coefficients**


| Feature | Weight |
|---|---|
| peak_reflectance | 0.010757 |
| peak_wavelength_nm | -0.001596 |
| auc | 0.014480 |
| continuum_depth | 0.005636 |
| fwhm_nm | -0.088410 |
| max_curvature | -0.004504 |
| fixed_peak_reflectance | 0.000000 |
| fixed_continuum_depth | 0.000000 |
| Intercept | 12.577186 |

**Notes**
- Source row: `analysis/phase9/tables/advanced_baseline_summary.csv` (target=dad_total, kind=Sigma, band=peak_band_auto).
- Features computed from `analysis/phase9/tables/advanced_band_features.csv`; no smoothing applied, fixed peaks captured via `fixed_peak_*` columns where applicable.
- Concentration targets drawn from `canonical_dataset/dose_level_canonical_summary.csv`.


## latent_total – Sigma / peak_band_auto (REPORTING)


| Field | Value |
|---|---|
| Modality | latent |
| Component | total |
| Estimator | ridge (ridge) |
| LOOCV RMSE | 0.289561 |
| LOOCV MAE | 0.273152 |
| In-sample RMSE | 0.104712 |
| In-sample MAE | 0.103560 |
| In-sample R² | 0.750894 |


**Feature Columns & Coefficients**


| Feature | Weight |
|---|---|
| peak_reflectance | 0.005399 |
| peak_wavelength_nm | -0.018848 |
| auc | 0.009366 |
| continuum_depth | 0.002845 |
| fwhm_nm | -0.042162 |
| max_curvature | -0.002355 |
| fixed_peak_reflectance | 0.000000 |
| fixed_continuum_depth | 0.000000 |
| Intercept | 16.631322 |

**Notes**
- Source row: `analysis/phase9/tables/advanced_baseline_summary.csv` (target=latent_total, kind=Sigma, band=peak_band_auto).
- Features computed from `analysis/phase9/tables/advanced_band_features.csv`; no smoothing applied, fixed peaks captured via `fixed_peak_*` columns where applicable.
- Concentration targets drawn from `canonical_dataset/dose_level_canonical_summary.csv`.


## chrom_reduced – Delta / 320_480nm (QA)


| Field | Value |
|---|---|
| Modality | chrom |
| Component | reduced |
| Estimator | ridge (ridge) |
| LOOCV RMSE | 0.440106 |
| LOOCV MAE | 0.382564 |
| In-sample RMSE | 0.144916 |
| In-sample MAE | 0.106254 |
| In-sample R² | 0.797925 |


**Feature Columns & Coefficients**


| Feature | Weight |
|---|---|
| peak_reflectance | 0.001506 |
| peak_wavelength_nm | -0.001220 |
| auc | 0.025308 |
| continuum_depth | 0.000324 |
| fwhm_nm | -0.003901 |
| max_curvature | 0.000509 |
| fixed_peak_reflectance | 0.000000 |
| fixed_continuum_depth | 0.000000 |
| Intercept | 1.985856 |

**Notes**
- Source row: `analysis/phase9/tables/advanced_baseline_summary.csv` (target=chrom_reduced, kind=Delta, band=320_480nm).
- Features computed from `analysis/phase9/tables/advanced_band_features.csv`; no smoothing applied, fixed peaks captured via `fixed_peak_*` columns where applicable.
- Concentration targets drawn from `canonical_dataset/dose_level_canonical_summary.csv`.


## dad_reduced – 6Oclock / peak_band_auto (QA)


| Field | Value |
|---|---|
| Modality | dad |
| Component | reduced |
| Estimator | ridge (ridge) |
| LOOCV RMSE | 0.384400 |
| LOOCV MAE | 0.305269 |
| In-sample RMSE | 0.144206 |
| In-sample MAE | 0.119815 |
| In-sample R² | 0.831377 |


**Feature Columns & Coefficients**


| Feature | Weight |
|---|---|
| peak_reflectance | 0.001344 |
| peak_wavelength_nm | -0.039653 |
| auc | 0.006902 |
| continuum_depth | -0.001260 |
| fwhm_nm | -0.081786 |
| max_curvature | -0.000772 |
| fixed_peak_reflectance | 0.000000 |
| fixed_continuum_depth | 0.000000 |
| Intercept | 35.211533 |

**Notes**
- Source row: `analysis/phase9/tables/advanced_baseline_summary.csv` (target=dad_reduced, kind=6Oclock, band=peak_band_auto).
- Features computed from `analysis/phase9/tables/advanced_band_features.csv`; no smoothing applied, fixed peaks captured via `fixed_peak_*` columns where applicable.
- Concentration targets drawn from `canonical_dataset/dose_level_canonical_summary.csv`.

