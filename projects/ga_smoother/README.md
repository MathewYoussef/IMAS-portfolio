# GA Smoother (spectral pre-processing)

Genetic search across smoothing families (Savitzky–Golay, moving median/average, wavelet denoising) to pick parameters that **maximize SNR** and **minimize CV** in a UVA band (λ 310–410 nm). Built on cyanobacteria/Nostoc spectra where a diagnostic dip is expected; windows adapt to your wavelength grid.

## How to run (local)
1) Put your CSVs under a folder with a subfolder `norm_out/` containing files like `*_CR.csv`.
2) Run with an environment variable pointing at that data directory, e.g.:

GA_SMOOTHER_DATA_DIR="/path/to/your/experiment" \
python projects/ga_smoother/ga_smoother.py

## Outputs
- norm_out/nostoc_CR_smoothed_OPT.csv (best smoothed set)
- norm_out/GA_smoother_diagnostics.png (quick plot)
- Console prints: SNR, CV, Δλ metrics

## Notes
- Default target band: 310–410 nm (hard-coded in this version).
- Dependencies: numpy, pandas, scipy, pywavelets, pygad, matplotlib.
- Rationale: let a GA compete across methods/params; select preprocessing by data, not hunch.
