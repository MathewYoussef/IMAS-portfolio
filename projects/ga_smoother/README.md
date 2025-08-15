# GA Smoother (spectral pre-processing)

Genetic search across smoothing families (Savitzky–Golay, moving median/average, wavelet denoising) to pick parameters that **maximize SNR** and **minimize CV** in a UVA band (λ 310–410 nm). Built on cyanobacteria/Nostoc spectra where a diagnostic dip is expected; windows adapt to your wavelength grid.

## Install
```bash
python -m pip install -r projects/ga_smoother/requirements.txt
Run (your data)

Your experiment folder must contain a subfolder norm_out/ with CSVs like *_CR.csv or *_norm*.csv.
GA_SMOOTHER_DATA_DIR="/path/to/your/experiment" \
python projects/ga_smoother/run_and_describe.py
GA_SMOOTHER_DATA_DIR="projects/ga_smoother/demo" \
python projects/ga_smoother/run_and_describe.py
GA_SMOOTHER_DATA_DIR="." \
GA_SMOOTHER_DECODE_ONLY='[0,3,5,1,4,1]' \
python projects/ga_smoother/run_and_describe.py
Outputs
•norm_out/nostoc_CR_smoothed_OPT.csv — best smoothed set
•norm_out/GA_smoother_diagnostics.png — quick plot
•norm_out/GA_smoother_best.json — chosen family + params + raw genome (+ mapping details)
•Console prints: SNR, CV, Δλ, and a human summary (“Chosen smoother: …”)

Notes
•Default target band: 310–410 nm (fixed in this version).
•Dependencies: numpy, pandas, scipy, pywavelets, pygad, matplotlib.
•Rationale: let a GA compete across methods/params; pick preprocessing by data, not hunch.

## Try the bundled demo

```bash
GA_SMOOTHER_DATA_DIR="projects/ga_smoother/demo" \
python projects/ga_smoother/run_and_describe.py
```
