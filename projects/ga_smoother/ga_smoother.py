import os
#!/usr/bin/env python3
# ------------------------------------------------------------------
#  GA_for_smoother.py   •  GA smoother optimiser (λ 310–410 nm only)
# ------------------------------------------------------------------
import pandas as pd, numpy as np, pywt, pygad
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.special import wofz
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", message="Mean of empty slice")

# ================================================================
# 1. AUTO-DISCOVER & LOAD  (no more ALL_CR.csv needed)
# ================================================================
ROOT = Path(os.environ.get("GA_SMOOTHER_DATA_DIR", Path.cwd()))  # set via env or uses CWD
CR_DIR = ROOT / "norm_out"

from typing import List
def _gather_csvs(base: Path) -> list[Path]:
    candidates = []
    cr_sub = base / "CR"
    if cr_sub.exists():
        candidates += sorted(cr_sub.glob("*.csv"))
    pats = ["*_CR.csv", "*_norm_*.csv", "*_norm.csv", "*_UVrange.csv", "*.csv"]
    for pat in pats:
        candidates += [pp for pp in base.glob(pat) if pp.is_file()]
    uniq, seen = [], set()
    for pp in candidates:
        name = pp.name
        # exclude our own outputs
        if any(x in name for x in ("_smoothed_", "_OPT")):
            continue
        if name.startswith("."):  # hidden files
            continue
        key = pp.resolve()
        if key in seen:
            continue
        seen.add(key)
        uniq.append(pp)
    return uniq

CR_DIR.mkdir(parents=True, exist_ok=True)

csv_list = _gather_csvs(CR_DIR)
print(f"[ga] found {len(csv_list)} csv(s) under {CR_DIR}")
if not csv_list:
    raise SystemExit(f"❌  No suitable CSVs found in {CR_DIR}. Expecting *_CR.csv or *_norm*.csv; found none.")

frames = []
for fp in csv_list:
    tmp = pd.read_csv(fp)                       # expects columns: wavelength, reflectance_CR
    # ---- rename & tag -------------------------------------------------------
    tmp = tmp.rename(columns={"reflectance_CR": "cr"})
    tmp["file_name"] = fp.stem.replace("_CR", "")
    frames.append(tmp)

df = pd.concat(frames, ignore_index=True)

# keep only wavelengths 310-410 nm -------------------------------------------
ROI_MIN, ROI_MAX = 310, 410
df = df[(df["wavelength"] >= ROI_MIN) & (df["wavelength"] <= ROI_MAX)]

# ----- pivot to spectra-matrix ----------------------------------------------
spec_mat = (df.pivot_table(index="file_name",
                           columns="wavelength",
                           values="cr")
              .sort_index(axis=1))             # λ ascending

λ       = spec_mat.columns.to_numpy(dtype=float)
spectra = spec_mat.to_numpy()

# (optional) metadata ---------------------------------------------------------
meta = (df[["file_name"]]
        .drop_duplicates("file_name")
        .set_index("file_name")
        .reindex(spec_mat.index))

# ================================================================
# 2. SMOOTHING OPERATIONS (unchanged)
# ================================================================
def sg(y, order, window):
    return savgol_filter(y, window_length=window, polyorder=order)

def wav(y, family, level):
    coeffs = pywt.wavedec(y, family, mode="periodization", level=level)
    σ      = np.median(np.abs(coeffs[-1])) / 0.6745
    uth    = σ * np.sqrt(2*np.log(len(y)))
    coeffs[1:] = [pywt.threshold(c, value=uth, mode='hard') for c in coeffs[1:]]
    return pywt.waverec(coeffs, family, mode="periodization")[:len(y)]

def cascade(y, sg_par=None, wav_par=None, order='sg->wav'):
    if order == 'sg->wav':
        y1 = sg(y, *sg_par)  if sg_par  else y
        return wav(y1, *wav_par) if wav_par else y1
    else:
        y1 = wav(y, *wav_par) if wav_par else y
        return sg(y1, *sg_par)  if sg_par  else y1

# ================================================================
# 3. METRICS (windows auto-adapt to λ array)
# ================================================================
mask_lo, mask_hi = (λ>=335)&(λ<=345), (λ>=400)&(λ<=410)
shoulder_idx     = mask_lo | mask_hi
dip_idx          = (λ>=355)&(λ<=395)

def voigt(x, amp, cen, sigma, gamma):
    z = ((x-cen) + 1j*gamma) / (sigma*np.sqrt(2))
    with np.errstate(over="ignore", invalid="ignore"):
        y = amp * np.real(wofz(z)) / (sigma*np.sqrt(2*np.pi))
    return np.nan_to_num(y, nan=0.0, posinf=0.0, neginf=0.0)
def eval_metrics(y_smooth):
    shoulder = y_smooth[:, shoulder_idx]
    snr      = shoulder.mean() / shoulder.std()

    dip_depth = 1 - y_smooth[:, dip_idx].min(axis=1)
    cv        = dip_depth.std() / dip_depth.mean()

    y_mean = y_smooth.mean(axis=0)
    try:
        popt, _ = curve_fit(voigt,
                            λ[dip_idx],
                            1 - y_mean[dip_idx],
                            p0=[0.4, 370, 5, 2],
                            maxfev=4000)
        delta = abs(popt[1] - 370)
    except RuntimeError:
        dip_centre = λ[dip_idx][np.argmin(y_mean[dip_idx])]
        delta      = abs(dip_centre - 370)
    return snr, cv, delta

def smooth_complexity(which, sg_win, wav_lvl):
    return sg_win if which==0 else wav_lvl*10 if which==1 else sg_win+wav_lvl*10

def rms_outside(y_raw, y_smooth):
    keep = shoulder_idx & ~dip_idx
    return np.sqrt(np.mean((y_raw[:, keep] - y_smooth[:, keep])**2))

def z(x): return x                    # simple proxy

families = ['sym4','sym6','sym8','coif3','db4','db5']

# ================================================================
# 4. FITNESS FUNCTION
# ================================================================
def fitness_func(_, sol, __):
    which, sg_ord, sg_win, fam_idx, lvl, order_flag = sol
    which, sg_ord, sg_win = map(lambda v: int(round(v)), (which, sg_ord, sg_win))
    fam_idx, lvl          = int(round(fam_idx)), int(round(lvl))
    order_flag            = float(order_flag)

    if which in (0,2) and ((sg_win%2==0) or (sg_win<=sg_ord)):        return -1e9
    if which in (1,2):
        fam = families[fam_idx]
        if lvl > pywt.dwt_max_level(len(λ), pywt.Wavelet(fam).dec_len): return -1e9

    sg_par  = (sg_ord, sg_win)            if which in (0,2) else None
    wav_par = (families[fam_idx], lvl)    if which in (1,2) else None

    y_sm = np.apply_along_axis(cascade, 1, spectra,
                               sg_par=sg_par, wav_par=wav_par,
                               order='sg->wav' if order_flag<0.5 else 'wav->sg')

    snr, cv, delta = eval_metrics(y_sm)
    rms  = rms_outside(spectra, y_sm)
    comp = smooth_complexity(which, sg_win, lvl)

    if not np.isfinite(snr+cv+delta): return -1e9
    return z(snr) - z(cv) - z(delta) - 0.7*z(rms) - 0.4*z(comp)

# ================================================================
# 5. GA CONFIGURATION
# ================================================================
gene_space = [
    [0,1,2],                # smoother type
    [2,3,4],                # SG order
    [5,7,9,11,13],          # SG window (odd)
    list(range(len(families))),  # wavelet family idx
    [2,3,4,5,6],            # wavelet level
    [0,1]                   # cascade order flag
]

ga = pygad.GA(num_generations        = 100,
              sol_per_pop            = 64,
              num_genes              = len(gene_space),
              parent_selection_type  = "tournament",
              K_tournament           = 4,
              num_parents_mating     = 24,
              keep_parents           = 8,
              crossover_type         = "uniform",
              crossover_probability  = 0.9,
              mutation_type          = "adaptive",
              mutation_percent_genes = [30,10],
              gene_space             = gene_space,
              fitness_func           = fitness_func,
              stop_criteria          = ["saturate_20"],
              random_seed            = 42)

ga.run()

# ================================================================
# 6. REPORT & EXPORT
# ================================================================
best_genome, best_fitness, _ = ga.best_solution()
print("\nBest genome :", best_genome)
print("Fitness      :", best_fitness)

which, sg_ord, sg_win, fam_idx, lvl, order_flag = best_genome
sg_ord, sg_win, fam_idx, lvl = map(int, map(round, (sg_ord, sg_win, fam_idx, lvl)))

sg_par  = (sg_ord, sg_win)                if which in (0,2) else None
wav_par = (families[fam_idx], lvl)        if which in (1,2) else None

opt_sm = np.apply_along_axis(cascade, 1, spectra,
                             sg_par=sg_par, wav_par=wav_par,
                             order='sg->wav' if order_flag<0.5 else 'wav->sg')

out_csv = ROOT / "norm_out" / "nostoc_CR_smoothed_OPT.csv"
pd.DataFrame(opt_sm, columns=λ).to_csv(out_csv, index=False)

snr, cv, delta = eval_metrics(opt_sm)
print(f"SNR={snr:.3f}   CV={cv:.4f}   Δλ={delta:.2f} nm  (λ 310–410 nm)")

plt.figure()
plt.plot(ga.best_solutions_fitness, marker='o')
plt.title("GA best fitness per generation")
plt.xlabel("Generation"); plt.ylabel("Fitness"); plt.grid(True)
plt.tight_layout(); plt.show()