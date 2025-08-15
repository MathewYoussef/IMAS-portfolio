import os, re, json, subprocess, sys
from pathlib import Path

# ---- Mapping of genome -> params (edit if your GA uses different search spaces) ----
FAMILIES = ["savitzky_golay", "median", "moving_average", "wavelet"]
WIN_OPTS  = [5,7,9,11,13,15,17]             # smoothing windows (odd)
POLY_OPTS = [2,3,4,5]                        # SG polyorder
DER_OPTS  = [0,1,2]                          # SG derivative order
WVL_NAMES = ["db2","db3","sym3","coif1","coif2"]
WVL_LVLS  = [1,2,3,4,5]
THRESH    = ["soft","hard"]

def _idx_detail(opts, i):
    raw = int(round(float(i)))
    mod = raw % len(opts)
    return opts[mod], raw, mod

def decode_with_details(genome):
    g = [float(x) for x in genome]
    fam_raw = int(round(g[0]))
    fam_mod = fam_raw % len(FAMILIES)
    fam_name = FAMILIES[fam_mod]
    details = {"family_idx_raw": fam_raw, "family_idx_mod": fam_mod, "family_name": fam_name}

    if fam_name == "savitzky_golay":
        win, wr, wm = _idx_detail(WIN_OPTS,  g[1])
        po,  pr, pm = _idx_detail(POLY_OPTS, g[2])
        der, dr, dm = _idx_detail(DER_OPTS,  g[3])
        decoded = {"family": fam_name, "params": {"window": win, "polyorder": po, "derivative": der}}
        details |= {"window":{"raw":wr,"mod":wm,"value":win},
                    "polyorder":{"raw":pr,"mod":pm,"value":po},
                    "derivative":{"raw":dr,"mod":dm,"value":der}}
    elif fam_name in {"median","moving_average"}:
        win, wr, wm = _idx_detail(WIN_OPTS, g[1])
        decoded = {"family": fam_name, "params": {"window": win}}
        details |= {"window":{"raw":wr,"mod":wm,"value":win}}
    else:  # wavelet
        wv, wr, wm = _idx_detail(WVL_NAMES, g[1])
        lv, lr, lm = _idx_detail(WVL_LVLS,  g[2])
        th, tr, tm = _idx_detail(THRESH,    g[3])
        decoded = {"family": fam_name, "params": {"wavelet": wv, "level": lv, "mode": th}}
        details |= {"wavelet":{"raw":wr,"mod":wm,"value":wv},
                    "level":{"raw":lr,"mod":lm,"value":lv},
                    "mode":{"raw":tr,"mod":tm,"value":th}}
    return decoded, details

def print_human(decoded, details, data_root, genes):
    print("\n=== HUMAN SUMMARY ===")
    print("RAW GENOME:", [float(x) for x in genes])
    print("DECODE DETAILS:", json.dumps(details, indent=2))
    fam = decoded["family"]
    params = ", ".join(f"{k}={v}" for k,v in decoded["params"].items())
    print(f"Chosen smoother: {fam}  |  {params}")
    out_dir = Path(data_root) / "norm_out"
    out_dir.mkdir(parents=True, exist_ok=True)
    info = {"genome": [float(x) for x in genes],
            "selection": decoded,
            "details": details,
            "note": "Mapping is heuristic; adjust arrays if your GA space differs."}
    (out_dir / "GA_smoother_best.json").write_text(json.dumps(info, indent=2))
    print(f"Wrote: {out_dir/'GA_smoother_best.json'}")
    print("======================")

def main():
    env = os.environ.copy()
    data_root = env.get("GA_SMOOTHER_DATA_DIR", "")
    if not data_root:
        print("ERROR: GA_SMOOTHER_DATA_DIR not set.")
        sys.exit(2)

    # Optional quick test: decode-only mode (skip running GA)
    if "GA_SMOOTHER_DECODE_ONLY" in env and env["GA_SMOOTHER_DECODE_ONLY"].strip():
        genes = json.loads(env["GA_SMOOTHER_DECODE_ONLY"])
        decoded, details = decode_with_details(genes)
        print_human(decoded, details, data_root, genes)
        return

    # Otherwise run the GA script and capture stdout
    proc = subprocess.run(
        [sys.executable, "projects/ga_smoother/ga_smoother.py"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env
    )
    print(proc.stdout)  # show GA stdout/stderr

    m = re.search(r"Best genome\s*:\s*\[([^\]]+)\]", proc.stdout)
    if not m:
        print("Could not locate 'Best genome' line. No decode produced.")
        sys.exit(1)

    genes = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", "[" + m.group(1) + "]")]
    decoded, details = decode_with_details(genes)
    print_human(decoded, details, data_root, genes)

if __name__ == "__main__":
    main()
