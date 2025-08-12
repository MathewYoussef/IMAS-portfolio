import argparse
from pathlib import Path
import pandas as pd
import xarray as xr

def scan_file(path: Path):
    info = {
        "file": str(path),
        "size_kb": round(path.stat().st_size / 1024, 1),
        "variables": None,
        "time_start": None,
        "time_end": None,
    }
    try:
        if path.suffix.lower() in [".nc", ".nc4", ".cdf"]:
            ds = xr.open_dataset(path, decode_times=False)
            info["variables"] = ",".join(list(ds.data_vars))
            if "time" in ds:
                tvals = ds["time"].values
                if len(tvals) > 0:
                    info["time_start"] = str(tvals[0])
                    info["time_end"] = str(tvals[-1])
        elif path.suffix.lower() == ".csv":
            df = pd.read_csv(path)
            info["variables"] = ",".join(df.columns)
            if "time" in df.columns:
                tvals = pd.to_datetime(df["time"], errors="coerce").dropna()
                if len(tvals) > 0:
                    info["time_start"] = str(tvals.min())
                    info["time_end"] = str(tvals.max())
    except Exception as e:
        info["variables"] = f"ERROR: {e}"
    return info

def main():
    p = argparse.ArgumentParser(description="Catalogue .csv/.nc files in a directory")
    p.add_argument("--root", type=str, default="data/sample", help="Root folder to scan")
    p.add_argument("--out", type=str, default="catalog.csv", help="Output CSV path")
    args = p.parse_args()

    root = Path(args.root)
    files = [f for f in root.rglob("*") if f.suffix.lower() in [".csv", ".nc", ".nc4", ".cdf"]]
    rows = [scan_file(f) for f in files]
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print(f"Wrote {args.out} with {len(df)} entries.")

if __name__ == "__main__":
    main()
