import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

def load_dataset(path: Path):
    if path.suffix.lower() in [".nc", ".nc4", ".cdf"]:
        return xr.open_dataset(path)
    elif path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
        return df
    else:
        raise ValueError(f"Unsupported file type: {path.suffix}")

def simple_qc_flags(series: pd.Series, z=4.0):
    """Return boolean mask of 'good' values using z-score outlier filter."""
    vals = series.to_numpy(dtype=float)
    mu, sigma = np.nanmean(vals), np.nanstd(vals)
    if sigma == 0 or np.isnan(sigma):
        return np.isfinite(vals)
    zscores = (vals - mu) / sigma
    return np.isfinite(vals) & (np.abs(zscores) < z)

def main():
    p = argparse.ArgumentParser(description="Demo: ingest, QC, and plot a small time series.")
    p.add_argument("--input", type=str, default="data/sample/sample.csv",
                   help="Path to a tiny sample file (.nc or .csv). Default uses CSV fallback.")
    p.add_argument("--var", type=str, default="value", help="Variable/column name to process.")
    p.add_argument("--time", type=str, default="time", help="Time dimension/column name.")
    p.add_argument("--out", type=str, default="demo_output.png", help="Output plot path.")
    args = p.parse_args()

    path = Path(args.input)
    if not path.exists():
        # create a tiny CSV time series on the fly (offline-safe demo)
        Path("data/sample").mkdir(parents=True, exist_ok=True)
        t = pd.date_range("2020-01-01", periods=200, freq="H")
        y = np.sin(np.linspace(0, 6*np.pi, len(t))) + 0.1*np.random.randn(len(t))
        # sprinkle a few outliers
        y[::53] += 3.5
        df = pd.DataFrame({ "time": t, "value": y })
        df.to_csv("data/sample/sample.csv", index=False)
        path = Path("data/sample/sample.csv")

    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path, parse_dates=[args.time])
        df = df.sort_values(args.time)
        mask = simple_qc_flags(df[args.var], z=4.0)
        df["good"] = mask
        fig = plt.figure(figsize=(10,4))
        plt.plot(df[args.time], df[args.var], alpha=0.3, label="raw")
        plt.plot(df.loc[mask, args.time], df.loc[mask, args.var], linewidth=1.5, label="QC-pass")
        plt.legend()
        plt.title("Demo ingest → QC → plot")
        plt.tight_layout()
        fig.savefig(args.out, dpi=150)
        print(f"Wrote {args.out} | rows: {len(df)} | QC-pass: {mask.sum()}")
    else:
        ds = load_dataset(path)
        if args.time not in ds.dims and args.time not in ds:
            print(f"Time key '{args.time}' not found in dataset.", file=sys.stderr)
            sys.exit(1)
        if args.var not in ds and args.var not in ds.data_vars:
            print(f"Var '{args.var}' not found in dataset.", file=sys.stderr)
            sys.exit(1)
        series = ds[args.var].to_pandas()
        if isinstance(series, pd.DataFrame):
            series = series.iloc[:,0]
        tindex = series.index
        mask = simple_qc_flags(series, z=4.0)
        fig = plt.figure(figsize=(10,4))
        plt.plot(tindex, series.values, alpha=0.3, label="raw")
        plt.plot(tindex[mask], series.values[mask], linewidth=1.5, label="QC-pass")
        plt.legend()
        plt.title("Demo ingest → QC → plot (NetCDF)")
        plt.tight_layout()
        fig.savefig(args.out, dpi=150)
        print(f"Wrote {args.out} | len: {len(series)} | QC-pass: {int(mask.sum())}")

if __name__ == "__main__":
    main()
