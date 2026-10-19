"""Plot Btotal vs time for each Physics Toolbox magnetometer CSV.

Usage (from the repo root):
    python analysis/plot_btotal.py                     # first file through 12-02-52
    python analysis/plot_btotal.py --start 12-36-58 --end 13-26-39
"""
import argparse
import filecmp
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "magnetic_field_data"
FIG_DIR = REPO / "figures"


def stamp(path):
    """'magnetometer_2026-09-19_11-09-51.csv' -> '11-09-51'."""
    return path.stem.split("_")[-1].split(" ")[0]


def select_files(start, end):
    files = sorted(DATA_DIR.glob("magnetometer_*.csv"))
    picked = []
    for f in files:
        if " (1)" in f.name:
            original = f.with_name(f.name.replace(" (1)", ""))
            if original.exists() and filecmp.cmp(f, original, shallow=False):
                continue  # exact duplicate download
        s = stamp(f)
        if (start is None or s >= start) and (end is None or s <= end):
            picked.append(f)
    return picked


def load(path):
    # First line is "# sensor:magnetometer"; columns: time,Bx,By,Bz,Btotal (µT)
    return pd.read_csv(path, comment="#")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--start", default=None, help="first timestamp, HH-MM-SS")
    p.add_argument("--end", default="12-02-52", help="last timestamp, HH-MM-SS")
    p.add_argument("--ncols", type=int, default=6)
    args = p.parse_args()

    files = select_files(args.start, args.end)
    if not files:
        raise SystemExit("No files in that time range.")

    ncols = args.ncols
    nrows = math.ceil(len(files) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.2 * nrows),
                             squeeze=False, layout="constrained")

    for ax, f in zip(axes.flat, files):
        df = load(f)
        ax.plot(df["time"], df["Btotal"], lw=0.8, color="tab:blue")
        ax.set_title(stamp(f), fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    for ax in axes.flat[len(files):]:
        ax.axis("off")

    fig.supxlabel("time (s)")
    fig.supylabel("Btotal (µT)")
    first, last = stamp(files[0]), stamp(files[-1])
    fig.suptitle(f"Btotal vs time, 2026-09-19, {first} to {last} "
                 f"({len(files)} files)")

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / f"btotal_{first}_to_{last}.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)} ({len(files)} files)")


if __name__ == "__main__":
    main()
