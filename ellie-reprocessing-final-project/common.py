"""Shared paths and helpers for the Bz-only reprocessing (no Btotal anywhere).

The original run files in named_data/ are only read, never changed. Everything
this folder produces goes into ellie-reprocessing-final-project/data/ (and
figures/ when a script is run with --figures).
"""
import csv
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
NAMED_DIR = REPO / "named_data"          # original data, read only
DATA_DIR = HERE / "data"
FIG_DIR = HERE / "figures"
GRID_CM = np.linspace(0, 90, 901)        # common 0.1 cm grid for template matching


def load(path):
    """One run file: '# sensor:magnetometer' line, then time, distance_cm, Bx, By, Bz, ...
    Only time, distance_cm and Bz are used here."""
    return pd.read_csv(path, comment="#")


def regular_runs(dike):
    """[(position, run, path)] for a dike's positions 1-20, runs 1-3, sorted."""
    with open(NAMED_DIR / "metadata.csv", newline="") as fh:
        rows = [r for r in csv.DictReader(fh)
                if r["dike"] == dike and r["experiment"] == "positions"]
    return sorted((int(r["position"]), int(r["run"]), NAMED_DIR / r["labeled_file"]) for r in rows)


def repeated_rows(df):
    """True where the logger wrote the same Bx, By, Bz as the previous row (a repeated
    sample, not a new measurement). These rows are left out of every calculation."""
    v = df[["Bx", "By", "Bz"]].to_numpy()
    return np.r_[False, ~np.diff(v, axis=0).any(axis=1)]


def moving_average(x, y, width):
    """Centred sliding mean of y over a window `width` wide in x (x sorted). Used only to
    prepare run shapes for template matching, never for the dike cutoffs."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    csum = np.concatenate([[0.0], np.cumsum(y)])
    lo = np.searchsorted(x, x - width / 2, side="left")
    hi = np.searchsorted(x, x + width / 2, side="right")
    return (csum[hi] - csum[lo]) / (hi - lo)


def robust_sd(d):
    """1.4826 x median absolute deviation: an SD that ignores outliers."""
    return 1.4826 * np.median(np.abs(d - np.median(d)))


def loess_slope(x, y, span):
    """Local-linear LOESS slope of y vs x: at each sample, a straight line fitted by
    weighted least squares to the samples within ±span/2, with tricube weights
    (1 - |d / (span/2)|^3)^3. nan where the span runs past either end."""
    h = span / 2
    lo = np.searchsorted(x, x - h)
    hi = np.searchsorted(x, x + h, side="right")
    out = np.full(len(x), np.nan)
    for i in range(len(x)):
        if x[i] - h < x[0] or x[i] + h > x[-1]:
            continue
        d = x[lo[i]:hi[i]] - x[i]
        v = y[lo[i]:hi[i]]
        w = (1 - np.abs(d / h) ** 3) ** 3
        mx = (w * d).sum() / w.sum()
        my = (w * v).sum() / w.sum()
        den = (w * (d - mx) ** 2).sum()
        if den > 0:
            out[i] = (w * (d - mx) * (v - my)).sum() / den
    return out


def run_name(dike, pos, run, path):
    """Output file name for a run, keeping the original time stamp: e.g.
    small-dike-position-01-run-1_11-09-51_bz.csv."""
    return f"{Path(path).stem}_bz.csv"
