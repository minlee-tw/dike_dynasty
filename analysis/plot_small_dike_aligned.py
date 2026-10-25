"""Overlay all small-dike runs (positions 1-20, runs 1-3) with the anomaly lined up.

Each run crosses the dike at a slightly different distance (start point and
hand speed vary), so the anomaly lands at different distance_cm values. This
shifts every run sideways so the anomalies coincide, and plots field vs
distance from the anomaly peak (0 cm = peak).

Alignment, on a common 0.1 cm grid with each run's median removed and a
SMOOTH_CM moving average applied:
  1. Template: runs with a clear peak (more than TEMPLATE_PEAK_UT above their
     median) are shifted so their maxima line up and averaged.
  2. Every run is cross-correlated with the template, using only the
     template's anomaly (±WINDOW_CM around its peak) and shifts of up to
     ±MAX_SHIFT_CM. The run is shifted to match if the best correlation is at
     least MIN_CORR and the matched stretch varies by at least MIN_RANGE_UT,
     or if the anomaly is large (at least BIG_RANGE_UT) and the correlation is
     at least BIG_MIN_CORR, which lets dip-only runs like position 3 through;
     otherwise there is no recognisable anomaly to align on, so the run
     is left where it is and drawn faint.
Cross-correlation uses the whole low-high shape, so runs with only a dip
(e.g. position 3) still line up, and the shift doesn't depend on one noisy
maximum.

Usage (from the repo root):
    python analysis/plot_small_dike_aligned.py
    python analysis/plot_small_dike_aligned.py --component Bz
    python analysis/plot_small_dike_aligned.py --dots
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import COMPONENTS, GRID_CM, RUN_STYLES, regular_runs
from plot_profiles_distance import moving_average

SMOOTH_CM = 2.0
TEMPLATE_PEAK_UT = 2.0
WINDOW_CM = 15.0
MAX_SHIFT_CM = 12.0
MIN_CORR = 0.7
MIN_RANGE_UT = 0.8
BIG_RANGE_UT = 3.0
BIG_MIN_CORR = 0.5
STEP_CM = GRID_CM[1] - GRID_CM[0]


def anomaly(df, component):
    y = np.interp(GRID_CM, df["distance_cm"], df[component])
    return moving_average(GRID_CM, y - np.median(y), SMOOTH_CM)


def shift(a, n):
    """Shift a by n samples (positive = later), padding with zeros."""
    out = np.zeros_like(a)
    if n >= 0:
        out[n:] = a[:len(a) - n]
    else:
        out[:n] = a[-n:]
    return out


def best_lag(a, template, peak_idx):
    """(lag, correlation, range) best matching a to the template's anomaly.

    lag is in samples (run's anomaly sits lag samples after the template's);
    range is the peak-to-peak size of the matched stretch of a, in µT.
    """
    half = int(round(WINDOW_CM / STEP_CM))
    max_lag = int(round(MAX_SHIFT_CM / STEP_CM))
    lo, hi = peak_idx - half, peak_idx + half
    t = template[lo:hi]
    best = (0, -np.inf, 0.0)
    for lag in range(-max_lag, max_lag + 1):
        if lo + lag < 0 or hi + lag > len(a):
            continue
        seg = a[lo + lag:hi + lag]
        r = np.corrcoef(seg, t)[0, 1]
        if r > best[1]:
            best = (lag, r, seg.max() - seg.min())
    return best


def align(runs, component):
    anomalies = [anomaly(load(path), component) for _, _, path in runs]

    # 1. Template from the clear peaks, lined up on their maxima.
    strong = [a for a in anomalies if a.max() > TEMPLATE_PEAK_UT]
    peaks = [int(np.argmax(a)) for a in strong]
    ref = int(np.median(peaks))
    template = np.mean([shift(a, ref - p) for a, p in zip(strong, peaks)], axis=0)
    peak_idx = int(np.argmax(template))
    peak_cm = GRID_CM[peak_idx]

    # 2. Cross-correlate every run with the template's anomaly.
    results = []
    for (pos, run, path), a in zip(runs, anomalies):
        lag, r, size = best_lag(a, template, peak_idx)
        aligned = ((r >= MIN_CORR and size >= MIN_RANGE_UT)
                   or (r >= BIG_MIN_CORR and size >= BIG_RANGE_UT))
        results.append(dict(pos=pos, run=run, path=path, r=r, size=size,
                            aligned=aligned, match_cm=peak_cm + lag * STEP_CM,
                            offset_cm=peak_cm + (lag * STEP_CM if aligned else 0)))
    return results, len(strong)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--component", default="Btotal", choices=COMPONENTS)
    p.add_argument("--dots", action="store_true",
                   help="draw each sample as a dot instead of a line")
    args = p.parse_args()

    runs = regular_runs("small-dike")
    results, n_template = align(runs, args.component)
    positions = sorted({r["pos"] for r in results})
    colors = dict(zip(positions, plt.cm.turbo(np.linspace(0, 1, len(positions)))))

    fig, ax = plt.subplots(figsize=(20, 11), layout="constrained")
    for r in results:
        df = load(r["path"])
        x = df["distance_cm"] - r["offset_cm"]
        label = f"Position {r['pos']:02d}, run {r['run']}"
        if not r["aligned"]:
            label += " (not aligned)"
        fade = dict(alpha=1.0) if r["aligned"] else dict(alpha=0.25)
        style = (dict(ls="none", marker="o", ms=1.6, mew=0) if args.dots else
                 dict(ls=RUN_STYLES.get(r["run"], "-."), lw=0.9 if r["aligned"] else 0.6))
        ax.plot(x, df[args.component], color=colors[r["pos"]], label=label,
                **style, **fade)
        print(f"position {r['pos']:2d} run {r['run']}: corr {r['r']:.2f}, "
              f"size {r['size']:4.1f} µT, "
              + (f"peak at {r['offset_cm']:5.1f} cm" if r["aligned"] else "not aligned"))

    n_aligned = sum(r["aligned"] for r in results)
    ax.axvline(0, color="k", lw=1, alpha=0.5)
    ax.set_xlim(-50, 50)
    ax.set_xlabel("Distance from anomaly peak (cm)", fontsize=13)
    ax.set_ylabel(f"{args.component} (µT)", fontsize=13)
    ax.set_title(f"Small dike: {args.component} vs distance, positions 1-20, runs 1-3, "
                 f"anomalies lined up ({n_aligned} of {len(results)} runs aligned, "
                 f"2026-09-19)\nShift from cross-correlation with a template of the "
                 f"{n_template} clearest peaks. Faint = no anomaly found, left unshifted",
                 fontsize=14)
    ax.tick_params(labelsize=11)
    ax.grid(alpha=0.3)
    leg = ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), ncols=3,
                    fontsize=10, frameon=False, handlelength=3,
                    title="Position, run (colour = position, line = run)",
                    title_fontsize=12)
    for line in leg.get_lines():
        line.set_linewidth(3)
        line.set_markersize(8)

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / (f"small-dike_{args.component.lower()}_distance_aligned"
                     f"{'_dots' if args.dots else ''}.png")
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
