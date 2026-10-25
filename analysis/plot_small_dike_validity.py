"""Plot all 60 small-dike runs and sort them into valid / marginal / no dike.

Every run gets the anomaly window found by cross-correlation with the
template of clear peaks (plot_small_dike_aligned.py), ±WINDOW_CM around the
best match, even runs the alignment rejected. Two tests on each run:
  1. Shape: correlation of that window with the template (corr).
  2. Stands out from the rest of the run: peak-to-peak of the field smoothed
     over SMOOTH_CM inside the window, divided by the largest peak-to-peak in
     any other window of the same width that doesn't overlap it (ratio).
     A run with no dike still has a "biggest bump"; this asks whether that
     bump is bigger than what the same run does elsewhere, so it uses the
     run's own noise, slow wobble included, as the yardstick.
Verdict:
  valid     ratio >= MIN_RATIO and either corr >= MIN_CORR or the window
            varies by at least BIG_PTP_UT whatever its shape (a big anomaly
            of a different shape, e.g. the dip-only position 3 runs; same
            idea as BIG_RANGE_UT in the alignment)
  marginal  only one of corr >= MIN_CORR, ratio >= MIN_RATIO
  no dike   neither
A high corr on its own is weak evidence: smooth wobble can match the
template once the best of ±12 cm of shifts is picked, so the ratio is the
deciding test.
Runs too short for a comparison window get ratio = nan and are judged on
corr alone (at best marginal).

Writes figures/small-dike_<component>_validity.png (positions down, runs
across; panel colour = verdict) and figures/small-dike_<component>_validity.csv
so other scripts can keep only the valid runs.

Usage (from the repo root):
    python analysis/plot_small_dike_validity.py
    python analysis/plot_small_dike_validity.py --component Bz
"""
import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import COMPONENTS, NAMED_DIR, regular_runs
from plot_profiles_distance import moving_average
from plot_small_dike_aligned import SMOOTH_CM, WINDOW_CM, align
from small_dike_noise import notes

MIN_CORR = 0.7
MIN_RATIO = 2.0
BIG_PTP_UT = 3.0
STEP_CM = 1.0
VERDICT_COLORS = {"valid": "tab:green", "marginal": "tab:orange", "no dike": "tab:red"}


def stand_out(x, y, centre):
    """(window peak-to-peak, largest other peak-to-peak) of y smoothed over SMOOTH_CM."""
    s = moving_average(x, y, SMOOTH_CM)

    def ptp(c):
        m = (x >= c - WINDOW_CM) & (x <= c + WINDOW_CM)
        return np.ptp(s[m]) if m.sum() > 10 else np.nan

    inside = ptp(centre)
    others = [ptp(c) for c in np.arange(x.min() + WINDOW_CM, x.max() - WINDOW_CM + 1e-9, STEP_CM)
              if abs(c - centre) >= 2 * WINDOW_CM]
    others = [o for o in others if np.isfinite(o)]
    return inside, (max(others) if others else np.nan)


def verdict(corr, ratio, inside):
    shape, big = corr >= MIN_CORR, np.isfinite(ratio) and ratio >= MIN_RATIO
    if big and (shape or inside >= BIG_PTP_UT):
        return "valid"
    return "marginal" if shape or big else "no dike"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--component", default="Btotal", choices=COMPONENTS)
    args = p.parse_args()
    comp = args.component

    results, _ = align(regular_runs("small-dike"), comp)
    note = notes()
    for r in results:
        df = load(r["path"])
        r["x"], r["y"] = df["distance_cm"].to_numpy(), df[comp].to_numpy()
        r["inside"], r["other"] = stand_out(r["x"], r["y"], r["match_cm"])
        r["ratio"] = r["inside"] / r["other"]
        r["verdict"] = verdict(r["r"], r["ratio"], r["inside"])
        r["note"] = note.get(r["path"], "")

    counts = {v: sum(r["verdict"] == v for r in results) for v in VERDICT_COLORS}
    print(", ".join(f"{n} {v}" for v, n in counts.items()))
    for r in results:
        print(f"position {r['pos']:2d} run {r['run']}: corr {r['r']:.2f}, "
              f"window {r['inside']:4.1f} µT vs elsewhere {r['other']:4.1f} µT "
              f"(ratio {r['ratio']:4.1f}) -> {r['verdict']:8s} {r['note']}")

    positions = sorted({r["pos"] for r in results})
    runs = sorted({r["run"] for r in results})
    fig, axes = plt.subplots(len(positions), len(runs), figsize=(15, 2.1 * len(positions)),
                             sharex=True, squeeze=False, layout="constrained")
    for r in results:
        ax = axes[positions.index(r["pos"]), runs.index(r["run"])]
        c = VERDICT_COLORS[r["verdict"]]
        ax.plot(r["x"], r["y"], color="0.6", lw=0.5)
        ax.plot(r["x"], moving_average(r["x"], r["y"], SMOOTH_CM), color="k", lw=1)
        ax.axvspan(r["match_cm"] - WINDOW_CM, r["match_cm"] + WINDOW_CM, color=c, alpha=0.15)
        ax.axvline(r["match_cm"], color=c, lw=1)
        for spine in ax.spines.values():
            spine.set_edgecolor(c)
            spine.set_linewidth(2.5)
        flag = " *friction" if "friction" in r["note"] else ""
        ax.set_title(f"P{r['pos']:02d} run {r['run']}: {r['verdict']}{flag}\n"
                     f"corr {r['r']:.2f}, ratio {r['ratio']:.1f} "
                     f"({r['inside']:.1f} vs {r['other']:.1f} µT)",
                     fontsize=9, color=c)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel("Distance (cm)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel(f"{comp} (µT)", fontsize=9)
    fig.suptitle(f"Small dike, all {len(results)} runs: is there a dike anomaly? "
                 f"({counts['valid']} valid, {counts['marginal']} marginal, "
                 f"{counts['no dike']} no dike)\n"
                 f"Shaded = best match to the template (±{WINDOW_CM:.0f} cm); valid = "
                 f"window ≥ {MIN_RATIO:.0f}× the biggest bump elsewhere in the run and "
                 f"corr ≥ {MIN_CORR} (or any shape if ≥ {BIG_PTP_UT:.0f} µT). Grey = raw, black = {SMOOTH_CM:.0f} cm smooth",
                 fontsize=12)

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / f"small-dike_{comp.lower()}_validity.png"
    fig.savefig(out, dpi=110)
    table = FIG_DIR / f"small-dike_{comp.lower()}_validity.csv"
    with open(table, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["position", "run", "labeled_file", "match_cm", "corr",
                    "window_ptp_uT", "elsewhere_ptp_uT", "ratio", "verdict", "note"])
        for r in results:
            w.writerow([r["pos"], r["run"], r["path"].relative_to(NAMED_DIR).as_posix(),
                        f"{r['match_cm']:.1f}", f"{r['r']:.2f}", f"{r['inside']:.2f}",
                        f"{r['other']:.2f}", f"{r['ratio']:.2f}", r["verdict"], r["note"]])
    print(f"Saved {out.relative_to(REPO)} and {table.relative_to(REPO)}")


if __name__ == "__main__":
    main()
