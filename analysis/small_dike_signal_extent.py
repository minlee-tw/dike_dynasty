"""Is there a dike signal in each small-dike run, and where does it begin and end?

Works on the anomaly's size rather than its slope: weak anomalies are broad,
so their slope never stands out from noise, but the field does sit clearly
above the background level.
  1. Anomaly a(x): Btotal minus a straight line fitted to the far field
     (at least FAR_CM from the dike match), smoothed over SMOOTH_CM.
     Repeated rows are dropped first.
  2. Noise sigma: robust spread of a(x) in the far field (at least NOISE_CM
     from the match) of the runs small_dike_quality.py calls good or rotated.
  3. Null: in the same far fields, the biggest |a| a run reaches by chance
     tops out at ALONE_SIGMA (6.5 sigma, the largest of 50 runs).
  4. Each run's peak is its largest |a| within PEAK_CM of the match. The dike
     is detected if
       - |peak| >= ALONE_SIGMA, or
       - |peak| >= PAIR_SIGMA and another run at the same position has a
         peak of the same sign, also >= PAIR_SIGMA, within PAIR_CM of it.
     Two runs both reaching PAIR_SIGMA by chance in the same place with the
     same sign is very unlikely, so repeatability stands in for size.
  5. Extent: starting from the peak, grow outward while |a| stays above
     EDGE_SIGMA, bridging dips below it shorter than GAP_CM (the anomaly
     crosses zero between its dip and its peak). The edges are where the
     signal begins and ends, in cm from the match.

Prints which runs have a detected dike and the interval needed to contain
every run's signal, with and without the runs small_dike_quality.py calls
bad (and structure, whose second feature can extend the extent), and draws each run's extent (figures/small-dike_signal_extent.png).

Usage (from the repo root):
    python analysis/small_dike_quality.py        # once, for the verdicts
    python analysis/small_dike_signal_extent.py
"""
import csv

import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import regular_runs
from plot_profiles_distance import moving_average
from plot_small_dike_aligned import align
from small_dike_noise import robust_sd
from small_dike_quality import FAR_CM

SMOOTH_CM = 2.0
NOISE_CM = 30.0
PEAK_CM = 15.0
ALONE_SIGMA = 6.5
PAIR_SIGMA = 4.0
PAIR_CM = 5.0
EDGE_SIGMA = 3.0
GAP_CM = 3.0
VERDICT_COLORS = {"good": "tab:blue", "rotated": "tab:green",
                  "structure": "tab:purple", "bad": "tab:red"}


def anomaly(r):
    df = load(r["path"])
    keep = np.r_[True, np.diff(df[["Bx", "By", "Bz"]].to_numpy(), axis=0).any(axis=1)]
    x = df["distance_cm"].to_numpy()[keep] - r["match_cm"]
    y = df["Btotal"].to_numpy()[keep]
    far = np.abs(x) >= FAR_CM
    base = np.polyval(np.polyfit(x[far], y[far], 1), x)
    return x, moving_average(x, y - base, SMOOTH_CM)


def extent(x, z, i):
    """(left, right) cm around index i where |z| > EDGE_SIGMA, bridging gaps < GAP_CM."""
    above = np.abs(z) > EDGE_SIGMA
    edges = []
    for step in (-1, 1):
        j, last = i, i
        while 0 <= j + step < len(x):
            j += step
            if above[j]:
                last = j
            elif abs(x[j] - x[last]) >= GAP_CM:
                break
        edges.append(x[last])
    return tuple(edges)


def main():
    with open(FIG_DIR / "small-dike_quality.csv", newline="") as fh:
        verdict = {(int(r["position"]), int(r["run"])): r["verdict"]
                   for r in csv.DictReader(fh)}
    results, _ = align(regular_runs("small-dike"), "Btotal")
    runs = []
    for r in results:
        x, a = anomaly(r)
        runs.append(dict(pos=r["pos"], run=r["run"], v=verdict[r["pos"], r["run"]], x=x, a=a))

    sigma = robust_sd(np.concatenate([r["a"][np.abs(r["x"]) >= NOISE_CM]
                                      for r in runs if r["v"] in ("good", "rotated")]))
    for r in runs:
        r["z"] = r["a"] / sigma
        near = np.flatnonzero(np.abs(r["x"]) <= PEAK_CM)
        r["i"] = near[np.argmax(np.abs(r["z"][near]))]
        r["peak"], r["x_peak"] = r["z"][r["i"]], r["x"][r["i"]]

    # Bad runs can't vouch for others: a run's repeat partners are the other
    # non-bad runs at its position (any run, for a bad one).
    for r in runs:
        partners = [o for o in runs if o["pos"] == r["pos"] and o is not r
                    and (o["v"] != "bad" or r["v"] == "bad")
                    and abs(o["peak"]) >= PAIR_SIGMA
                    and np.sign(o["peak"]) == np.sign(r["peak"])
                    and abs(o["x_peak"] - r["x_peak"]) <= PAIR_CM]
        r["how"] = ("alone" if abs(r["peak"]) >= ALONE_SIGMA else
                    "repeat" if abs(r["peak"]) >= PAIR_SIGMA and partners else None)
        r["extent"] = extent(r["x"], r["z"], r["i"]) if r["how"] else None

    for name, keep in (("all runs", lambda r: True),
                       ("without bad runs", lambda r: r["v"] != "bad"),
                       ("without bad or structure runs", lambda r: r["v"] not in ("bad", "structure"))):
        sel = [r for r in runs if keep(r)]
        found = [r for r in sel if r["extent"]]
        left = np.array([r["extent"][0] for r in found])
        right = np.array([r["extent"][1] for r in found])
        print(f"\n{name}: dike detected in {len(found)} of {len(sel)} runs "
              f"({sum(r['how'] == 'alone' for r in found)} alone, "
              f"{sum(r['how'] == 'repeat' for r in found)} by repeat)")
        print("  not detected: " + ", ".join(
            f"P{r['pos']:02d} r{r['run']} ({r['peak']:+.1f}σ)" for r in sel if not r["extent"]))
        print(f"  signal begins: median {np.median(left):+.1f} cm, "
              f"furthest {left.min():+.1f} cm (P{found[left.argmin()]['pos']:02d} "
              f"r{found[left.argmin()]['run']})")
        print(f"  signal ends:   median {np.median(right):+.1f} cm, "
              f"furthest {right.max():+.1f} cm (P{found[right.argmax()]['pos']:02d} "
              f"r{found[right.argmax()]['run']})")
        print(f"  interval containing every run's signal: [{left.min():+.1f}, {right.max():+.1f}] cm; "
              f"90 %: [{np.percentile(left, 10):+.1f}, {np.percentile(right, 90):+.1f}] cm")
    print(f"\nnoise sigma of the {SMOOTH_CM:.0f} cm-smoothed anomaly: {sigma:.3f} µT")

    fig, ax = plt.subplots(figsize=(13, 14), layout="constrained")
    for k, r in enumerate(runs):
        y = len(runs) - k
        c = VERDICT_COLORS[r["v"]]
        ax.plot([r["x"].min(), r["x"].max()], [y, y], color="0.85", lw=0.8, zorder=0)
        if r["extent"]:
            ax.plot(r["extent"], [y, y], color=c, lw=5, solid_capstyle="butt")
            ax.plot(r["x_peak"], y, "k|", ms=9)
        else:
            ax.plot(r["x_peak"], y, "x", color=c, ms=6)
        ax.annotate(f"{r['peak']:+.1f}σ" + (" (repeat)" if r["how"] == "repeat" else ""),
                    (52, y), fontsize=7, va="center", color=c)
    ax.set_yticks(range(len(runs), 0, -1),
                  [f"P{r['pos']:02d} r{r['run']}" for r in runs], fontsize=7)
    ax.axvline(0, color="k", lw=0.8)
    for v, c in VERDICT_COLORS.items():
        ax.plot([], [], color=c, lw=5, label=v)
    ax.plot([], [], "kx", label="peak, not detected")
    ax.legend(title="Quality verdict", loc="upper left", fontsize=9)
    ax.set_xlim(-60, 66)
    ax.set_xlabel("Distance from best match to the dike template (cm)", fontsize=12)
    ax.set_title(f"Small dike: where the dike signal begins and ends in each run\n"
                 f"detected if peak ≥ {ALONE_SIGMA}σ, or ≥ {PAIR_SIGMA:.0f}σ and repeated "
                 f"at the same position; bar = where |anomaly| > {EDGE_SIGMA:.0f}σ "
                 f"(σ = {sigma:.3f} µT after {SMOOTH_CM:.0f} cm smoothing)", fontsize=12)
    ax.grid(alpha=0.3, axis="x")
    out = FIG_DIR / "small-dike_signal_extent.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
