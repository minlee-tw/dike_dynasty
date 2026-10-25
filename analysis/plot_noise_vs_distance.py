"""How far from the small dike does the field settle down to noise?

Runs are lined up on the anomaly peak (cross-correlation from
plot_small_dike_aligned.py; runs without a recognisable anomaly are left out)
and each run's median is removed, so what's left is anomaly + noise. Binning
by distance from the peak, BIN_CM wide:
  - top: every aligned run, with the median and the 16-84 % band per bin;
  - bottom: the RMS of all samples in each bin, for the strong positions
    (3-9) and the weak ones (the rest) separately. Near the dike the RMS is
    dominated by the anomaly; far away it flattens to the noise level. Where
    it flattens is the distance at which the far field starts.
The far-field sigma from small_dike_noise.py and the cutoff it uses (FAR_CM)
are drawn for comparison.

Usage (from the repo root):
    python analysis/plot_noise_vs_distance.py
    python analysis/plot_noise_vs_distance.py --component Bz
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import COMPONENTS, regular_runs
from plot_small_dike_aligned import align
from small_dike_noise import FAR_CM

BIN_CM = 2.0
STRONG = range(3, 10)
EDGES = np.arange(-50, 50 + BIN_CM, BIN_CM)
CENTRES = (EDGES[:-1] + EDGES[1:]) / 2


def binned(x, y, stat):
    idx = np.digitize(x, EDGES) - 1
    return np.array([stat(y[idx == i]) if (idx == i).sum() >= 20 else np.nan
                     for i in range(len(CENTRES))])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--component", default="Btotal", choices=COMPONENTS)
    p.add_argument("--sigma", type=float, default=None,
                   help="far-field sigma to draw, µT (default: 0.155 for "
                        "Btotal, 0.244 for Bz, from small_dike_noise.py)")
    args = p.parse_args()
    comp = args.component
    sigma = args.sigma or {"Btotal": 0.155, "Bz": 0.244}.get(comp)

    results, _ = align(regular_runs("small-dike"), comp)
    groups = {"strong (positions 3-9)": ([], []), "weak (positions 1-2, 10-20)": ([], [])}
    fig, (ax_top, ax_rms) = plt.subplots(2, 1, figsize=(16, 12), sharex=True,
                                         layout="constrained")
    for r in results:
        if not r["aligned"]:
            continue
        df = load(r["path"])
        x = df["distance_cm"].to_numpy() - r["offset_cm"]
        y = df[comp].to_numpy()
        y = y - np.median(y)
        key = list(groups)[0 if r["pos"] in STRONG else 1]
        groups[key][0].append(x)
        groups[key][1].append(y)
        ax_top.plot(x, y, color="tab:blue" if r["pos"] in STRONG else "tab:orange",
                    lw=0.5, alpha=0.35)

    colors = ["tab:blue", "tab:orange"]
    for (label, (xs, ys)), c in zip(groups.items(), colors):
        x, y = np.concatenate(xs), np.concatenate(ys)
        med = binned(x, y, np.median)
        lo = binned(x, y, lambda v: np.percentile(v, 16))
        hi = binned(x, y, lambda v: np.percentile(v, 84))
        ax_top.plot(CENTRES, med, color=c, lw=2.5, label=f"{label}: median")
        ax_top.fill_between(CENTRES, lo, hi, color=c, alpha=0.2,
                            label=f"{label}: 16-84 %")
        rms = binned(x, y, lambda v: np.sqrt(np.mean(v ** 2)))
        ax_rms.plot(CENTRES, rms, "o-", color=c, lw=2,
                    label=f"{label}, {len(xs)} runs")

    for ax in (ax_top, ax_rms):
        for edge in (-FAR_CM, FAR_CM):
            ax.axvline(edge, color="tab:pink", lw=1.5, ls="--")
        ax.axvspan(-60, -FAR_CM, color="tab:pink", alpha=0.08)
        ax.axvspan(FAR_CM, 60, color="tab:pink", alpha=0.08)
        ax.grid(alpha=0.3)
    if sigma:
        ax_rms.axhline(sigma, color="k", ls=":", lw=1.5,
                       label=f"far-field noise σ = {sigma:.3f} µT")

    ax_top.set_ylabel(f"{comp} minus run median (µT)", fontsize=12)
    ax_top.set_title("Aligned runs, each run's median removed "
                     f"(lines = runs, thick = median per {BIN_CM:.0f} cm bin)",
                     fontsize=13)
    ax_top.legend(fontsize=10, loc="lower right")
    ax_rms.set_yscale("log")
    ax_rms.set_xlim(-50, 50)
    ax_rms.set_xlabel("Distance from anomaly peak (cm)", fontsize=12)
    ax_rms.set_ylabel(f"RMS of {comp} per bin (µT, log)", fontsize=12)
    ax_rms.set_title("RMS vs distance: high near the dike (signal), flattens to "
                     f"noise in the far field (dashed = ±{FAR_CM:.0f} cm cutoff)",
                     fontsize=13)
    ax_rms.legend(fontsize=11)
    fig.suptitle(f"Small dike: where does the {comp} anomaly fade into noise? "
                 "(2026-09-19)", fontsize=15)

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / f"small-dike_{comp.lower()}_noise_vs_distance.png"
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
