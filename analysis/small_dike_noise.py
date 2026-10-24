"""Estimate the far-field noise of the small-dike runs and how much of the
dike signature stands above it.

Following the lecture 8 strategy:
  1. Isolate the parts of each transect that should reflect only noise: the
     far field, more than FAR_CM from the anomaly peak on either side (runs
     are lined up with the cross-correlation from plot_small_dike_aligned.py;
     runs with no recognisable anomaly are skipped).
  2. For each far-field sample, compute a deviation attributable to random
     noise: the sample minus a straight line fitted to that side of that run.
     The line soaks up each run's baseline level and any slow drift (phone
     orientation, regional gradient), so what's left is the scatter.
  3. Plot the distribution of these deviations; its width sigma (robust:
     1.4826 x median absolute deviation) is the random noise per sample.
Then, for the near field (within FAR_CM of the peak), the anomaly is the field
minus the baseline interpolated between the two far-field lines, and:
  - amplitude = peak-to-peak of the anomaly smoothed over SMOOTH_CM,
  - SNR = amplitude / sigma,
  - signal fraction = 1 - sigma^2 / mean(anomaly^2), the share of the
    near-field variation that is not random noise.
The spread of the far-field levels between runs at the same position is a
separate, systematic error (bias): orientation, height above the rock, other
phones nearby.

Writes figures/small-dike_<component>_noise.png and a per-run table
figures/small-dike_<component>_noise.csv.

Usage (from the repo root):
    python analysis/small_dike_noise.py
    python analysis/small_dike_noise.py --component Bz
"""
import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import COMPONENTS, NAMED_DIR, regular_runs
from plot_profiles_distance import moving_average
from plot_small_dike_aligned import SMOOTH_CM, align

FAR_CM = 25.0
MIN_SIDE_CM = 3.0


def notes():
    """{labeled path: field note} for the small-dike runs."""
    with open(NAMED_DIR / "metadata.csv", newline="") as fh:
        return {NAMED_DIR / r["labeled_file"]: r["note"] for r in csv.DictReader(fh)}


def robust_sd(d):
    return 1.4826 * np.median(np.abs(d - np.median(d)))


def analyse(r, component):
    """Far-field deviations and near-field signal measures for one aligned run."""
    df = load(r["path"])
    x = df["distance_cm"].to_numpy() - r["offset_cm"]
    y = df[component].to_numpy()
    fits, devs = {}, []
    for side, mask in (("left", x <= -FAR_CM), ("right", x >= FAR_CM)):
        if mask.sum() < 2 or np.ptp(x[mask]) < MIN_SIDE_CM:
            continue
        coef = np.polyfit(x[mask], y[mask], 1)
        fits[side] = (coef, x[mask])
        devs.append(y[mask] - np.polyval(coef, x[mask]))
    if not fits:
        return None

    near = np.abs(x) < FAR_CM
    if len(fits) == 2:  # baseline: straight line between the two far-field fits
        (cl, _), (cr, _) = fits["left"], fits["right"]
        yl, yr = np.polyval(cl, -FAR_CM), np.polyval(cr, FAR_CM)
        base = yl + (yr - yl) * (x[near] + FAR_CM) / (2 * FAR_CM)
    else:  # one side only: hold its level
        coef, xs = next(iter(fits.values()))
        base = np.full(near.sum(), np.polyval(coef, xs).mean())
    anom = y[near] - base
    levels = [np.polyval(c, xs).mean() for c, xs in fits.values()]
    return dict(r, x=x, y=y, fits=fits, dev=np.concatenate(devs),
                near_x=x[near], anom=anom,
                amplitude=np.ptp(moving_average(x[near], anom, SMOOTH_CM)),
                ms=np.mean(anom ** 2), level=np.mean(levels),
                far_cm=sum(np.ptp(xs) for _, xs in fits.values()))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--component", default="Btotal", choices=COMPONENTS)
    args = p.parse_args()
    comp = args.component

    results, _ = align(regular_runs("small-dike"), comp)
    note = notes()
    runs = [a for a in (analyse(r, comp) for r in results if r["aligned"]) if a]
    skipped = len(results) - len(runs)

    dev = np.concatenate([a["dev"] for a in runs])
    # Robust sigma (1.4826 x median absolute deviation): the tails come from a
    # few runs with steps or leftover structure in the far field (position 1
    # run 1, positions 4-5) and would otherwise inflate the noise.
    sigma = robust_sd(dev)
    print(f"plain standard deviation {dev.std(ddof=1):.3f} µT (tails included)")
    for a in runs:
        a["sigma_run"] = robust_sd(a["dev"])
        a["snr"] = a["amplitude"] / sigma
        a["signal_frac"] = max(0.0, 1 - sigma ** 2 / a["ms"])
        a["note"] = note.get(a["path"], "")

    # Systematic: spread of far-field levels between runs at one position.
    by_pos = {}
    for a in runs:
        by_pos.setdefault(a["pos"], []).append(a["level"])
    level_sd = [np.std(v, ddof=1) for v in by_pos.values() if len(v) > 1]

    print(f"{len(runs)} runs used ({skipped} skipped: no anomaly or no far field)")
    print(f"random noise sigma = {sigma:.3f} µT from {len(dev)} far-field samples")
    print(f"far-field level spread between runs at one position: "
          f"median {np.median(level_sd):.2f} µT, max {np.max(level_sd):.2f} µT")
    print(f"median SNR {np.median([a['snr'] for a in runs]):.1f}, "
          f"median signal fraction {np.median([a['signal_frac'] for a in runs]):.0%}")
    for a in runs:
        print(f"position {a['pos']:2d} run {a['run']}: sigma {a['sigma_run']:.3f} µT, "
              f"amplitude {a['amplitude']:4.1f} µT, SNR {a['snr']:5.1f}, "
              f"signal {a['signal_frac']:4.0%}  {a['note']}")

    fig = plt.figure(figsize=(20, 13), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1])
    ax_runs, ax_hist, ax_snr = (fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]),
                                fig.add_subplot(gs[1, 1]))
    positions = sorted(by_pos)
    colors = dict(zip(positions, plt.cm.turbo(np.linspace(0, 1, len(positions)))))

    # 1. Aligned runs, far field shaded, with the fitted baselines.
    for a in runs:
        ax_runs.plot(a["x"], a["y"], color=colors[a["pos"]], lw=0.7, alpha=0.7)
        for coef, xs in a["fits"].values():
            ax_runs.plot(xs, np.polyval(coef, xs), color="k", lw=0.8)
    for lo, hi in ((-60, -FAR_CM), (FAR_CM, 60)):
        ax_runs.axvspan(lo, hi, color="tab:pink", alpha=0.12)
    ax_runs.set_xlim(-55, 55)
    ax_runs.set_xlabel("Distance from anomaly peak (cm)", fontsize=12)
    ax_runs.set_ylabel(f"{comp} (µT)", fontsize=12)
    ax_runs.set_title(f"Aligned runs; shaded = far field (|x| ≥ {FAR_CM:.0f} cm) "
                      "used for noise, black = straight-line baseline per side",
                      fontsize=13)
    ax_runs.grid(alpha=0.3)

    # 2. Distribution of far-field deviations.
    lim = 4 * sigma
    ax_hist.hist(dev, bins=np.linspace(-lim, lim, 61), density=True,
                 color="tab:blue", alpha=0.7, label="far-field deviations")
    g = np.linspace(-lim, lim, 400)
    ax_hist.plot(g, np.exp(-g ** 2 / (2 * sigma ** 2)) / (sigma * np.sqrt(2 * np.pi)),
                 color="k", lw=1.5, label=f"Gaussian, σ = {sigma:.3f} µT")
    ax_hist.set_xlabel(f"{comp} minus far-field baseline (µT)", fontsize=12)
    ax_hist.set_ylabel("Probability density", fontsize=12)
    ax_hist.set_title(f"Random noise: {len(dev)} samples from {len(runs)} runs",
                      fontsize=13)
    ax_hist.legend(fontsize=11)
    ax_hist.grid(alpha=0.3)

    # 3. Per-run anomaly amplitude against the noise.
    xs = np.arange(len(runs))
    flagged = ["friction" in a["note"] for a in runs]
    ax_snr.bar(xs, [a["snr"] for a in runs],
               color=[colors[a["pos"]] for a in runs],
               edgecolor=["k" if f else "none" for f in flagged], linewidth=1.5)
    ax_snr.axhline(3, color="k", ls="--", lw=1, label="SNR = 3")
    ax_snr.set_xticks(xs, [f"{a['pos']}.{a['run']}" for a in runs], rotation=90,
                      fontsize=8)
    ax_snr.set_xlabel("Position.run (black outline = friction on rock in field notes)",
                      fontsize=12)
    ax_snr.set_ylabel("Anomaly peak-to-peak / σ", fontsize=12)
    ax_snr.set_title(f"Signal-to-noise per run (median "
                     f"{np.median([a['snr'] for a in runs]):.0f}; median signal "
                     f"fraction {np.median([a['signal_frac'] for a in runs]):.0%})",
                     fontsize=13)
    ax_snr.legend(fontsize=11)
    ax_snr.grid(alpha=0.3, axis="y")

    fig.suptitle(f"Small dike: far-field noise in {comp} (2026-09-19). Run-to-run "
                 f"far-field level spread at one position: median "
                 f"{np.median(level_sd):.2f} µT (systematic, not in σ)", fontsize=15)

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / f"small-dike_{comp.lower()}_noise.png"
    fig.savefig(out, dpi=150)
    table = out.with_suffix(".csv")
    with open(table, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["position", "run", "peak_offset_cm", "far_field_cm",
                    "far_field_level_uT", "sigma_run_uT", "amplitude_uT", "snr",
                    "signal_fraction", "note"])
        for a in runs:
            w.writerow([a["pos"], a["run"], f"{a['offset_cm']:.1f}", f"{a['far_cm']:.1f}",
                        f"{a['level']:.2f}", f"{a['sigma_run']:.3f}",
                        f"{a['amplitude']:.2f}", f"{a['snr']:.1f}",
                        f"{a['signal_frac']:.3f}", a["note"]])
    print(f"Saved {out.relative_to(REPO)} and {table.relative_to(REPO)}")


if __name__ == "__main__":
    main()
