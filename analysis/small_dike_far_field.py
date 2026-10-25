"""Far-field cutoffs for each small-dike run, from the slope of the field.

The far field is where the field isn't changing faster than noise. For each
run (repeated rows dropped):
  1. Slope = least-squares slope of Btotal vs distance over a centred window
     WINDOW_CM wide.
  2. Slope noise = robust spread of that slope at least NOISE_CM from the
     dike match, pooled over the runs small_dike_quality.py calls good or
     rotated.
  3. Signal = |slope| > K_SIGMA x noise for a stretch at least PERSIST_CM
     long. The dike's signal is the stretch nearest the template match
     (within SEARCH_CM), joined with neighbours less than GAP_CM away (the
     slope passes through zero at the dip and the peak).
  4. Cutoffs: the far field is everything more than WINDOW_CM / 2 beyond the
     dike's signal (a slope at x uses data out to x ± WINDOW_CM / 2). Each run
     gets its own left and right cutoff. A run with no dike signal is far
     field throughout.

All 60 runs are included and plotted alike, the bad ones too;
small_dike_quality.py's verdicts are only used to choose which far fields set
the slope noise (and, in --scan, to leave out known far-field structure).

--scan tries a grid of windows and thresholds and reports how many runs each
one finds the dike in, and how many runs have a false alarm (a signal stretch
at least NOISE_CM from the match, outside the positions with known far-field
structure).

Also prints, and adds to the CSV, each run's far-field mean Btotal on each
side of the dike (with its SD and length); a run with no dike
is split at the template match. The far-field value of each run is the mean of
all its far-field samples, both sides together, with the SD of those samples
and the standard error of the mean from BLOCK_CM block means (neighbouring
samples are correlated, so SD / sqrt(n) would be too small).

Writes figures/small-dike_far_field_cutoffs.csv (cutoffs in cm along the run,
as in distance_cm, and far-field means), figures/small-dike_far_field_cutoffs.png (cutoffs as bars, one row
per run), figures/small-dike_far_field_cutoffs_overlay.png (every run on one
axis aligned on the dike match, cutoffs marked), the same at the recorded
distances (..._overlay_unaligned.png), figures/small-dike_far_field_means.png
(far-field means per side), figures/small-dike_far_field_value.png (one far-field
value per run), figures/small-dike_far_field_means_grid.png (one
panel per run: samples, dike shaded, far-field mean ± SD per side) and
figures/small-dike_far_field_cutoffs_grid.png (one panel per run).

Usage (from the repo root):
    python analysis/small_dike_quality.py        # once, for the verdicts
    python analysis/small_dike_far_field.py --scan
    python analysis/small_dike_far_field.py
"""
import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import regular_runs
from plot_profiles_distance import moving_average
from plot_small_dike_aligned import align
from small_dike_noise import robust_sd

# Fixed parameters, agreed with the group (2026-10-24): use these for all
# far-field work; --scan is only for showing why they were chosen.
WINDOW_CM = 4.0
K_SIGMA = 2.5
PERSIST_CM = 0.5
NOISE_CM = 30.0
SEARCH_CM = 10.0
GAP_CM = 5.0
BLOCK_CM = 4.0  # far-field samples are averaged in blocks this long for standard errors
DIKE_COLOR = "tab:orange"


def slope(x, y, width):
    """Least-squares slope of y vs x over a centred window `width` wide (nan near the ends)."""
    lo = np.searchsorted(x, x - width / 2)
    hi = np.searchsorted(x, x + width / 2, side="right")

    def window_sum(v):
        c = np.concatenate([[0.0], np.cumsum(v)])
        return c[hi] - c[lo]

    n = hi - lo
    sx, sy, sxx, sxy = window_sum(x), window_sum(y), window_sum(x * x), window_sum(x * y)
    den = n * sxx - sx * sx
    out = np.full(len(x), np.nan)
    ok = (den > 0) & (x - width / 2 >= x[0]) & (x + width / 2 <= x[-1])
    out[ok] = (n * sxy - sx * sy)[ok] / den[ok]
    return out


def stretches(mask, x, min_len):
    """[(x_start, x_end)] of contiguous True runs in mask at least min_len cm long."""
    out, i = [], 0
    while i < len(mask):
        if mask[i]:
            j = i
            while j < len(mask) and mask[j]:
                j += 1
            if x[j - 1] - x[i] >= min_len:
                out.append((x[i], x[j - 1]))
            i = j
        else:
            i += 1
    return out


def dike_extent(segs):
    """(left, right) of the stretch nearest 0 within SEARCH_CM, merged across gaps < GAP_CM."""
    near = [s for s in segs if s[1] >= -SEARCH_CM and s[0] <= SEARCH_CM]
    if not near:
        return None
    left, right = min(near, key=lambda s: 0 if s[0] <= 0 <= s[1] else min(abs(s[0]), abs(s[1])))
    for a, b in sorted(segs):
        if a > right and a - right < GAP_CM:
            right = b
    for a, b in sorted(segs, reverse=True):
        if b < left and left - b < GAP_CM:
            left = a
    return left, right


def load_runs():
    with open(FIG_DIR / "small-dike_quality.csv", newline="") as fh:
        verdict = {(int(r["position"]), int(r["run"])): r["verdict"]
                   for r in csv.DictReader(fh)}
    results, _ = align(regular_runs("small-dike"), "Btotal")
    runs = []
    for r in results:
        df = load(r["path"])
        keep = np.r_[True, np.diff(df[["Bx", "By", "Bz"]].to_numpy(), axis=0).any(axis=1)]
        runs.append(dict(pos=r["pos"], run=r["run"], v=verdict[r["pos"], r["run"]],
                         match_cm=r["match_cm"],
                         x=df["distance_cm"].to_numpy()[keep] - r["match_cm"],
                         y=df["Btotal"].to_numpy()[keep]))
    return runs


def detect(runs, width, k):
    """Fill each run's signal stretches and dike extent; return the slope noise."""
    for r in runs:
        r["s"] = slope(r["x"], r["y"], width)
    noise = robust_sd(np.concatenate(
        [r["s"][(np.abs(r["x"]) >= NOISE_CM) & np.isfinite(r["s"])]
         for r in runs if r["v"] in ("good", "rotated")]))
    for r in runs:
        r["segs"] = stretches(np.nan_to_num(np.abs(r["s"])) > k * noise, r["x"], PERSIST_CM)
        r["extent"] = dike_extent(r["segs"])
    return noise


def scan(runs):
    print(f"{len(runs)} runs; false alarms counted outside the positions with far-field structure")
    print("window  k    threshold   dike found     false alarms (of them in bad runs)")
    for width in (2, 3, 4, 5, 6):
        for k in (2.5, 3, 3.5):
            noise = detect(runs, width, k)
            found = sum(r["extent"] is not None for r in runs)
            alarm = [r for r in runs if r["v"] != "structure"
                     and any(min(abs(a), abs(b)) >= NOISE_CM for a, b in r["segs"])]
            print(f"{width:4.1f} cm {k:3.1f}  {k * noise:.2f} µT/cm  {found:3d}/{len(runs)} "
                  f"({found / len(runs):4.0%})   {len(alarm)} ({sum(r['v'] == 'bad' for r in alarm)})")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--scan", action="store_true", help="try a grid of windows and thresholds")
    args = p.parse_args()
    runs = load_runs()
    if args.scan:
        scan(runs)
        return

    noise = detect(runs, WINDOW_CM, K_SIGMA)
    half = WINDOW_CM / 2
    for r in runs:
        r["cut"] = (r["extent"][0] - half, r["extent"][1] + half) if r["extent"] else None
    found = [r for r in runs if r["cut"]]
    print(f"window {WINDOW_CM} cm, {K_SIGMA}σ: slope noise {noise:.3f} µT/cm, "
          f"threshold {K_SIGMA * noise:.2f} µT/cm")
    print(f"dike found in {len(found)} of {len(runs)} runs; "
          "no signal (far field throughout): "
          + ", ".join(f"P{r['pos']:02d} r{r['run']}" for r in runs if not r["cut"]))
    left = np.array([r["cut"][0] for r in found])
    right = np.array([r["cut"][1] for r in found])
    print(f"far field starts, from the dike match: left median {np.median(left):+.1f} cm "
          f"(range {left.min():+.1f} to {left.max():+.1f}), right median {np.median(right):+.1f} cm "
          f"(range {right.min():+.1f} to {right.max():+.1f})")

    # Far-field average on each side. With no dike found, the run is split at
    # the template match.
    print("\nfar-field mean Btotal (µT) on each side; ± = SD of the samples, [cm of far field]")
    for r in runs:
        a, b = r["cut"] if r["cut"] else (0.0, 0.0)
        r["sides"] = {}
        for side, m in (("left", r["x"] < a), ("right", r["x"] > b)):
            y = r["y"][m]
            r["sides"][side] = (y.mean(), y.std(ddof=1), np.ptp(r["x"][m])) if len(y) > 1 else None
        text = []
        for side in ("left", "right"):
            s = r["sides"][side]
            text.append(f"{side} {s[0]:6.2f} ± {s[1]:.2f} [{s[2]:4.1f}]" if s else f"{side}  none")
        diff = (r["sides"]["left"][0] - r["sides"]["right"][0]
                if r["sides"]["left"] and r["sides"]["right"] else np.nan)
        print(f"  P{r['pos']:02d} r{r['run']}: {text[0]}   {text[1]}   left - right {diff:+.2f}"
              + ("" if r["cut"] else "   (no dike: split at match)"))

    # One far-field value per run: the mean of every far-field sample, both
    # sides together (so the longer side counts for more). SD is the spread of
    # those samples; SE is the uncertainty of the mean, from BLOCK_CM block
    # means, because neighbouring samples are correlated.
    print("\nfar-field value per run, all far-field samples together (µT): mean ± SE; SD of samples")
    for r in runs:
        a, b = r["cut"] if r["cut"] else (0.0, 0.0)
        far = (r["x"] < a) | (r["x"] > b)
        x, y = r["x"][far], r["y"][far]
        # blocks never straddle the dike: number them separately on each side
        idx = np.where(x < a, np.floor((x - x.min()) / BLOCK_CM),
                       1000 + np.floor((x - b) / BLOCK_CM)).astype(int)
        means = np.array([y[idx == i].mean() for i in np.unique(idx) if (idx == i).sum() >= 3])
        r["far"] = (y.mean(), y.std(ddof=1), means.std(ddof=1) / np.sqrt(len(means)),
                    np.ptp(x[x < a]) + np.ptp(x[x > b]) if (x < a).any() and (x > b).any()
                    else np.ptp(x))
        print(f"  P{r['pos']:02d} r{r['run']}: {r['far'][0]:6.2f} ± {r['far'][2]:.3f} (SE)   "
              f"SD {r['far'][1]:.2f}   [{r['far'][3]:.0f} cm, {far.sum()} samples]")

    out = FIG_DIR / "small-dike_far_field_cutoffs.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["position", "run", "verdict", "match_cm", "dike_found", "far_left_below_cm",
                    "far_right_above_cm", "left_from_match_cm", "right_from_match_cm",
                    "far_left_mean_uT", "far_left_sd_uT", "far_left_length_cm",
                    "far_right_mean_uT", "far_right_sd_uT", "far_right_length_cm",
                    "far_mean_uT", "far_sd_uT", "far_se_uT", "far_length_cm"])
        for r in runs:
            if r["cut"]:
                a, b = r["cut"]
                cut = [f"{a + r['match_cm']:.2f}", f"{b + r['match_cm']:.2f}", f"{a:.2f}", f"{b:.2f}"]
            else:
                cut = ["", "", "", ""]
            sides = []
            for side in ("left", "right"):
                s = r["sides"][side]
                sides += [f"{s[0]:.3f}", f"{s[1]:.3f}", f"{s[2]:.1f}"] if s else ["", "", ""]
            far = [f"{v:.3f}" for v in r["far"][:3]] + [f"{r['far'][3]:.1f}"]
            w.writerow([r["pos"], r["run"], r["v"], f"{r['match_cm']:.2f}",
                        "yes" if r["cut"] else "no"] + cut + sides + far)
    print(f"Saved {out.relative_to(REPO)}")

    fig, ax = plt.subplots(figsize=(13, 14), layout="constrained")
    for i, r in enumerate(runs):
        y = len(runs) - i
        c = DIKE_COLOR
        x0, x1 = r["x"].min(), r["x"].max()
        if r["cut"]:
            ax.plot([x0, r["cut"][0]], [y, y], color="0.5", lw=1.5)
            ax.plot([r["cut"][1], x1], [y, y], color="0.5", lw=1.5)
            ax.plot(r["cut"], [y, y], color=c, lw=6, solid_capstyle="butt")
        else:
            ax.plot([x0, x1], [y, y], color="0.5", lw=1.5, ls=":")
        for a, b in r["segs"]:
            if not r["extent"] or b < r["extent"][0] or a > r["extent"][1]:
                ax.plot([a, b], [y, y], color="k", lw=3)
    ax.set_yticks(range(len(runs), 0, -1),
                  [f"P{r['pos']:02d} r{r['run']}" for r in runs], fontsize=7)
    ax.axvline(0, color="k", lw=0.8)
    ax.plot([], [], color=DIKE_COLOR, lw=6, label=f"dike signal (+ {half:.0f} cm margin)")
    ax.plot([], [], color="0.5", lw=1.5, label="far field")
    ax.plot([], [], color="0.5", lw=1.5, ls=":", label="no dike found: all far field")
    ax.plot([], [], color="k", lw=3, label="other steep stretches (in far field)")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_xlim(-60, 60)
    ax.set_xlabel("Distance from best match to the dike template (cm)", fontsize=12)
    ax.set_title(f"Small dike: far-field cutoffs per run from the slope over {WINDOW_CM:.0f} cm\n"
                 f"signal where |slope| > {K_SIGMA:g}σ = {K_SIGMA * noise:.2f} µT/cm for at "
                 f"least {PERSIST_CM:g} cm; dike found in {len(found)} of {len(runs)} runs",
                 fontsize=12)
    ax.grid(alpha=0.3, axis="x")
    out = FIG_DIR / "small-dike_far_field_cutoffs.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)}")
    plot_means(runs)
    plot_far_value(runs)
    plot_means_grid(runs)
    plot_overlay(runs, noise, aligned=True)
    plot_overlay(runs, noise, aligned=False)
    plot_grid(runs, noise)


def plot_means(runs):
    """Far-field mean on each side per run, and the left - right difference."""
    left_c, right_c = "#2a78d6", "#eb6834"
    fig, (ax, ax_d) = plt.subplots(2, 1, figsize=(17, 10), sharex=True,
                                   height_ratios=[2, 1], layout="constrained")
    idx = np.arange(len(runs))
    lm = np.array([r["sides"]["left"][0] for r in runs])
    ls = np.array([r["sides"]["left"][1] for r in runs])
    rm = np.array([r["sides"]["right"][0] for r in runs])
    rs = np.array([r["sides"]["right"][1] for r in runs])
    for i in idx:
        ax.plot([i - 0.15, i + 0.15], [lm[i], rm[i]], color="0.75", lw=1, zorder=1)
    ax.errorbar(idx - 0.15, lm, yerr=ls, fmt="o", ms=6, color=left_c, ecolor=left_c,
                elinewidth=1, capsize=0, label="left of the dike (mean ± SD)", zorder=2)
    ax.errorbar(idx + 0.15, rm, yerr=rs, fmt="s", ms=6, color=right_c, ecolor=right_c,
                elinewidth=1, capsize=0, label="right of the dike (mean ± SD)", zorder=2)
    igrf = 51.18
    ax.axhline(igrf, color="0.3", ls="--", lw=1, label=f"IGRF model for the site ({igrf} µT)")
    ax.set_ylabel("Far-field mean Btotal (µT)", fontsize=12)
    ax.set_title("Far-field mean on each side of the dike, per run", fontsize=13)
    ax.legend(loc="upper left", fontsize=10, frameon=False)

    diff = lm - rm
    ax_d.bar(idx, diff, width=0.6, color="0.45")
    ax_d.axhline(0, color="0.2", lw=0.8)
    med = np.median(diff)
    ax_d.axhline(med, color="0.2", ls=":", lw=1)
    ax_d.annotate(f"median {med:+.2f} µT", (len(runs) - 0.5, med), xytext=(0, -12),
                  textcoords="offset points", ha="right", fontsize=9, color="0.2")
    ax_d.set_ylabel("Left − right (µT)", fontsize=12)
    ax_d.set_title(f"Difference between the two sides: right higher in "
                   f"{int((diff < 0).sum())} of {len(runs)} runs", fontsize=13)
    ax_d.set_xticks(idx, [f"P{r['pos']:02d} r{r['run']}" for r in runs], rotation=90, fontsize=7)
    for a in (ax, ax_d):
        a.grid(alpha=0.25, axis="y")
        for i in idx[::3][1:]:  # separate positions
            a.axvline(i - 0.5, color="0.9", lw=0.8, zorder=0)
        a.spines[["top", "right"]].set_visible(False)
    ax_d.set_xlim(-0.7, len(runs) - 0.3)
    fig.suptitle(f"Small dike, all {len(runs)} runs: far-field averages (cutoffs from the "
                 f"{WINDOW_CM:.0f} cm slope at {K_SIGMA:g}σ; runs with no dike split at the match)",
                 fontsize=14)
    out = FIG_DIR / "small-dike_far_field_means.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)}")


def plot_overlay(runs, noise, aligned):
    """All runs on one axis, Btotal as recorded; aligned on the dike match, or at
    their recorded distance along the run."""
    fig, ax = plt.subplots(figsize=(16, 9), layout="constrained")
    for r in runs:
        far = np.ones(len(r["x"]), bool) if not r["cut"] else (r["x"] < r["cut"][0]) | (r["x"] > r["cut"][1])
        shift = 0.0 if aligned else r["match_cm"]
        x = r["x"] + shift
        y = r["y"]
        ax.plot(np.where(far, x, np.nan), np.where(far, y, np.nan), color="0.6", lw=0.5, alpha=0.5)
        if r["cut"]:
            dike = ~far
            ax.plot(np.where(dike, x, np.nan), np.where(dike, y, np.nan),
                    color=DIKE_COLOR, lw=1.1, alpha=0.8)
            for c in r["cut"]:
                i = np.argmin(np.abs(x - c - shift))
                ax.plot(x[i], y[i], "o", ms=4, mfc=DIKE_COLOR, mec="k", mew=0.5)
    n = sum(r["cut"] is not None for r in runs)
    ax.plot([], [], color=DIKE_COLOR, lw=2, label=f"dike signal ({n} runs)")
    ax.plot([], [], color="0.6", lw=1, label="far field")
    ax.plot([], [], "o", ms=5, mfc="w", mec="k", label="far-field cutoff")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=10, frameon=False)
    if aligned:
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlim(-60, 60)
        ax.set_xlabel("Distance from best match to the dike template (cm)", fontsize=12)
    else:
        ax.set_xlabel("Distance along run (cm, as recorded)", fontsize=12)
    ax.set_ylabel("Btotal (µT)", fontsize=12)
    ax.set_title(f"Small dike, all {len(runs)} runs: far-field cutoffs from the slope over "
                 f"{WINDOW_CM:.0f} cm (|slope| > {K_SIGMA:g}σ = {K_SIGMA * noise:.2f} µT/cm "
                 f"for ≥ {PERSIST_CM:g} cm)\ndike found in "
                 f"{sum(r['cut'] is not None for r in runs)} runs; grey = far field "
                 "(the whole run where no dike was found)", fontsize=13)
    ax.grid(alpha=0.3)
    out = FIG_DIR / ("small-dike_far_field_cutoffs_overlay"
                     f"{'' if aligned else '_unaligned'}.png")
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)}")


def plot_grid(runs, noise):
    """One panel per run: Btotal along the run, far field shaded, cutoffs marked."""
    positions = sorted({r["pos"] for r in runs})
    cols = sorted({r["run"] for r in runs})
    fig, axes = plt.subplots(len(positions), len(cols), figsize=(15, 2.2 * len(positions)),
                             squeeze=False, layout="constrained")
    for r in runs:
        ax = axes[positions.index(r["pos"]), cols.index(r["run"])]
        d = r["x"] + r["match_cm"]  # back to the run's own distance_cm
        ax.plot(d, r["y"], color="0.55", lw=0.5)
        ax.plot(d, moving_average(d, r["y"], WINDOW_CM), color="k", lw=1)
        if r["cut"]:
            a, b = (v + r["match_cm"] for v in r["cut"])
            ax.axvspan(d.min(), a, color="tab:cyan", alpha=0.15, lw=0)
            ax.axvspan(b, d.max(), color="tab:cyan", alpha=0.15, lw=0)
            ax.axvspan(a, b, color="tab:orange", alpha=0.2, lw=0)
            for v in (a, b):
                ax.axvline(v, color="tab:orange", lw=1.5)
            cut = f"far field < {a:.1f} cm and > {b:.1f} cm"
        else:
            ax.axvspan(d.min(), d.max(), color="tab:cyan", alpha=0.15, lw=0)
            cut = "no dike found: all far field"
        lo, hi = ax.get_ylim()
        for s0, s1 in r["segs"]:  # every steep stretch, along the bottom
            ax.plot([s0 + r["match_cm"], s1 + r["match_cm"]], [lo, lo], color="tab:red",
                    lw=4, solid_capstyle="butt", clip_on=False)
        ax.set_ylim(lo, hi)
        ax.set_title(f"P{r['pos']:02d} run {r['run']}\n{cut}", fontsize=8.5)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel("Distance along run (cm)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("Btotal (µT)", fontsize=9)
    fig.suptitle(f"Small dike, all {len(runs)} runs: far-field cutoffs from the slope over "
                 f"{WINDOW_CM:.0f} cm (|slope| > {K_SIGMA:g}σ = {K_SIGMA * noise:.2f} µT/cm)\n"
                 "blue shading = far field, orange = dike signal + margin, red bars = steep "
                 f"stretches; grey = raw, black = {WINDOW_CM:.0f} cm average", fontsize=12)
    out = FIG_DIR / "small-dike_far_field_cutoffs_grid.png"
    fig.savefig(out, dpi=110)
    print(f"Saved {out.relative_to(REPO)}")


def plot_far_value(runs):
    """One far-field value per run: mean of all far-field samples, SD and SE bars."""
    c = "#2a78d6"
    igrf = 51.18
    idx = np.arange(len(runs))
    mean = np.array([r["far"][0] for r in runs])
    sd = np.array([r["far"][1] for r in runs])
    se = np.array([r["far"][2] for r in runs])
    fig, ax = plt.subplots(figsize=(17, 7), layout="constrained")
    mu, sig = mean.mean(), mean.std(ddof=1)  # across the runs' far-field values
    ax.axhspan(mu - sig, mu + sig, color="0.85", lw=0, zorder=0)
    ax.axhline(mu, color="0.2", lw=1.2)
    for k in (-2, 2):
        ax.axhline(mu + k * sig, color="0.4", ls="--", lw=1)
    ax.errorbar(idx, mean, yerr=sd, fmt="none", ecolor=c, alpha=0.3, elinewidth=6, capsize=0)
    ax.errorbar(idx, mean, yerr=se, fmt="o", ms=5, color=c, ecolor="k", elinewidth=1.2, capsize=3)
    ax.text(0.99, 0.97, f"mean of all {len(runs)} runs: {mu:.2f} µT\n"
            f"SD across runs: {sig:.2f} µT\n2σ band: {mu - 2 * sig:.2f} – {mu + 2 * sig:.2f} µT\n"
            f"IGRF model for the site: {igrf} µT",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.6"))
    ax.grid(alpha=0.25, axis="y")
    for i in idx[::3][1:]:  # separate positions
        ax.axvline(i - 0.5, color="0.9", lw=0.8, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.plot([], [], "o", color=c, label="mean of all far-field samples")
    ax.plot([], [], color="k", lw=1.2, label="± standard error of the mean (4 cm blocks)")
    ax.plot([], [], color=c, alpha=0.3, lw=6, label="± SD of the samples")
    ax.plot([], [], color="0.2", lw=1.2, label="mean of all runs")
    ax.plot([], [], color="0.85", lw=8, label="± 1 SD across runs")
    ax.plot([], [], color="0.4", ls="--", label="± 2 SD across runs")
    ax.legend(loc="upper left", fontsize=9, frameon=False, ncols=2)
    ax.set_ylabel("Far-field Btotal (µT)", fontsize=11)
    ax.set_xticks(idx, [f"P{r['pos']:02d} r{r['run']}" for r in runs], rotation=90, fontsize=7)
    ax.set_xlim(-0.7, len(runs) - 0.3)
    ax.set_title(f"Small dike, all {len(runs)} runs: far-field value per run (every sample "
                 f"outside the dike, both sides together; cutoffs from the {WINDOW_CM:.0f} cm "
                 f"slope at {K_SIGMA:g}σ)", fontsize=13)
    out = FIG_DIR / "small-dike_far_field_value.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)}")


def plot_means_grid(runs):
    """One panel per run: every sample, the dike shaded, far-field mean ± SD per side."""
    left_c, right_c = "#2a78d6", "#eb6834"
    positions = sorted({r["pos"] for r in runs})
    cols = sorted({r["run"] for r in runs})
    fig, axes = plt.subplots(len(positions), len(cols), figsize=(15, 2.4 * len(positions)),
                             squeeze=False, layout="constrained")
    for r in runs:
        ax = axes[positions.index(r["pos"]), cols.index(r["run"])]
        d = r["x"] + r["match_cm"]  # back to the run's own distance_cm
        a, b = (v + r["match_cm"] for v in (r["cut"] or (0.0, 0.0)))
        ax.plot(d, r["y"], ".", ms=1.2, color="0.35", mew=0)
        if r["cut"]:
            ax.axvspan(a, b, color="0.2", alpha=0.15, lw=0)
        else:
            ax.axvline(a, color="0.5", ls=":", lw=1)
        boxes = []
        for side, c, xs, corner in (("left", left_c, (d.min(), a), (0.02, "left")),
                                    ("right", right_c, (b, d.max()), (0.98, "right"))):
            s = r["sides"][side]
            if not s:
                continue
            ax.hlines(s[0], *xs, color=c, lw=1.5, ls="--")
            ax.text(corner[0], 0.95, f"{side}: {s[0]:.2f} ± {s[1]:.2f} µT", transform=ax.transAxes,
                    ha=corner[1], va="top", fontsize=7.5, color="k",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=c, lw=1.2, alpha=0.9))
        hi = r["y"].max()
        lo = r["y"].min()
        ax.set_ylim(lo - 0.08 * (hi - lo), hi + 0.35 * (hi - lo))  # headroom for the boxes
        ax.set_title(f"P{r['pos']:02d} run {r['run']}"
                     + ("" if r["cut"] else " (no dike found: split at match)"), fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel("Distance along run (cm)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("Btotal (µT)", fontsize=9)
    fig.suptitle(f"Small dike, all {len(runs)} runs: every sample, dike shaded (cutoffs from the "
                 f"{WINDOW_CM:.0f} cm slope at {K_SIGMA:g}σ)\nboxes and dashed lines = far-field "
                 "mean ± SD on each side (blue = left, orange = right)", fontsize=12)
    out = FIG_DIR / "small-dike_far_field_means_grid.png"
    fig.savefig(out, dpi=110)
    print(f"Saved {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
