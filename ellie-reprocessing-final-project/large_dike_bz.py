"""Remove the phone's calibration steps from Bz in the large-dike runs.

In the large-dike data the magnetometer re-calibrates while recording, which
shows up as sudden steps (every ~0.41 s or a multiple of it, with no gap in
the recording). These are not real field changes.

Only Bz is used. stitch() walks through each run and, at every sample where
the Bz step is both larger than `thr` µT and much larger than the recent
sample-to-sample noise (10x the median of the previous `k` steps, plus 5 µT),
it shifts that sample and everything after it so Bz continues along the slope
of the previous `k` samples. A recalibration that moves Bz by less than that
is left in (a few µT at most).

The corrected Bz keeps each run's shape, not its absolute level: the level
before the first recalibration is itself off (-171 to +305 µT across runs).
So each run is then shifted (anchored) so that its far end, the median of Bz
over ANCHOR_CM (80-90 cm along the run), equals the small-dike far-field Bz
(read from data/small-dike/small-dike_bz_far_field_value.csv, written by
small_dike_bz.py, and rounded to 0.01 µT: -42.15 µT). The far end is used for
every run because it is about twice as far from the dike as the start (the
dike's steepest drop is at a median 31 cm along the run) and is the flattest
stretch in 54 of 60 runs (median LOESS slope 1.6 µT/cm, against 5 µT/cm over
the first 10 cm). Anchoring the start, or the mean of both ends, gives the
same agreement between repeat runs to within about 1 µT, but would leave the
far, flattest end 70-145 µT off background. Some dike field still remains at
80-90 cm, so the anchored levels are "relative to the far end, taken as
background" (possibly off by tens of µT at P02, P03, P08 and P09).

The original run files in named_data/ are only read, never changed.

Outputs (data, always), in ellie-reprocessing-final-project/data/large-dike/:
  - runs/<run>_bz.csv, one per run: every sample, with time, distance_cm,
    Bz_raw (as recorded), recalibration_step (yes on the first sample after
    each removed step), Bz_step_corrected (levels arbitrary), Bz_anchored
    (shifted so the 80-90 cm median is the small-dike far field) and
    anchor_window (yes for the 80-90 cm samples used to set the level).
  - large-dike_bz_runs_summary.csv: one row per run: steps removed, the
    corrected median over 80-90 cm, the shift applied, the anchored range, and
    the LOESS slope at the far end (how flat the anchor stretch is).
No far field is determined for the large dike itself: its level comes from
the small-dike far-field Bz.

With --figures, also writes the figures to ellie-reprocessing-final-project/figures/:
large-dike_bz_steps_removed.png (raw, steps removed, anchored),
large-dike_bz_steps_removed_grid.png, large-dike_bz_anchored_grid.png and
large-dike_bz_corrected_to_anchored.png.

Usage (from the repo root, after small_dike_bz.py):
    python ellie-reprocessing-final-project/large_dike_bz.py
    python ellie-reprocessing-final-project/large_dike_bz.py --figures
"""
import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np

from common import DATA_DIR, FIG_DIR, REPO, load, loess_slope, regular_runs, run_name

OUT_DIR = DATA_DIR / "large-dike"
SMALL_DIKE_VALUE = DATA_DIR / "small-dike" / "small-dike_bz_far_field_value.csv"
ANCHOR_CM = (80.0, 90.0)     # far end of each run, set to the small-dike far field (median)


def small_dike_far_field():
    """The small-dike far-field Bz (µT), rounded to 0.01 µT, from small_dike_bz.py's output."""
    with open(SMALL_DIKE_VALUE, newline="") as fh:
        rows = {r["quantity"]: r["value"] for r in csv.DictReader(fh)}
    return round(float(rows["far_field_bz_mean"]), 2)


FAR_FIELD_BZ_UT = small_dike_far_field()


def stitch(d, thr=20, k=5):
    """Return (corrected Bz array, list of jump times in s)."""
    z = d["Bz"].to_numpy(dtype=float).copy()
    t = d["time"].to_numpy(dtype=float)
    jumps = []
    for i in range(k + 1, len(z)):
        step = abs(z[i] - z[i - 1])
        prev = np.median(np.abs(np.diff(z[i - k - 1:i])))
        if step > thr and step > 10 * prev + 5:
            slope = (z[i - 1] - z[i - k - 1]) / (t[i - 1] - t[i - k - 1])
            pred = z[i - 1] + slope * (t[i] - t[i - 1])
            z[i:] += pred - z[i]
            jumps.append(t[i])
    return z, jumps


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--figures", action="store_true", help="also write the figures")
    args = p.parse_args()

    runs = []
    for pos, run, f in regular_runs("large-dike"):
        d = load(f)
        z, jumps = stitch(d)
        t = d["time"].to_numpy()
        x = d["distance_cm"].to_numpy()
        anchor = (x >= ANCHOR_CM[0]) & (x <= ANCHOR_CM[1])
        shift = FAR_FIELD_BZ_UT - np.median(z[anchor])
        runs.append(dict(pos=pos, run=run, path=f, t=t, x=x, raw=d["Bz"].to_numpy(), z=z,
                         za=z + shift, shift=shift, anchor=anchor,
                         jump_i=np.searchsorted(t, jumps), jump_x=x[np.searchsorted(t, jumps)]))
        first = f"{jumps[0]:5.2f}s" if jumps else "  -  "
        print(f"pos {pos:02d} run {run}: {len(jumps):2d} steps, first {first}, "
              f"corrected min {z.min():6.0f} max {z.max():6.0f}, anchored "
              f"{(z + shift).min():6.0f} to {(z + shift).max():6.0f}")
    write_outputs(runs)
    if args.figures:
        FIG_DIR.mkdir(exist_ok=True)
        plot_overview(runs)
        plot_anchored_grid(runs)
        plot_corrected_to_anchored(runs)


def write_outputs(runs):
    """Per-run CSVs with every sample, and a one-row-per-run summary."""
    (OUT_DIR / "runs").mkdir(parents=True, exist_ok=True)
    for r in runs:
        step = np.zeros(len(r["x"]), bool)
        step[r["jump_i"]] = True
        out = OUT_DIR / "runs" / run_name("large-dike", r["pos"], r["run"], r["path"])
        with open(out, "w", newline="") as fh:
            fh.write(f"# large dike position {r['pos']} run {r['run']} (Bz only): "
                     f"{len(r['jump_i'])} recalibration steps removed; anchored so the median "
                     f"over {ANCHOR_CM[0]:g}-{ANCHOR_CM[1]:g} cm is {FAR_FIELD_BZ_UT} uT "
                     f"(shift {r['shift']:+.2f} uT)\n")
            w = csv.writer(fh)
            w.writerow(["time", "distance_cm", "Bz_raw", "recalibration_step", "Bz_step_corrected",
                        "Bz_anchored", "anchor_window"])
            for i in range(len(r["x"])):
                w.writerow([f"{r['t'][i]:.3f}", f"{r['x'][i]:.2f}", f"{r['raw'][i]:.2f}",
                            "yes" if step[i] else "no", f"{r['z'][i]:.2f}", f"{r['za'][i]:.2f}",
                            "yes" if r["anchor"][i] else "no"])
    print(f"Saved {len(runs)} run files to {(OUT_DIR / 'runs').relative_to(REPO)}")

    out = OUT_DIR / "large-dike_bz_runs_summary.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["position", "run", "steps_removed", "corrected_median_80_90cm_uT", "shift_uT",
                    "anchored_min_uT", "anchored_max_uT", "far_end_loess_slope_uT_per_cm"])
        for r in runs:
            s = loess_slope(r["x"], r["z"], 6.0)
            w.writerow([r["pos"], r["run"], len(r["jump_i"]),
                        f"{np.median(r['z'][r['anchor']]):.2f}", f"{r['shift']:.2f}",
                        f"{r['za'].min():.2f}", f"{r['za'].max():.2f}",
                        f"{np.nanmedian(np.abs(s[(r['x'] >= 77) & (r['x'] <= 87)])):.2f}"])
    print(f"Saved {out.relative_to(REPO)}")


def plot_overview(runs):
    """Raw, steps removed, anchored; and one panel per run, steps marked."""
    # 1. All runs: raw, corrected, anchored.
    cmap = plt.get_cmap("turbo")
    fig, ax = plt.subplots(1, 3, figsize=(26, 7.5), sharex=True)
    for r in runs:
        c = cmap((r["pos"] - 1) / 19 * 0.9 + 0.05)
        ax[0].plot(r["x"], r["raw"], ".", ms=1, color=c)
        ax[1].plot(r["x"], r["z"], ".", ms=1, color=c)
        ax[2].plot(r["x"], r["za"], ".", ms=1, color=c,
                   label=f"Pos {r['pos']:02d}" if r["run"] == 1 else None)
    ax[2].axvspan(*ANCHOR_CM, color="0.5", alpha=0.15, lw=0)
    ax[2].axhline(FAR_FIELD_BZ_UT, color="k", ls="--", lw=1)
    ax[0].set_title("1. Bz, raw")
    ax[1].set_title("2. Calibration steps removed (levels arbitrary)")
    ax[2].set_title(f"3. Anchored: median over {ANCHOR_CM[0]:g}-{ANCHOR_CM[1]:g} cm (grey) set to the "
                    f"small-dike far field, {FAR_FIELD_BZ_UT} µT (dashed)")
    for a in ax:
        a.set_xlabel("Distance along profile (cm)")
        a.set_ylabel("Bz (µT)")
        a.grid(alpha=0.3)
        a.set_xlim(0, 90)
    ax[2].legend(fontsize=8, markerscale=10, ncols=2, loc="upper right")
    fig.suptitle(f"Large dike: positions 1-20, runs 1-3 ({len(runs)} runs), "
                 "colour = position", fontsize=14)
    out = FIG_DIR / "large-dike_bz_steps_removed.png"
    fig.savefig(out, dpi=110, bbox_inches="tight")
    print(f"Saved {out.relative_to(REPO)}")

    # 2. One panel per run.
    positions = sorted({r["pos"] for r in runs})
    cols = sorted({r["run"] for r in runs})
    fig, axes = plt.subplots(len(positions), len(cols), figsize=(15, 2.1 * len(positions)),
                             squeeze=False, layout="constrained")
    for r in runs:
        ax = axes[positions.index(r["pos"]), cols.index(r["run"])]
        ax.plot(r["x"], r["z"], color="tab:blue", lw=0.6)
        lo, hi = ax.get_ylim()
        ax.plot(r["jump_x"], np.full(len(r["jump_x"]), lo), "|", color="tab:red", ms=6,
                clip_on=False)
        ax.set_ylim(lo, hi)
        ax.set_title(f"P{r['pos']:02d} run {r['run']}: {len(r['jump_x'])} steps removed",
                     fontsize=8.5)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel("Distance along profile (cm)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("Bz (µT)", fontsize=9)
    fig.suptitle(f"Large dike, all {len(runs)} runs: Bz with calibration steps removed "
                 "(each panel on its own scale; red ticks = removed steps)", fontsize=13)
    out = FIG_DIR / "large-dike_bz_steps_removed_grid.png"
    fig.savefig(out, dpi=110)
    print(f"Saved {out.relative_to(REPO)}")


def plot_anchored_grid(runs):
    """One panel per run: anchored Bz, anchor stretch shaded, small-dike far field dashed."""
    positions = sorted({r["pos"] for r in runs})
    cols = sorted({r["run"] for r in runs})
    fig, axes = plt.subplots(len(positions), len(cols), figsize=(15, 2.1 * len(positions)),
                             squeeze=False, layout="constrained")
    for r in runs:
        ax = axes[positions.index(r["pos"]), cols.index(r["run"])]
        ax.axvspan(*ANCHOR_CM, color="0.5", alpha=0.15, lw=0)
        ax.axhline(FAR_FIELD_BZ_UT, color="k", ls="--", lw=0.8)
        ax.plot(r["x"], r["za"], color="tab:blue", lw=0.6)
        ax.set_title(f"P{r['pos']:02d} run {r['run']}: shifted {r['shift']:+.0f} µT", fontsize=8.5)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
        ax.set_xlim(0, 90)
    for ax in axes[-1]:
        ax.set_xlabel("Distance along profile (cm)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("Bz (µT)", fontsize=9)
    fig.suptitle(f"Large dike, all {len(runs)} runs: Bz with calibration steps removed, anchored so "
                 f"the median over {ANCHOR_CM[0]:g}-{ANCHOR_CM[1]:g} cm (grey) is the small-dike "
                 f"far field, {FAR_FIELD_BZ_UT} µT (dashed)", fontsize=12)
    out = FIG_DIR / "large-dike_bz_anchored_grid.png"
    fig.savefig(out, dpi=110)
    print(f"Saved {out.relative_to(REPO)}")


def plot_corrected_to_anchored(runs):
    """Two panels: steps removed (levels arbitrary), then anchored at the far end."""
    cmap = plt.get_cmap("turbo")
    fig, ax = plt.subplots(1, 2, figsize=(20, 7.5), sharex=True, layout="constrained")
    for r in runs:
        c = cmap((r["pos"] - 1) / 19 * 0.9 + 0.05)
        ax[0].plot(r["x"], r["z"], color=c, lw=0.7)
        ax[1].plot(r["x"], r["za"], color=c, lw=0.7,
                   label=f"Pos {r['pos']:02d}" if r["run"] == 1 else None)
    for a in ax:
        a.axvspan(*ANCHOR_CM, color="0.5", alpha=0.18, lw=0)
        a.set_xlabel("Distance along profile (cm)", fontsize=12)
        a.set_ylabel("Bz (µT)", fontsize=12)
        a.grid(alpha=0.3)
        a.set_xlim(0, 90)
    ax[1].axhline(FAR_FIELD_BZ_UT, color="k", ls="--", lw=1.2,
                  label=f"small-dike far field, {FAR_FIELD_BZ_UT} µT")
    spread = [np.median(r["z"][(r["x"] >= ANCHOR_CM[0]) & (r["x"] <= ANCHOR_CM[1])]) for r in runs]
    ax[0].set_title(f"Calibration steps removed: levels arbitrary\n(far ends spread over "
                    f"{min(spread):.0f} to {max(spread):.0f} µT)", fontsize=13)
    ax[1].set_title(f"Anchored: each run shifted so its median over {ANCHOR_CM[0]:g}-{ANCHOR_CM[1]:g} cm "
                    f"(grey)\nequals the small-dike far-field Bz, {FAR_FIELD_BZ_UT} µT (dashed)",
                    fontsize=13)
    ax[1].legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    fig.suptitle(f"Large dike, all {len(runs)} runs: Bz from step-corrected to anchored "
                 "(colour = position)", fontsize=14)
    out = FIG_DIR / "large-dike_bz_corrected_to_anchored.png"
    fig.savefig(out, dpi=120)
    print(f"Saved {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
