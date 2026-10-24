"""Overlay every regular run (positions 1-20, runs 1-3) of each dike on one graph.

Reads named_data/<dike>/1_positions-01-to-20/ and plots the chosen field
component against the distance_cm column, one figure per dike. Colour marks
the position and line style marks the run (solid 1, dashed 2, dotted 3).

--average collapses the runs at each position into one curve: every run is
interpolated onto a common 0.1 cm distance grid and the runs are averaged
point by point, giving 20 lines per figure.

Usage (from the repo root):
    python analysis/plot_positions_overlay.py                  # Btotal, both dikes
    python analysis/plot_positions_overlay.py --component Bz
    python analysis/plot_positions_overlay.py --smooth 5       # 5 cm moving average
    python analysis/plot_positions_overlay.py --dike small-dike
    python analysis/plot_positions_overlay.py --average        # mean of runs 1-3
    python analysis/plot_positions_overlay.py --all-components # Btotal, Bx, By, Bz in a 2x2
    python analysis/plot_positions_overlay.py --dots           # one dot per sample
"""
import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_profiles_distance import moving_average

NAMED_DIR = REPO / "named_data"
RUN_STYLES = {1: "-", 2: "--", 3: ":"}
COMPONENTS = ["Btotal", "Bx", "By", "Bz"]
GRID_CM = np.linspace(0, 90, 901)  # 0.1 cm spacing for --average


def regular_runs(dike):
    """[(position, run, path)] for a dike's positions 1-20 runs, sorted."""
    with open(NAMED_DIR / "metadata.csv", newline="") as fh:
        rows = [r for r in csv.DictReader(fh)
                if r["dike"] == dike and r["experiment"] == "positions"]
    runs = [(int(r["position"]), int(r["run"]), NAMED_DIR / r["labeled_file"])
            for r in rows]
    return sorted(runs)


def averaged_curves(runs, component):
    """[(position, n_runs, x, mean y)] with each position's runs averaged."""
    by_pos = {}
    for pos, _, path in runs:
        df = load(path)
        by_pos.setdefault(pos, []).append(
            np.interp(GRID_CM, df["distance_cm"], df[component]))
    return [(pos, len(ys), GRID_CM, np.mean(ys, axis=0))
            for pos, ys in sorted(by_pos.items())]


def draw_component(ax, runs, colors, component, args):
    if args.average:
        for pos, n, x, y in averaged_curves(runs, component):
            if args.smooth:
                y = moving_average(x, y, args.smooth)
            ax.plot(x, y, color=colors[pos], lw=1.4,
                    label=f"Position {pos:02d} (mean of {n})")
    else:
        for pos, run, path in runs:
            df = load(path)
            x, y = df["distance_cm"], df[component]
            if args.smooth:
                y = moving_average(x, y, args.smooth)
            if args.dots:
                ax.plot(x, y, color=colors[pos], ls="none", marker="o", ms=1.6,
                        mew=0, label=f"Position {pos:02d}, run {run}")
            else:
                ax.plot(x, y, color=colors[pos], ls=RUN_STYLES.get(run, "-."),
                        lw=1.0 if args.smooth else 0.7,
                        label=f"Position {pos:02d}, run {run}")
    ax.set_xlim(0, 90)
    ax.set_ylabel(f"{component} (µT)", fontsize=13)
    ax.tick_params(labelsize=11)
    ax.grid(alpha=0.3)


def sampling_note(runs):
    """Median sample spacing over all runs, in time and in distance."""
    dts, dxs = [], []
    for _, _, path in runs:
        df = load(path)
        dts.append(df["time"].diff().median())
        dxs.append(df["distance_cm"].diff().median())
    dt, dx = np.median(dts), np.median(dxs)
    return f"one dot per sample: median {dt * 1000:.0f} ms ({1 / dt:.0f} Hz), {dx:.2f} cm apart"


def plot_dike(dike, args):
    runs = regular_runs(dike)
    if not runs:
        raise SystemExit(f"No position runs found for {dike}.")
    positions = sorted({p for p, _, _ in runs})
    colors = dict(zip(positions, plt.cm.turbo(np.linspace(0, 1, len(positions)))))
    components = COMPONENTS if args.all_components else [args.component]

    if args.all_components:
        fig, axes = plt.subplots(2, 2, figsize=(30, 18), sharex=True,
                                 layout="constrained")
        axes = axes.flat
    else:
        fig, ax = plt.subplots(figsize=(20, 11), layout="constrained")
        axes = [ax]
    for ax, component in zip(axes, components):
        draw_component(ax, runs, colors, component, args)
        if args.all_components:
            ax.set_title(component, fontsize=17)
    for ax in axes[-2:] if args.all_components else axes:
        ax.set_xlabel("Distance along profile (cm)", fontsize=13)

    what = ("mean of runs 1-3, 2026-09-19" if args.average
            else f"runs 1-3 ({len(runs)} runs, 2026-09-19)")
    name = "Btotal, Bx, By, Bz" if args.all_components else args.component
    title = (f"{dike.replace('-', ' ').capitalize()}: {name} vs distance, "
             f"positions 1-20, {what}")
    if args.smooth:
        title += f"\n{args.smooth:g} cm overlapping moving average"
    if args.dots and not args.average:
        title += f"\n{sampling_note(runs)}"

    leg_title = ("Position (runs 1-3 averaged)" if args.average else
                 ("Position, run (colour = position)" if args.dots
                 else "Position, run (colour = position, line = run)"))
    leg_kw = dict(ncols=1 if args.average else 3, frameon=False, handlelength=3,
                  fontsize=13 if args.average else 11, title=leg_title,
                  title_fontsize=13 if args.average else 12)
    if args.all_components:
        fig.suptitle(title, fontsize=20)
        handles, labels = axes[0].get_legend_handles_labels()
        leg_kw.update(fontsize=15 if args.average else 13,
                      title_fontsize=16 if args.average else 14)
        leg = fig.legend(handles, labels, loc="outside right upper", **leg_kw)
    else:
        axes[0].set_title(title, fontsize=15)
        leg = axes[0].legend(loc="upper left", bbox_to_anchor=(1.01, 1), **leg_kw)
    for line in leg.get_lines():
        line.set_linewidth(3)
        line.set_markersize(8)

    FIG_DIR.mkdir(exist_ok=True)
    tag = "_averaged" if args.average else ""
    tag += f"_smooth{args.smooth:g}cm" if args.smooth else ""
    tag += "_dots" if args.dots and not args.average else ""
    stem = "all-components" if args.all_components else args.component.lower()
    out = FIG_DIR / f"{dike}_{stem}_distance_positions_overlay{tag}.png"
    fig.savefig(out, dpi=120 if args.all_components else 150)
    plt.close(fig)
    print(f"Saved {out.relative_to(REPO)} ({len(runs)} runs"
          f"{f', {len(positions)} averaged curves' if args.average else ''})")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dike", choices=["small-dike", "large-dike"], default=None,
                   help="plot only one dike (default: both)")
    p.add_argument("--component", default="Btotal", choices=COMPONENTS)
    p.add_argument("--all-components", action="store_true",
                   help="2x2 figure with Btotal, Bx, By and Bz")
    p.add_argument("--smooth", type=float, default=None,
                   help="moving-average window width in cm, e.g. 5")
    p.add_argument("--average", action="store_true",
                   help="average the runs at each position into one curve")
    p.add_argument("--dots", action="store_true",
                   help="draw each sample as a dot instead of a line")
    args = p.parse_args()

    for dike in [args.dike] if args.dike else ["small-dike", "large-dike"]:
        plot_dike(dike, args)


if __name__ == "__main__":
    main()
