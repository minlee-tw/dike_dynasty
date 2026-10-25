"""Plot each position separately with its 3 runs on the same axes.

For each dike this saves:
  - one grid figure with 20 panels (positions 1-20), and
  - 20 individual figures in figures/<dike>_positions/, one per position.
Each panel shows runs 1-3 against distance_cm, coloured by run. Panels have
their own y-axis range so small anomalies stay visible.

Usage (from the repo root):
    python analysis/plot_positions_panels.py                   # Btotal, both dikes
    python analysis/plot_positions_panels.py --component Bz
    python analysis/plot_positions_panels.py --smooth 3        # 3 cm moving average
    python analysis/plot_positions_panels.py --dike small-dike
    python analysis/plot_positions_panels.py --shared-y        # same y-range on every panel
"""
import argparse
import math

import matplotlib.pyplot as plt

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import regular_runs
from plot_profiles_distance import moving_average

RUN_COLORS = {1: "tab:blue", 2: "tab:orange", 3: "tab:green"}


def runs_by_position(dike):
    by_pos = {}
    for pos, run, path in regular_runs(dike):
        by_pos.setdefault(pos, []).append((run, path))
    if not by_pos:
        raise SystemExit(f"No position runs found for {dike}.")
    return by_pos


def draw_position(ax, runs, args, lw):
    for run, path in runs:
        df = load(path)
        x, y = df["distance_cm"], df[args.component]
        if args.smooth:
            y = moving_average(x, y, args.smooth)
        ax.plot(x, y, color=RUN_COLORS.get(run, "0.4"), lw=lw, label=f"Run {run}")
    ax.set_xlim(0, 90)
    ax.grid(alpha=0.3)


def smooth_note(args):
    return f", {args.smooth:g} cm moving average" if args.smooth else ""


def tag(args):
    return f"_smooth{args.smooth:g}cm" if args.smooth else ""


def dike_name(dike):
    return dike.replace("-", " ").capitalize()


def plot_grid(dike, by_pos, args):
    ncols = 4
    nrows = math.ceil(len(by_pos) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.2 * nrows),
                             sharex=True, sharey=args.shared_y, squeeze=False,
                             layout="constrained")
    for ax, (pos, runs) in zip(axes.flat, sorted(by_pos.items())):
        draw_position(ax, runs, args, lw=0.8)
        ax.set_title(f"Position {pos:02d}", fontsize=12)
        ax.tick_params(labelsize=9)
    for ax in axes.flat[len(by_pos):]:
        ax.axis("off")

    handles, labels = axes.flat[0].get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc="outside upper right", ncols=3,
                     fontsize=12, frameon=False)
    for line in leg.get_lines():
        line.set_linewidth(3)
    fig.supxlabel("Distance along profile (cm)", fontsize=13)
    fig.supylabel(f"{args.component} (µT)", fontsize=13)
    fig.suptitle(f"{dike_name(dike)}: {args.component} vs distance by position, "
                 f"runs 1-3, 2026-09-19{smooth_note(args)}", fontsize=15)

    out = FIG_DIR / f"{dike}_{args.component.lower()}_distance_positions_grid{tag(args)}.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Saved {out.relative_to(REPO)}")


def plot_individual(dike, by_pos, args):
    out_dir = FIG_DIR / f"{dike}_positions"
    out_dir.mkdir(parents=True, exist_ok=True)
    for pos, runs in sorted(by_pos.items()):
        fig, ax = plt.subplots(figsize=(10, 5.5), layout="constrained")
        draw_position(ax, runs, args, lw=1.0)
        ax.set_xlabel("Distance along profile (cm)", fontsize=12)
        ax.set_ylabel(f"{args.component} (µT)", fontsize=12)
        ax.set_title(f"{dike_name(dike)}, position {pos:02d}: {args.component} "
                     f"vs distance, runs 1-3{smooth_note(args)}", fontsize=13)
        leg = ax.legend(fontsize=12, frameon=False)
        for line in leg.get_lines():
            line.set_linewidth(3)
        fig.savefig(out_dir / f"{dike}_position-{pos:02d}_"
                    f"{args.component.lower()}{tag(args)}.png", dpi=130)
        plt.close(fig)
    print(f"Saved {len(by_pos)} figures in {out_dir.relative_to(REPO)}/")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dike", choices=["small-dike", "large-dike"], default=None,
                   help="plot only one dike (default: both)")
    p.add_argument("--component", default="Btotal",
                   choices=["Btotal", "Bx", "By", "Bz"])
    p.add_argument("--smooth", type=float, default=None,
                   help="moving-average window width in cm, e.g. 3")
    p.add_argument("--shared-y", action="store_true",
                   help="use the same y-axis range on every panel of the grid")
    args = p.parse_args()

    FIG_DIR.mkdir(exist_ok=True)
    for dike in [args.dike] if args.dike else ["small-dike", "large-dike"]:
        by_pos = runs_by_position(dike)
        plot_grid(dike, by_pos, args)
        plot_individual(dike, by_pos, args)


if __name__ == "__main__":
    main()
