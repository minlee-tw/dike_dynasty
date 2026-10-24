"""Plot every large-dike position 10 run on one graph.

By default this plots the 15 runs taken under the same conditions: runs 1-3
from the positions survey plus repeats 4-15. --all adds the 9 scan-speed runs
and the 11 phone-height runs (35 in total), coloured by experiment.

--speeds compares scan speeds with 3 runs per speed: all 3 runs of the
3, 5 and 15 s test passes, and 3 representative runs of the ~10 s baseline
(the 15 same-condition runs). The representatives are the 3 runs whose
profiles, on a common 0.1 cm grid, have the smallest RMS difference from the
baseline's point-by-point median profile.

--heights compares the 11 phone-height runs (1.1 to 100 cm above the rock),
with a linear y-axis on top and a log y-axis below so the weak high runs
stay readable.

Usage (from the repo root):
    python analysis/plot_position10.py                  # 15 same-condition runs, Btotal
    python analysis/plot_position10.py --all            # all 35 position 10 runs
    python analysis/plot_position10.py --component Bz
    python analysis/plot_position10.py --dots           # one dot per sample
    python analysis/plot_position10.py --speeds         # 3 runs per scan speed
    python analysis/plot_position10.py --heights        # one run per phone height
"""
import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import COMPONENTS, GRID_CM, NAMED_DIR

# Colour map per experiment, so the groups stay distinguishable with --all.
GROUP_CMAPS = {"same": plt.cm.turbo, "speeds": plt.cm.Greys,
               "heights": plt.cm.PuRd}


def position10_runs(include_all):
    """[(group, label, path)] in plotting order."""
    with open(NAMED_DIR / "metadata.csv", newline="") as fh:
        rows = [r for r in csv.DictReader(fh)
                if r["dike"] == "large-dike" and r["position"] == "10"]
    same = sorted((r for r in rows if r["experiment"] in ("positions", "repeats")),
                  key=lambda r: int(r["run"]))
    runs = [("same", f"Run {int(r['run']):02d}", r) for r in same]
    if include_all:
        speeds = sorted((r for r in rows if r["experiment"] == "speeds"),
                        key=lambda r: (int(r["scan_time_s"]), int(r["run"])))
        runs += [("speeds", f"{int(r['scan_time_s']):2d} s pass, run {r['run']}", r)
                 for r in speeds]
        heights = sorted((r for r in rows if r["experiment"] == "heights"),
                         key=lambda r: float(r["height_cm"]))
        runs += [("heights", f"Height {float(r['height_cm']):g} cm", r)
                 for r in heights]
    return [(g, label, NAMED_DIR / r["labeled_file"]) for g, label, r in runs]


def speed_groups():
    """{nominal seconds per pass: [(run, path)]}, with the baseline as 10 s."""
    groups = {}
    with open(NAMED_DIR / "metadata.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["dike"] != "large-dike" or r["position"] != "10":
                continue
            if r["experiment"] in ("positions", "repeats"):
                secs = 10
            elif r["experiment"] == "speeds":
                secs = int(r["scan_time_s"])
            else:
                continue
            groups.setdefault(secs, []).append(
                (int(r["run"]), NAMED_DIR / r["labeled_file"]))
    return dict(sorted(groups.items()))


def most_representative(runs, component, n):
    """The n (run, path) closest (RMS) to the group median, sorted by run."""
    if len(runs) <= n:
        return sorted(runs)
    curves = []
    for _, path in runs:
        df = load(path)
        curves.append(np.interp(GRID_CM, df["distance_cm"], df[component]))
    curves = np.array(curves)
    rms = np.sqrt(((curves - np.median(curves, axis=0)) ** 2).mean(axis=1))
    return sorted(runs[i] for i in np.argsort(rms)[:n])


def plot_speeds(args):
    groups = speed_groups()
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(groups)))
    fig, ax = plt.subplots(figsize=(20, 11), layout="constrained")
    for (secs, runs), color in zip(groups.items(), colors):
        picked = most_representative(runs, args.component, 3)
        for (run, path), ls in zip(picked, ["-", "--", ":"]):
            df = load(path)
            group = " of 15 (baseline)" if secs == 10 else ""
            label = (f"{secs} s pass, run {run}{group}: "
                     f"{df['time'].iloc[-1]:.1f} s recorded")
            style = (dict(ls="none", marker="o", ms=2.5, mew=0) if args.dots
                     else dict(lw=1.4, ls=ls))
            ax.plot(df["distance_cm"], df[args.component], color=color,
                    label=label, **style)
        print(f"{secs:2d} s: runs {', '.join(str(r) for r, _ in picked)}")

    ax.set_xlim(0, 90)
    ax.set_xlabel("Distance along profile (cm)", fontsize=13)
    ax.set_ylabel(f"{args.component} (µT)", fontsize=13)
    ax.set_title(f"Large dike, position 10: {args.component} vs distance at "
                 f"different scan speeds (3 runs per speed, 2026-09-19)", fontsize=15)
    ax.tick_params(labelsize=11)
    ax.grid(alpha=0.3)
    leg = ax.legend(loc="upper left", fontsize=13, frameon=True, framealpha=0.9,
                    title="Time to cross 90 cm (10 s: the 3 runs closest "
                    "to the baseline median)", title_fontsize=13, handlelength=3)
    for line in leg.get_lines():
        line.set_linewidth(4)
        line.set_markersize(10)

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / (f"large-dike_position-10_{args.component.lower()}_speeds"
                     f"{'_dots' if args.dots else ''}.png")
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(REPO)}")


def plot_heights(args):
    runs = [(label, path) for g, label, path in position10_runs(True)
            if g == "heights"]
    colors = plt.cm.plasma(np.linspace(0, 0.9, len(runs)))
    fig, (ax_lin, ax_log) = plt.subplots(2, 1, figsize=(20, 14), sharex=True,
                                         layout="constrained")
    for (label, path), color in zip(runs, colors):
        df = load(path)
        style = (dict(ls="none", marker="o", ms=2.0, mew=0) if args.dots
                 else dict(lw=1.4))
        for ax in (ax_lin, ax_log):
            ax.plot(df["distance_cm"], df[args.component].abs() if ax is ax_log
                    else df[args.component], color=color, label=label, **style)

    ax_log.set_yscale("log")
    ax_lin.set_ylabel(f"{args.component} (µT)", fontsize=13)
    ylog = args.component if args.component == "Btotal" else f"|{args.component}|"
    ax_log.set_ylabel(f"{ylog} (µT, log scale)", fontsize=13)
    ax_log.set_xlabel("Distance along profile (cm)", fontsize=13)
    ax_lin.set_title("Linear scale", fontsize=13)
    ax_log.set_title("Log scale", fontsize=13)
    for ax in (ax_lin, ax_log):
        ax.set_xlim(0, 90)
        ax.tick_params(labelsize=11)
        ax.grid(alpha=0.3, which="both")
    fig.suptitle(f"Large dike, position 10: {args.component} vs distance at "
                 f"different phone heights above the rock ({len(runs)} runs, "
                 f"2026-09-19)", fontsize=16)
    handles, labels = ax_lin.get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc="outside right upper", fontsize=13,
                     frameon=False, title="Phone height", title_fontsize=14,
                     handlelength=3)
    for line in leg.get_lines():
        line.set_linewidth(4)
        line.set_markersize(10)

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / (f"large-dike_position-10_{args.component.lower()}_heights"
                     f"{'_dots' if args.dots else ''}.png")
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(REPO)}")


def colours(runs):
    out = []
    for group, cmap in GROUP_CMAPS.items():
        n = sum(g == group for g, _, _ in runs)
        lo = 0.0 if group == "same" else 0.35  # skip the palest shades
        out += list(cmap(np.linspace(lo, 1, n))) if n else []
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--all", action="store_true",
                   help="also plot the scan-speed and phone-height runs")
    p.add_argument("--component", default="Btotal", choices=COMPONENTS)
    p.add_argument("--dots", action="store_true",
                   help="draw each sample as a dot instead of a line")
    p.add_argument("--speeds", action="store_true",
                   help="one representative run per scan speed (3, 5, 10, 15 s)")
    p.add_argument("--heights", action="store_true",
                   help="one run per phone height above the rock")
    args = p.parse_args()

    if args.heights:
        plot_heights(args)
        return
    if args.speeds:
        plot_speeds(args)
        return
    runs = position10_runs(args.all)
    fig, ax = plt.subplots(figsize=(20, 11), layout="constrained")
    for (group, label, path), color in zip(runs, colours(runs)):
        df = load(path)
        style = (dict(ls="none", marker="o", ms=1.8, mew=0) if args.dots
                 else dict(lw=1.0, ls="-" if group == "same" else "--"))
        ax.plot(df["distance_cm"], df[args.component], color=color,
                label=label, **style)

    ax.set_xlim(0, 90)
    ax.set_xlabel("Distance along profile (cm)", fontsize=13)
    ax.set_ylabel(f"{args.component} (µT)", fontsize=13)
    what = ("all position 10 runs (same-condition, speed and height tests)"
            if args.all else "runs 1-15 under the same conditions")
    ax.set_title(f"Large dike, position 10: {args.component} vs distance, {what} "
                 f"({len(runs)} runs, 2026-09-19)", fontsize=15)
    ax.tick_params(labelsize=11)
    ax.grid(alpha=0.3)
    leg = ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1),
                    ncols=2 if args.all else 1, fontsize=13, frameon=False,
                    title="Run", title_fontsize=14, handlelength=3)
    for line in leg.get_lines():
        line.set_linewidth(3)
        line.set_markersize(8)

    FIG_DIR.mkdir(exist_ok=True)
    tag = ("_all" if args.all else "_same-conditions") + ("_dots" if args.dots else "")
    out = FIG_DIR / f"large-dike_position-10_{args.component.lower()}{tag}.png"
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(REPO)} ({len(runs)} runs)")


if __name__ == "__main__":
    main()
