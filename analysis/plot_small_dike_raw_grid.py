"""All 60 small-dike runs as recorded, for looking over by eye. No verdicts.

Two figures, positions down and runs 1-3 across, time on the x-axis (the
distance column assumes constant speed, so time shows how long each run
really took):
  - figures/small-dike_raw_components_grid.png: Bx, By, Bz and Btotal, each
    minus its run median, on each panel's own scale. Shape and tilt show up
    here: if the phone rotates, the components trade off while Btotal stays
    flat.
  - figures/small-dike_raw_btotal_shared.png: Btotal as recorded, every panel
    on the same y-scale, so baseline levels can be compared.
Panel titles give the duration, number of samples, median Btotal, and the
field note if there is one.

Usage (from the repo root):
    python analysis/plot_small_dike_raw_grid.py
"""
import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import COMPONENTS, regular_runs
from small_dike_noise import notes

COLORS = {"Btotal": "k", "Bx": "tab:red", "By": "tab:green", "Bz": "tab:blue"}


def grid(runs):
    positions = sorted({p for p, _, _ in runs})
    cols = sorted({r for _, r, _ in runs})
    fig, axes = plt.subplots(len(positions), len(cols), figsize=(15, 2.1 * len(positions)),
                             squeeze=False, layout="constrained")
    return fig, {(p, r): axes[positions.index(p), cols.index(r)] for p, r, _ in runs}, axes


def title(ax, pos, run, df, note):
    flag = f"\n{note}" if note else ""
    ax.set_title(f"P{pos:02d} run {run}: {df['time'].iloc[-1]:.1f} s, {len(df)} samples, "
                 f"Btotal {df['Btotal'].median():.1f} µT{flag}", fontsize=8.5,
                 color="tab:purple" if note else "k")
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)


def main():
    runs = regular_runs("small-dike")
    note = notes()
    data = {(p, r): load(path) for p, r, path in runs}

    # 1. All four components, run median removed, own scale.
    fig, ax_of, axes = grid(runs)
    for p, r, path in runs:
        df, ax = data[p, r], ax_of[p, r]
        for comp in COMPONENTS:
            ax.plot(df["time"], df[comp] - df[comp].median(), color=COLORS[comp],
                    lw=1.1 if comp == "Btotal" else 0.6, label=comp)
        title(ax, p, r, df, note.get(path, ""))
    for ax in axes[-1]:
        ax.set_xlabel("Time (s)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("minus median (µT)", fontsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper right", ncols=4, fontsize=10)
    fig.suptitle(f"Small dike, all {len(runs)} runs as recorded: Bx, By, Bz, Btotal "
                 "minus each run's median, own scale per panel (2026-09-19)", fontsize=13)
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "small-dike_raw_components_grid.png"
    fig.savefig(out, dpi=110)
    print(f"Saved {out.relative_to(REPO)}")

    # 2. Btotal as recorded, shared scale.
    fig, ax_of, axes = grid(runs)
    lo = min(df["Btotal"].min() for df in data.values())
    hi = max(df["Btotal"].max() for df in data.values())
    typical = np.median([df["Btotal"].median() for df in data.values()])
    for p, r, path in runs:
        df, ax = data[p, r], ax_of[p, r]
        ax.axhline(typical, color="tab:orange", lw=1, ls="--")
        ax.plot(df["time"], df["Btotal"], color="k", lw=0.7)
        ax.set_ylim(lo - 0.5, hi + 0.5)
        title(ax, p, r, df, note.get(path, ""))
    for ax in axes[-1]:
        ax.set_xlabel("Time (s)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("Btotal (µT)", fontsize=8)
    fig.suptitle(f"Small dike, all {len(runs)} runs: Btotal as recorded, same y-scale "
                 f"everywhere (dashed = median level of all runs, {typical:.1f} µT)",
                 fontsize=13)
    out = FIG_DIR / "small-dike_raw_btotal_shared.png"
    fig.savefig(out, dpi=110)
    print(f"Saved {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
