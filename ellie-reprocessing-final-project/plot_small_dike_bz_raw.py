"""All 60 small-dike runs, Bz as recorded vs distance along the run. No
alignment, far field or verdicts: just the data.

Two figures, in ellie-reprocessing-final-project/figures/:
  - small-dike_bz_grid.png: one panel per run (positions down, runs
    1-3 across), each on its own y-scale so the shape shows.
  - small-dike_bz_all_runs.png: every run on one axis, coloured by
    position, so levels can be compared.

Usage (from the repo root):
    python ellie-reprocessing-final-project/plot_small_dike_bz_raw.py
"""
import matplotlib.pyplot as plt
import numpy as np

from common import FIG_DIR, REPO, load, regular_runs

RUN_STYLES = {1: "-", 2: "--", 3: ":"}


def main():
    runs = regular_runs("small-dike")
    data = {(p, r): load(path) for p, r, path in runs}
    positions = sorted({p for p, _ in data})
    cols = sorted({r for _, r in data})
    FIG_DIR.mkdir(exist_ok=True)

    # 1. One panel per run, own scale.
    fig, axes = plt.subplots(len(positions), len(cols), figsize=(15, 2.1 * len(positions)),
                             squeeze=False, layout="constrained")
    for (p, r), df in data.items():
        ax = axes[positions.index(p), cols.index(r)]
        ax.plot(df["distance_cm"], df["Bz"], color="tab:blue", lw=0.6)
        ax.set_title(f"P{p:02d} run {r}: median Bz {df['Bz'].median():.1f} µT", fontsize=8.5)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel("Distance along run (cm)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("Bz (µT)", fontsize=9)
    fig.suptitle(f"Small dike, all {len(data)} runs: Bz as recorded (each panel on its own scale)",
                 fontsize=13)
    out = FIG_DIR / "small-dike_bz_grid.png"
    fig.savefig(out, dpi=110)
    print(f"Saved {out.relative_to(REPO)}")

    # 2. Every run on one axis.
    colors = dict(zip(positions, plt.cm.turbo(np.linspace(0, 1, len(positions)))))
    fig, ax = plt.subplots(figsize=(18, 10), layout="constrained")
    for (p, r), df in data.items():
        ax.plot(df["distance_cm"], df["Bz"], color=colors[p], lw=0.8,
                ls=RUN_STYLES.get(r, "-."), label=f"P{p:02d} run {r}")
    ax.set_xlabel("Distance along run (cm)", fontsize=12)
    ax.set_ylabel("Bz (µT)", fontsize=12)
    ax.set_title(f"Small dike, all {len(data)} runs: Bz as recorded "
                 "(colour = position, line style = run)", fontsize=13)
    ax.grid(alpha=0.3)
    leg = ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), ncols=3, fontsize=8,
                    frameon=False, handlelength=3)
    for line in leg.get_lines():
        line.set_linewidth(2.5)
    out = FIG_DIR / "small-dike_bz_all_runs.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
