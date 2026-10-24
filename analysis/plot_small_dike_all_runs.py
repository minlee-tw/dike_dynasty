"""All 60 small-dike runs on one graph, coloured by far-field quality verdict.

Each run is shifted so its best match to the dike template (match_cm from
plot_small_dike_aligned.py, which every run gets, aligned or not) sits at
0 cm. Verdicts come from figures/small-dike_quality.csv (run
small_dike_quality.py first). Top: Btotal as recorded, with the IGRF value
for the site. Bottom: each run's far-field median removed, so shapes can be
compared.

Usage (from the repo root):
    python analysis/plot_small_dike_all_runs.py
"""
import csv

import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import regular_runs
from plot_small_dike_aligned import align
from small_dike_quality import FAR_CM, IGRF_UT

STYLES = {  # drawn in this order, so bad runs end up on top
    "good": dict(color="tab:blue", lw=0.6, alpha=0.35),
    "rotated": dict(color="tab:green", lw=0.8, alpha=0.6),
    "structure": dict(color="tab:purple", lw=0.9, alpha=0.7),
    "bad": dict(color="tab:red", lw=1.3, alpha=0.9),
}


def main():
    with open(FIG_DIR / "small-dike_quality.csv", newline="") as fh:
        verdict = {(int(r["position"]), int(r["run"])): r["verdict"]
                   for r in csv.DictReader(fh)}
    results, _ = align(regular_runs("small-dike"), "Btotal")

    fig, (ax_raw, ax_rel) = plt.subplots(2, 1, figsize=(16, 12), sharex=True,
                                         layout="constrained")
    for v, style in STYLES.items():
        runs = [r for r in results if verdict[r["pos"], r["run"]] == v]
        for i, r in enumerate(runs):
            df = load(r["path"])
            x = df["distance_cm"].to_numpy() - r["match_cm"]
            y = df["Btotal"].to_numpy()
            far = np.abs(x) >= FAR_CM
            label = f"{v} ({len(runs)})" if i == 0 else None
            ax_raw.plot(x, y, label=label, **style)
            ax_rel.plot(x, y - np.median(y[far]), **style)
            if v == "bad":
                end = np.argmax(x)
                ax_raw.annotate(f"P{r['pos']:02d} r{r['run']}", (x[end], y[end]),
                                xytext=(4, 0), textcoords="offset points",
                                color=style["color"], fontsize=9, va="center")

    ax_raw.axhline(IGRF_UT, color="k", ls="--", lw=1.2,
                   label=f"IGRF model, Middlesex Fells ({IGRF_UT} µT)")
    ax_raw.set_ylabel("Btotal (µT)", fontsize=12)
    ax_raw.set_title("Btotal as recorded", fontsize=13)
    ax_raw.legend(fontsize=11, loc="upper left", bbox_to_anchor=(1.01, 1),
                  title="Far-field quality verdict", frameon=False)
    ax_rel.set_ylabel("Btotal minus far-field median (µT)", fontsize=12)
    ax_rel.set_title("Far-field level removed", fontsize=13)
    ax_rel.set_xlabel("Distance from best match to the dike template (cm)", fontsize=12)
    for ax in (ax_raw, ax_rel):
        for edge in (-FAR_CM, FAR_CM):
            ax.axvline(edge, color="0.5", ls=":", lw=1)
        ax.grid(alpha=0.3)
    ax_rel.set_xlim(-60, 65)
    fig.suptitle(f"Small dike: all {len(results)} runs on one graph (2026-09-19). "
                 f"Dotted lines = ±{FAR_CM:.0f} cm far-field cutoff", fontsize=14)

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "small-dike_btotal_all_runs.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
