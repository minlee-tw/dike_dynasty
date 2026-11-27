"""Find the small dike in each run from Bz alone, and split dike from far field.

Everything here uses Bz only: no Btotal values, Btotal-based quality
verdicts or Btotal results feed into any step.

  1. Template. Runs with a clear Bz anomaly (peak-to-peak > STRONG_UT after a
     SMOOTH_CM average, median removed) are lined up on their dip, averaged,
     then re-aligned by cross-correlation with that average (ITER times). The
     template's centre is the steepest descent between its peak and its dip,
     which is where the dike's edge is, not the maximum.
  2. Match. The template (±HALF_CM around its centre) slides along the whole
     run; the best shift is the one with the highest correlation. The fitted
     scale gives the anomaly size in µT. Match quality:
       strong  correlation >= 0.7 and size >= 2 µT
       weak    correlation >= 0.5 and size >= 0.5 µT
       none    otherwise
  3. Dike vs far field: the slope finds the dike, the template confirms it.
     The slope is a LOESS slope: at each sample, a straight line fitted to
     the samples within ±LOESS_SPAN_CM / 2 with tricube weights (local-linear
     LOESS, no robustness passes), span 6 cm.
     Steep stretches: |slope| > K_SIGMA σ for ≥ PERSIST_CM, GAP_CM gap
     bridging, margin = half the span (3 cm). σ is the robust SD of the
     slope at least NOISE_CM from the centre, pooled over all 60 runs (the
     robust SD keeps bad stretches from inflating it). The dike is confirmed if a steep
     stretch overlaps the template's footprint (where the template is more
     than FOOTPRINT_FRAC of its largest value from zero, about -4 to +12.5 cm
     from the centre). The dike then runs from the confirming stretches
     through the template centre (so a stretch on the recovery side still
     covers the drop), joined with neighbours less than 5 cm away that start
     within the template window (±HALF_CM of the centre), plus the
     3 cm margin. A run with no confirmed dike is far field throughout.
     Runs in NOT_CONFIRMED are set to not confirmed by hand (see there).

Outputs (data, always), in ellie-reprocessing-final-project/data/small-dike/:
  - runs/<run>_bz.csv, one per run: every sample, with time, distance_cm, Bz,
    repeated_sample (the logger wrote the previous Bx, By, Bz again; left out
    of all calculations), distance_from_template_centre_cm, the LOESS slope,
    steep (|slope| above threshold) and region: far_field_left, dike or
    far_field_right. Runs with no confirmed dike are far field throughout,
    split into left and right at the template centre.
  - small-dike_bz_runs_summary.csv: one row per run: template match (class,
    correlation, dip size, centre), dike confirmed, far-field cutoffs (cm along
    the run: far field below the first and above the second), far-field mean,
    SD, SE and number of samples, and whether the run is kept.
  - small-dike_bz_far_field_value.csv: the far-field Bz from the kept runs.
With --figures, also writes the figures used in the presentation to
ellie-reprocessing-final-project/figures/.

Far-field value of a run: the mean of every sample outside the dike (both
sides together), with the SD of those samples and the standard error of the
mean from BLOCK_CM block means (blocks never straddle the dike). A run
with no confirmed dike is far field
throughout, so its value is the mean of the whole run. For comparison only,
the script also prints those runs' values with the template footprint (+ the
margin) cut out instead. Outlier runs: far-field value more than EXCLUDE_SD
SD from the mean of all runs, one pass; the far-field Bz is the mean of the
rest.

Usage (from the repo root):
    python ellie-reprocessing-final-project/small_dike_bz.py
    python ellie-reprocessing-final-project/small_dike_bz.py --figures
"""
import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np

from common import (DATA_DIR, FIG_DIR, GRID_CM, REPO, load, loess_slope, moving_average,
                    regular_runs, repeated_rows, robust_sd, run_name)

OUT_DIR = DATA_DIR / "small-dike"

SMOOTH_CM = 2.0
STRONG_UT = 3.0
ITER = 3
HALF_CM = 15.0
FOOTPRINT_FRAC = 0.2
EXCLUDE_SD = 1.0  # far-field outlier runs: outside 1 SD of all runs, one pass (2026-11-25)
LOESS_SPAN_CM = 6.0  # LOESS slope span (2026-11-25)
K_SIGMA = 2.5        # steep = |slope| above K_SIGMA x slope noise ...
PERSIST_CM = 0.5     # ... for at least this long
NOISE_CM = 30.0      # slope noise from samples at least this far from the centre
GAP_CM = 5.0         # steep stretches closer than this are joined
BLOCK_CM = 4.0       # far-field samples averaged in blocks this long for the SE
# Set to not confirmed by hand (decided 2026-11-25), with the reason.
NOT_CONFIRMED = {
    (18, 2): "weak template match (r 0.58, 0.5 µT); only a steep stretch at the very "
             "edge of the footprint",
}
STEP_CM = GRID_CM[1] - GRID_CM[0]
CLASS_COLORS = {"strong": "tab:green", "weak": "tab:orange", "none": "tab:red"}


def anomaly(x, y):
    """Bz on the common grid, run median removed, SMOOTH_CM average."""
    g = np.interp(GRID_CM, x, y)
    return moving_average(GRID_CM, g - np.median(g), SMOOTH_CM)


def shifted(a, n):
    """a shifted by n samples (positive = later), padded with nan."""
    out = np.full_like(a, np.nan)
    if n >= 0:
        out[n:] = a[:len(a) - n]
    else:
        out[:n] = a[-n:]
    return out


def best_match(a, t):
    """(centre index, correlation, scale) of template window t sliding along a."""
    h = len(t) // 2
    tz = t - t.mean()
    best = (None, -np.inf, 0.0)
    for c in range(h, len(a) - h):
        seg = a[c - h:c + h + 1]
        sz = seg - seg.mean()
        den = np.sqrt((sz ** 2).sum() * (tz ** 2).sum())
        if den == 0:
            continue
        r = (sz * tz).sum() / den
        if r > best[1]:
            best = (c, r, (sz * tz).sum() / (tz ** 2).sum())
    return best


def template_centre(t):
    """Index of the steepest descent between the template's peak and its dip."""
    dip = int(np.nanargmin(t))
    peak = int(np.nanargmax(t[:dip])) if dip > 0 else 0
    d = np.gradient(t)
    return peak + int(np.nanargmin(d[peak:dip + 1]))


def build_template(anoms):
    strong = [a for a in anoms if np.ptp(a) > STRONG_UT]
    ref = int(np.median([np.argmin(a) for a in strong]))
    stack = [shifted(a, ref - int(np.argmin(a))) for a in strong]
    t = np.nanmean(stack, axis=0)
    h = int(round(HALF_CM / STEP_CM))
    for _ in range(ITER):
        c = template_centre(t)
        win = t[c - h:c + h + 1]
        stack = []
        for a in strong:
            ca, _, _ = best_match(a, win)
            stack.append(shifted(a, c - ca))
        t = np.nanmean(stack, axis=0)
    c = template_centre(t)
    return t[c - h:c + h + 1], len(strong)


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


def footprint(t, t_x):
    """(lo, hi) offsets from the centre where the template is clearly away from zero."""
    a = t - np.median(np.r_[t[:20], t[-20:]])
    big = t_x[np.abs(a) > FOOTPRINT_FRAC * np.abs(a).max()]
    return big.min(), big.max()


def confirmed_extent(segs, lo, hi):
    """(left, right) from the stretches overlapping [lo, hi] and the centre, or None.

    Neighbouring stretches less than GAP_CM away are joined on, but only ones
    that start within the template window (±HALF_CM of the centre), so a chain
    of steep stretches can't carry the dike out to unrelated features, such as
    a dip at the start of the run.
    """
    hits = [s for s in segs if s[1] >= lo and s[0] <= hi]
    if not hits:
        return None
    left = min(min(a for a, _ in hits), 0.0)
    right = max(max(b for _, b in hits), 0.0)
    for a, b in sorted(segs):
        if a > right and a - right < GAP_CM and a <= HALF_CM:
            right = b
    for a, b in sorted(segs, reverse=True):
        if b < left and left - b < GAP_CM and b >= -HALF_CM:
            left = a
    return left, right


def far_value(d, y, a, b):
    """(mean, SD, SE, n) of the samples outside [a, b]; SE from BLOCK_CM block means.
    For a run with no dike, pass a == b (any split point)."""
    far = (d < a) | (d > b)
    x, v = d[far], y[far]
    idx = np.where(x < a, np.floor((x - d.min()) / BLOCK_CM),
                   1000 + np.floor((x - b) / BLOCK_CM)).astype(int)
    means = np.array([v[idx == i].mean() for i in np.unique(idx) if (idx == i).sum() >= 3])
    return v.mean(), v.std(ddof=1), means.std(ddof=1) / np.sqrt(len(means)), far.sum()


def classify(r, size):
    if r >= 0.7 and size >= 2.0:
        return "strong"
    if r >= 0.5 and size >= 0.5:
        return "weak"
    return "none"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--figures", action="store_true", help="also write the figures")
    args = p.parse_args()

    runs = []
    for pos, run, path in regular_runs("small-dike"):
        df = load(path)
        rep = repeated_rows(df)
        keep = ~rep
        d = df["distance_cm"].to_numpy()[keep]
        y = df["Bz"].to_numpy()[keep]
        runs.append(dict(pos=pos, run=run, path=path, df=df, rep=rep, d=d, y=y,
                         a=anomaly(d, y)))

    tmpl, n_strong = build_template([r["a"] for r in runs])
    h = len(tmpl) // 2
    t_x = (np.arange(len(tmpl)) - h) * STEP_CM
    lo, hi = footprint(tmpl, t_x)
    print(f"template from {n_strong} runs with Bz peak-to-peak > {STRONG_UT} µT; "
          f"size {np.ptp(tmpl):.1f} µT; footprint {lo:+.1f} to {hi:+.1f} cm from the centre")

    for r in runs:
        c, corr, scale = best_match(r["a"], tmpl)
        r["centre"] = GRID_CM[c]
        r["corr"] = corr
        r["size"] = scale * np.ptp(tmpl)
        r["cls"] = classify(corr, r["size"])
        r["x"] = r["d"] - r["centre"]  # distance from the match centre
        r["fit"] = (scale, np.mean(r["a"][c - h:c + h + 1]) - scale * tmpl.mean())

    # slope noise (all runs) and dike extent
    for r in runs:
        r["s"] = loess_slope(r["x"], r["y"], LOESS_SPAN_CM)
    noise = robust_sd(np.concatenate([r["s"][(np.abs(r["x"]) >= NOISE_CM) & np.isfinite(r["s"])]
                                      for r in runs]))
    half = LOESS_SPAN_CM / 2
    for r in runs:
        r["segs"] = stretches(np.nan_to_num(np.abs(r["s"])) > K_SIGMA * noise, r["x"],
                                 PERSIST_CM)
        ext = confirmed_extent(r["segs"], lo, hi)
        if (r["pos"], r["run"]) in NOT_CONFIRMED:
            ext = None
        r["cut"] = (ext[0] - half + r["centre"], ext[1] + half + r["centre"]) if ext else None

    print(f"slope noise {noise:.3f} µT/cm, threshold {K_SIGMA * noise:.2f} µT/cm")
    for cls in CLASS_COLORS:
        sel = [r for r in runs if r["cls"] == cls]
        found = sum(r["cut"] is not None for r in sel)
        print(f"  {cls:6s} match: {len(sel):2d} runs, dike confirmed in {found:2d}")
    print(f"dike confirmed in {sum(r['cut'] is not None for r in runs)} of {len(runs)} runs")
    print("\nrun     match   corr  size µT  centre   dike (cm along run)")
    for r in runs:
        cut = f"{r['cut'][0]:5.1f}–{r['cut'][1]:5.1f}" if r["cut"] else "not confirmed (all far field)"
        print(f"P{r['pos']:02d} r{r['run']}  {r['cls']:6s} {r['corr']:5.2f}  "
              f"{r['size']:6.2f}  {r['centre']:6.1f}    {cut}")

    # far-field value per run
    for r in runs:
        a, b = r["cut"] if r["cut"] else (r["centre"], r["centre"])
        r["far"] = far_value(r["d"], r["y"], a, b)
    print("\nruns with no confirmed dike: far-field Bz (µT) from the whole run vs with the "
          f"template footprint ({lo:+.1f} to {hi:+.1f} cm) + {half:g} cm margin cut out")
    for r in runs:
        if r["cut"]:
            continue
        alt = far_value(r["d"], r["y"], r["centre"] + lo - half, r["centre"] + hi + half)
        print(f"  P{r['pos']:02d} r{r['run']}: whole run {r['far'][0]:7.2f} ± {r['far'][2]:.3f}   "
              f"cut out {alt[0]:7.2f} ± {alt[2]:.3f}   difference {r['far'][0] - alt[0]:+.3f}")
    vals = np.array([r["far"][0] for r in runs])
    print(f"far-field Bz, all {len(runs)} runs: mean {vals.mean():.2f} µT, SD {vals.std(ddof=1):.2f}, "
          f"median {np.median(vals):.2f}")
    # Outlier runs: far-field value more than EXCLUDE_SD SD from the mean of all runs,
    # in one pass (decided 2026-11-25). These are mostly tilted runs.
    mu, sig = vals.mean(), vals.std(ddof=1)
    for r in runs:
        r["kept"] = abs(r["far"][0] - mu) <= EXCLUDE_SD * sig
    kept = np.array([r["far"][0] for r in runs if r["kept"]])
    print(f"excluded (outside {EXCLUDE_SD:g} SD, {mu - EXCLUDE_SD * sig:.2f} to "
          f"{mu + EXCLUDE_SD * sig:.2f} µT): "
          + ", ".join(f"P{r['pos']:02d} r{r['run']} ({r['far'][0]:.2f})" for r in runs if not r["kept"]))
    print(f"FAR-FIELD Bz, {len(kept)} kept runs: mean {kept.mean():.2f} µT, SD across runs "
          f"{kept.std(ddof=1):.2f}, SE {kept.std(ddof=1) / np.sqrt(len(kept)):.2f}, "
          f"median {np.median(kept):.2f}")

    (OUT_DIR / "runs").mkdir(parents=True, exist_ok=True)
    write_run_files(runs, noise)
    out = OUT_DIR / "small-dike_bz_far_field_value.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "value", "unit", "note"])
        w.writerow(["far_field_bz_mean", f"{kept.mean():.3f}", "uT",
                    f"mean of the far-field values of the {len(kept)} kept runs"])
        w.writerow(["far_field_bz_se", f"{kept.std(ddof=1) / np.sqrt(len(kept)):.3f}", "uT",
                    "standard error of that mean"])
        w.writerow(["far_field_bz_sd_across_runs", f"{kept.std(ddof=1):.3f}", "uT", ""])
        w.writerow(["far_field_bz_median", f"{np.median(kept):.3f}", "uT", ""])
        w.writerow(["runs_kept", len(kept), "", f"outside {EXCLUDE_SD:g} SD of all "
                    f"{len(runs)} runs excluded, one pass"])
        w.writerow(["dikes_confirmed", sum(r["cut"] is not None for r in runs), "", ""])
        w.writerow(["slope_noise_sigma", f"{noise:.4f}", "uT/cm",
                    f"robust SD of the {LOESS_SPAN_CM:g} cm LOESS slope >= {NOISE_CM:g} cm "
                    "from the template centre, all runs"])
    print(f"Saved {out.relative_to(REPO)}")

    out = OUT_DIR / "small-dike_bz_runs_summary.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["position", "run", "match", "correlation", "size_uT", "centre_cm",
                    "dike_confirmed", "far_left_below_cm",
                    "far_right_above_cm", "far_mean_uT", "far_sd_uT", "far_se_uT",
                    "far_samples", "kept"])
        for r in runs:
            cut = [f"{v:.2f}" for v in r["cut"]] if r["cut"] else ["", ""]
            w.writerow([r["pos"], r["run"], r["cls"], f"{r['corr']:.3f}", f"{r['size']:.2f}",
                        f"{r['centre']:.2f}", "yes" if r["cut"] else "no"] + cut
                       + [f"{r['far'][0]:.3f}", f"{r['far'][1]:.3f}", f"{r['far'][2]:.3f}",
                          r["far"][3], "yes" if r["kept"] else "no"])
    print(f"Saved {out.relative_to(REPO)}")

    if args.figures:
        FIG_DIR.mkdir(exist_ok=True)
        plot_grid(runs, tmpl, t_x, noise, (lo, hi))
        plot_aligned(runs, tmpl, t_x)
        plot_template_method(runs, tmpl, t_x, (lo, hi), n_strong)
        plot_template_grid(runs, tmpl, t_x, (lo, hi))
        plot_overlay(runs, noise, aligned=True)
        plot_overlay(runs, noise, aligned=False)
        plot_far_value(runs)


def write_run_files(runs, noise):
    """One CSV per run with every sample (repeated samples flagged), its LOESS slope,
    whether it is steep, and its region."""
    for r in runs:
        df, rep = r["df"], r["rep"]
        dist = df["distance_cm"].to_numpy()
        slope = np.full(len(df), np.nan)
        slope[~rep] = r["s"]
        steep = np.abs(np.nan_to_num(slope)) > K_SIGMA * noise
        if r["cut"]:
            a, b = r["cut"]
            region = np.where(dist < a, "far_field_left",
                              np.where(dist > b, "far_field_right", "dike"))
        else:
            region = np.where(dist < r["centre"], "far_field_left", "far_field_right")
        out = OUT_DIR / "runs" / run_name("small-dike", r["pos"], r["run"], r["path"])
        with open(out, "w", newline="") as fh:
            fh.write(f"# small dike position {r['pos']} run {r['run']} (Bz only): dike "
                     + (f"confirmed; far field below {r['cut'][0]:.2f} cm and above "
                        f"{r['cut'][1]:.2f} cm" if r["cut"] else
                        "not confirmed; far field throughout")
                     + f"; template centre {r['centre']:.2f} cm\n")
            w = csv.writer(fh)
            w.writerow(["time", "distance_cm", "Bz", "repeated_sample",
                        "distance_from_template_centre_cm", "loess_slope_uT_per_cm", "steep",
                        "region"])
            for i in range(len(df)):
                w.writerow([f"{df['time'].iat[i]:.3f}", f"{dist[i]:.2f}", f"{df['Bz'].iat[i]:.2f}",
                            "yes" if rep[i] else "no", f"{dist[i] - r['centre']:.2f}",
                            "" if np.isnan(slope[i]) else f"{slope[i]:.4f}",
                            "yes" if steep[i] else "no", region[i]])
    print(f"Saved {len(runs)} run files to {(OUT_DIR / 'runs').relative_to(REPO)}")


def plot_grid(runs, tmpl, t_x, noise, fp):
    positions = sorted({r["pos"] for r in runs})
    cols = sorted({r["run"] for r in runs})
    fig, axes = plt.subplots(len(positions), len(cols), figsize=(15, 2.3 * len(positions)),
                             squeeze=False, layout="constrained")
    for r in runs:
        ax = axes[positions.index(r["pos"]), cols.index(r["run"])]
        d0, d1 = r["d"].min(), r["d"].max()
        if r["cut"]:
            a, b = r["cut"]
            ax.axvspan(d0, a, color="tab:cyan", alpha=0.15, lw=0)
            ax.axvspan(b, d1, color="tab:cyan", alpha=0.15, lw=0)
            ax.axvspan(a, b, color="tab:orange", alpha=0.25, lw=0)
        else:
            ax.axvspan(d0, d1, color="tab:cyan", alpha=0.15, lw=0)
        ax.plot(r["d"], r["y"], color="0.35", lw=0.5)
        scale, offset = r["fit"]
        base = np.median(np.interp(GRID_CM, r["d"], r["y"]))
        ax.plot(t_x + r["centre"], base + offset + scale * tmpl, color="tab:red", lw=1.2)
        ax.axvline(r["centre"], color=CLASS_COLORS[r["cls"]], lw=1.2, ls="--")
        lo_, hi_ = ax.get_ylim()
        for s0, s1 in r["segs"]:  # every steep stretch, along the bottom
            ax.plot([s0 + r["centre"], s1 + r["centre"]], [lo_, lo_], color="tab:red", lw=4,
                    solid_capstyle="butt", clip_on=False)
        ax.plot([fp[0] + r["centre"], fp[1] + r["centre"]], [hi_, hi_], color="k", lw=3,
                solid_capstyle="butt", clip_on=False)
        ax.set_ylim(lo_, hi_)
        ax.set_title(f"P{r['pos']:02d} run {r['run']}: {r['cls']} match "
                     f"(r {r['corr']:.2f}, {r['size']:.1f} µT)"
                     + ("" if r["cut"] else "; not confirmed, all far field")
                     + (" (by hand)" if (r["pos"], r["run"]) in NOT_CONFIRMED else ""),
                     fontsize=8.5,
                     color=CLASS_COLORS[r["cls"]] if r["cls"] != "strong" else "k")
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel("Distance along run (cm)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("Bz (µT)", fontsize=9)
    n_conf = sum(r["cut"] is not None for r in runs)
    fig.suptitle(f"Small dike, all {len(runs)} runs: dike from Bz alone, slope finds it and the "
                 f"template confirms it ({n_conf} of {len(runs)} confirmed)\n"
                 "grey = Bz as recorded, red = fitted Bz template, dashed = match centre "
                 "(green strong, orange weak, red none); black bar on top = template footprint, "
                 "red bars on the bottom = steep stretches\n"
                 f"orange shading = dike + {LOESS_SPAN_CM / 2:g} cm margin, "
                 f"blue = far field (LOESS slope over {LOESS_SPAN_CM:g} cm > {K_SIGMA:g}σ = "
                 f"{K_SIGMA * noise:.2f} µT/cm)", fontsize=11)
    out = FIG_DIR / "small-dike_bz_dike_search_grid.png"
    fig.savefig(out, dpi=110)
    print(f"Saved {out.relative_to(REPO)}")


def plot_template_method(runs, tmpl, t_x, fp, n_strong):
    """Two panels for presenting the method: (a) the clear runs that build the template,
    lined up, with the template; (b) every run lined up on its best match, with the
    template scaled to each match class's median size. Footprint shaded, centre marked."""
    fig, axes = plt.subplots(1, 2, figsize=(20, 7.5), sharex=True, layout="constrained")
    strong = [r for r in runs if np.ptp(r["a"]) > STRONG_UT]
    panels = (
        (axes[0], strong, f"(a) Building the template: the {n_strong} runs with Bz peak-to-peak "
                          f"> {STRONG_UT:g} µT, lined up"),
        (axes[1], runs, f"(b) Matching: all {len(runs)} runs, each lined up where it best fits "
                        "the template"),
    )
    for ax, sel, title in panels:
        ax.axvspan(*fp, color="tab:orange", alpha=0.12, lw=0)
        for r in sel:
            color = "0.55" if ax is axes[0] else CLASS_COLORS[r["cls"]]
            ax.plot(GRID_CM - r["centre"], r["a"], color=color, lw=0.7, alpha=0.55)
        ax.axvline(0, color="k", lw=1, ls="--")
        ax.axhline(0, color="0.3", lw=0.6)
        ax.set_title(title, fontsize=12, loc="left")
        ax.set_xlabel("Distance from the template centre (cm)", fontsize=11)
        ax.grid(alpha=0.3)
        ax.set_xlim(-30, 30)
    axes[0].plot(t_x, tmpl, color="k", lw=3, label=f"template: mean of the {n_strong} runs")
    axes[0].plot([], [], color="0.55", lw=1.5, label="one run (Bz minus its median, 2 cm average)")
    for cls in CLASS_COLORS:
        sel = [r for r in runs if r["cls"] == cls]
        size = np.median([r["size"] for r in sel])
        axes[1].plot(t_x, tmpl / np.ptp(tmpl) * size, color=CLASS_COLORS[cls], lw=3,
                     path_effects=[__import__("matplotlib.patheffects", fromlist=["x"]).withStroke(
                         linewidth=5, foreground="k")],
                     label=f"{cls} match ({len(sel)} runs): template at their median size, "
                           f"{size:.1f} µT")
    for ax in axes:
        ax.plot([], [], color="k", lw=1, ls="--", label="template centre: steepest point of the drop")
        ax.plot([], [], color="tab:orange", alpha=0.3, lw=8,
                label=f"footprint ({fp[0]:+.1f} to {fp[1]:+.1f} cm): a steep stretch here "
                      "confirms the dike")
        ax.legend(loc="lower left", fontsize=9, frameon=True, framealpha=0.9)
    axes[0].set_ylabel("Bz minus run median, 2 cm average (µT)", fontsize=11)
    fig.suptitle("Small dike Bz: the dike template. Built from the clearest runs (a), then slid "
                 "along every run to find where the dike is (b)", fontsize=14)
    out = FIG_DIR / "small-dike_bz_template_method.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)}")


def plot_template_grid(runs, tmpl, t_x, fp):
    """One panel per run: Bz as recorded with the template fitted where it matches best,
    the footprint shaded and the centre marked. Template matching only, no dike cutoffs."""
    positions = sorted({r["pos"] for r in runs})
    cols = sorted({r["run"] for r in runs})
    fig, axes = plt.subplots(len(positions), len(cols), figsize=(15, 2.3 * len(positions)),
                             squeeze=False, layout="constrained")
    for r in runs:
        ax = axes[positions.index(r["pos"]), cols.index(r["run"])]
        ax.axvspan(r["centre"] + fp[0], r["centre"] + fp[1], color="tab:orange", alpha=0.15, lw=0)
        ax.plot(r["d"], r["y"], color="0.45", lw=0.5)
        scale, offset = r["fit"]
        base = np.median(np.interp(GRID_CM, r["d"], r["y"]))
        ax.plot(t_x + r["centre"], base + offset + scale * tmpl, color=CLASS_COLORS[r["cls"]],
                lw=2)
        ax.axvline(r["centre"], color="k", lw=0.9, ls="--")
        ax.set_title(f"P{r['pos']:02d} run {r['run']}: {r['cls']} match, r = {r['corr']:.2f}, "
                     f"dip {r['size']:.1f} µT", fontsize=8.5,
                     color="k" if r["cls"] == "strong" else CLASS_COLORS[r["cls"]])
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
        ax.set_xlim(0, 90)
    for ax in axes[-1]:
        ax.set_xlabel("Distance along run (cm)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("Bz (µT)", fontsize=9)
    fig.suptitle(f"Small dike, all {len(runs)} runs: where the Bz dike template fits best\n"
                 "grey = Bz as recorded; coloured = template, scaled to the run "
                 "(green strong, orange weak, red none)\ndashed = template centre "
                 f"(steepest drop); orange band = footprint ({fp[0]:+.1f} to {fp[1]:+.1f} cm)",
                 fontsize=12)
    out = FIG_DIR / "small-dike_bz_template_grid.png"
    fig.savefig(out, dpi=110)
    print(f"Saved {out.relative_to(REPO)}")


def plot_aligned(runs, tmpl, t_x):
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5), sharex=True, layout="constrained")
    for ax, cls in zip(axes, CLASS_COLORS):
        sel = [r for r in runs if r["cls"] == cls]
        for r in sel:
            ax.plot(GRID_CM - r["centre"], r["a"], color=CLASS_COLORS[cls], lw=0.7, alpha=0.6)
        if sel:
            size = np.median([r["size"] for r in sel])
            ax.plot(t_x, tmpl / np.ptp(tmpl) * size, color="k", lw=2,
                    label=f"template scaled to the median size ({size:.1f} µT)")
            ax.legend(loc="lower left", fontsize=9, frameon=False)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_title(f"{cls} match: {len(sel)} runs", fontsize=12)
        ax.set_xlabel("Distance from the Bz match centre (cm)", fontsize=11)
        ax.grid(alpha=0.3)
        ax.set_xlim(-45, 45)
    axes[0].set_ylabel(f"Bz minus run median, {SMOOTH_CM:g} cm average (µT)", fontsize=11)
    fig.suptitle("Small dike: every run lined up on its Bz template match", fontsize=13)
    out = FIG_DIR / "small-dike_bz_dike_search_aligned.png"
    fig.savefig(out, dpi=120)
    print(f"Saved {out.relative_to(REPO)}")


def plot_overlay(runs, noise, aligned):
    """Every run's Bz as recorded on one axis: dike orange, far field grey, cutoffs as
    dots; lined up on the template centre, or at the recorded distance."""
    n_conf = sum(r["cut"] is not None for r in runs)
    fig, ax = plt.subplots(figsize=(16, 9), layout="constrained")
    for r in runs:
        shift = r["centre"] if aligned else 0.0
        x = r["d"] - shift
        if r["cut"]:
            dike = (r["d"] >= r["cut"][0]) & (r["d"] <= r["cut"][1])
            ax.plot(np.where(dike, np.nan, x), r["y"], color="0.6", lw=0.5, alpha=0.5)
            ax.plot(np.where(dike, x, np.nan), r["y"], color="tab:orange", lw=1.1, alpha=0.8)
            for c in r["cut"]:
                i = np.argmin(np.abs(r["d"] - c))
                ax.plot(x[i], r["y"][i], "o", ms=4, mfc="tab:orange", mec="k", mew=0.5)
        else:
            ax.plot(x, r["y"], color="0.6", lw=0.5, alpha=0.5)
    ax.plot([], [], color="tab:orange", lw=2, label=f"dike signal ({n_conf} runs)")
    ax.plot([], [], color="0.6", lw=1, label="far field")
    ax.plot([], [], "o", ms=5, mfc="w", mec="k", label="far-field cutoff")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=10, frameon=False)
    if aligned:
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlim(-60, 60)
        ax.set_xlabel("Distance from the Bz template centre (cm)", fontsize=12)
    else:
        ax.set_xlabel("Distance along run (cm, as recorded)", fontsize=12)
    ax.set_ylabel("Bz (µT)", fontsize=12)
    ax.set_title(f"Small dike Bz, all {len(runs)} runs: far-field cutoffs from the LOESS slope over "
                 f"{LOESS_SPAN_CM:g} cm (|slope| > {K_SIGMA:g}σ = {K_SIGMA * noise:.2f} µT/cm "
                 f"for ≥ {PERSIST_CM:g} cm), confirmed by the Bz template\ndike confirmed in "
                 f"{n_conf} runs; grey = far field (the whole run where the dike was not confirmed)",
                 fontsize=13)
    ax.grid(alpha=0.3)
    out = FIG_DIR / f"small-dike_bz_dike_search_overlay{'' if aligned else '_unaligned'}.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)}")


def plot_far_value(runs):
    """One far-field Bz value per run: mean of all far-field samples, SD and SE bars."""
    c = "#2a78d6"
    noaa_z = -46.86  # WMM vertical component at the site, 2026-09-19 (z up = -down)
    idx = np.arange(len(runs))
    mean = np.array([r["far"][0] for r in runs])
    sd = np.array([r["far"][1] for r in runs])
    se = np.array([r["far"][2] for r in runs])
    conf = np.array([r["cut"] is not None for r in runs])
    fig, ax = plt.subplots(figsize=(17, 7), layout="constrained")
    mu, sig = mean.mean(), mean.std(ddof=1)  # across the runs' far-field values
    ax.axhspan(mu - sig, mu + sig, color="0.85", lw=0, zorder=0)
    ax.axhline(mu, color="0.2", lw=1.2)
    for k in (-2, 2):
        ax.axhline(mu + k * sig, color="0.4", ls="--", lw=1)
    ax.axhline(noaa_z, color="tab:red", ls=":", lw=1.2)
    ax.errorbar(idx, mean, yerr=sd, fmt="none", ecolor=c, alpha=0.3, elinewidth=6, capsize=0)
    ax.errorbar(idx, mean, yerr=se, fmt="none", ecolor="k", elinewidth=1.2, capsize=3)
    ax.plot(idx[conf], mean[conf], "o", ms=5, color=c)
    ax.plot(idx[~conf], mean[~conf], "o", ms=5, mfc="white", mec=c, mew=1.5)
    kept = np.array([r["kept"] for r in runs])
    ax.plot(idx[~kept], mean[~kept], "x", ms=11, mew=2, color="tab:red", zorder=5)
    km = mean[kept]
    k_mu, k_sd = km.mean(), km.std(ddof=1)
    ax.axhline(k_mu, color="tab:green", lw=1.6)
    ax.text(0.01, 0.03, f"mean of all {len(runs)} runs: {mu:.2f} µT, SD across runs {sig:.2f} µT\n"
            f"far-field Bz, {kept.sum()} kept runs: {k_mu:.2f} ± {k_sd / np.sqrt(kept.sum()):.2f} µT "
            f"(SE); SD {k_sd:.2f} µT\n"
            f"NOAA vertical field (level phone): {noaa_z} µT",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.6"))
    ax.grid(alpha=0.25, axis="y")
    for i in idx[::3][1:]:  # separate positions
        ax.axvline(i - 0.5, color="0.9", lw=0.8, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.plot([], [], "o", color=c, label="mean of all far-field samples (dike confirmed)")
    ax.plot([], [], "o", mfc="white", mec=c, mew=1.5,
            label="dike not confirmed: mean of the whole run")
    ax.plot([], [], color="k", lw=1.2, label=f"± standard error of the mean ({BLOCK_CM:g} cm blocks)")
    ax.plot([], [], color=c, alpha=0.3, lw=6, label="± SD of the samples")
    ax.plot([], [], color="0.2", lw=1.2, label="mean of all runs")
    ax.plot([], [], color="0.85", lw=8, label="± 1 SD across runs")
    ax.plot([], [], color="0.4", ls="--", label="± 2 SD across runs")
    ax.plot([], [], color="tab:red", ls=":", label="NOAA vertical field, level phone")
    ax.plot([], [], "x", ms=9, mew=2, color="tab:red",
            label=f"excluded: outside {EXCLUDE_SD:g} SD of all runs ({(~kept).sum()} runs)")
    ax.plot([], [], color="tab:green", lw=1.6, label=f"far-field Bz: mean of the {kept.sum()} kept runs")
    ax.legend(loc="upper left", bbox_to_anchor=(0.12, 1), fontsize=9, frameon=False, ncols=2)
    ax.set_ylabel("Far-field Bz (µT)", fontsize=11)
    ax.set_xticks(idx, [f"P{r['pos']:02d} r{r['run']}" for r in runs], rotation=90, fontsize=7)
    ax.set_xlim(-0.7, len(runs) - 0.3)
    ax.set_title(f"Small dike Bz, all {len(runs)} runs: far-field value per run (every sample "
                 f"outside the dike, both sides together; cutoffs from the {LOESS_SPAN_CM:g} cm LOESS "
                 f"slope at {K_SIGMA:g}σ, confirmed by the Bz template)", fontsize=12)
    out = FIG_DIR / "small-dike_bz_far_field_value.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
