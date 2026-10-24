"""Far-field quality checks for all 60 small-dike runs, independent of
whether a run shows the dike.

The far field of a run is everything at least FAR_CM from its best match to
the dike template (match_cm from plot_small_dike_aligned.py, which every run
gets, aligned or not). Repeated rows (the logger writing the same Bx, By, Bz
twice) are dropped first. Checks, each on the far field only:
  level     far-field median Btotal minus the IGRF model value at the
            Middlesex Fells on 2026-09-19 (51.18 µT). The model leaves out
            crustal fields and the phone's own offset, so up to LEVEL_UT is
            normal.
  sides     left far-field level minus right. A run whose two ends disagree
            changed conditions somewhere in the middle.
  sigma     scatter about a straight line fitted to each far-field side
            (robust), compared with the median over all runs.
  step      biggest change in Btotal between the STEP_S before and the
            STEP_S after any far-field sample (medians, so single spikes
            don't count).
  rotation  biggest change in the field vector (Bx, By, Bz) over the same
            STEP_S window where Btotal doesn't change: the phone turned.
            Btotal is unaffected; Bz analyses are.
A run is flagged "bad" for Btotal if level, sides, sigma or step fails, and
"rotated" if only the rotation check fails. A step or sigma failure that
repeats in at least 2 of the 3 runs at a position (and nothing else fails)
is labelled "structure": a real magnetic feature in the far field, which
means the far field there starts further out than FAR_CM, not a bad run.

Usage (from the repo root):
    python analysis/small_dike_quality.py
"""
import csv

import numpy as np

from plot_btotal import FIG_DIR, REPO, load
from plot_positions_overlay import regular_runs
from plot_small_dike_aligned import align
from small_dike_noise import notes, robust_sd

IGRF_UT = 51.18
FAR_CM = 25.0
LEVEL_UT = 1.5
SIDES_UT = 1.0
SIGMA_X = 2.0
STEP_S = 0.2
STEP_UT = 1.0
ROTATION_UT = 2.0


def far_field(r):
    df = load(r["path"])
    keep = np.r_[True, np.diff(df[["Bx", "By", "Bz"]].to_numpy(), axis=0).any(axis=1)]
    df = df[keep]
    x = df["distance_cm"].to_numpy() - r["match_cm"]
    return df, x, np.abs(x) >= FAR_CM


def windowed_change(t, v, far):
    """Largest |median(v after) - median(v before)| over STEP_S, at far-field samples."""
    best = 0.0
    for i in np.flatnonzero(far):
        before = (t >= t[i] - STEP_S) & (t < t[i])
        after = (t > t[i]) & (t <= t[i] + STEP_S)
        if before.sum() < 5 or after.sum() < 5:
            continue
        d = np.median(v[after], axis=0) - np.median(v[before], axis=0)
        best = max(best, float(np.linalg.norm(d)))
    return best


def check(r):
    df, x, far = far_field(r)
    t, b = df["time"].to_numpy(), df["Btotal"].to_numpy()
    vec = df[["Bx", "By", "Bz"]].to_numpy()
    sides, devs = {}, []
    for side, m in (("left", far & (x < 0)), ("right", far & (x > 0))):
        if m.sum() < 20:
            continue
        coef = np.polyfit(x[m], b[m], 1)
        devs.append(b[m] - np.polyval(coef, x[m]))
        sides[side] = np.median(b[m])
    step = windowed_change(t, b, far)
    # Rotation: vector change not matched by a Btotal change.
    rot = 0.0
    for i in np.flatnonzero(far)[::5]:
        before = (t >= t[i] - STEP_S) & (t < t[i])
        after = (t > t[i]) & (t <= t[i] + STEP_S)
        if before.sum() < 5 or after.sum() < 5:
            continue
        dv = np.linalg.norm(np.median(vec[after], 0) - np.median(vec[before], 0))
        db = abs(np.median(b[after]) - np.median(b[before]))
        rot = max(rot, dv - db)
    return dict(level=np.median(b[far]) - IGRF_UT,
                sides=(sides["left"] - sides["right"]) if len(sides) == 2 else np.nan,
                sigma=robust_sd(np.concatenate(devs)) if devs else np.nan,
                step=step, rotation=rot, far_cm=np.ptp(x[far]) if far.any() else 0)


def main():
    results, _ = align(regular_runs("small-dike"), "Btotal")
    note = notes()
    rows = []
    for r in results:
        c = check(r)
        c.update(pos=r["pos"], run=r["run"], path=r["path"], note=note.get(r["path"], ""))
        rows.append(c)
    typical = np.nanmedian([c["sigma"] for c in rows])
    for c in rows:
        fails = []
        if abs(c["level"]) > LEVEL_UT:
            fails.append("level")
        if abs(np.nan_to_num(c["sides"])) > SIDES_UT:
            fails.append("sides")
        if c["sigma"] > SIGMA_X * typical:
            fails.append("sigma")
        if c["step"] > STEP_UT:
            fails.append("step")
        c["fails"] = fails
    # A step or extra scatter that repeats in most runs at one position is a
    # real magnetic feature in the far field (positions 4-5), not a bad run.
    for c in rows:
        shape_fails = {"step", "sigma"} & set(c["fails"])
        repeats = sum(bool({"step", "sigma"} & set(o["fails"]))
                      for o in rows if o["pos"] == c["pos"])
        if shape_fails and set(c["fails"]) <= {"step", "sigma"} and repeats >= 2:
            c["verdict"] = "structure"
        else:
            c["verdict"] = ("bad" if c["fails"] else
                            "rotated" if c["rotation"] > ROTATION_UT else "good")

    print(f"far-field sigma: median {typical:.3f} µT (flag above {SIGMA_X * typical:.3f})")
    print(", ".join(f"{sum(c['verdict'] == v for c in rows)} {v}"
                    for v in ("good", "rotated", "structure", "bad")))
    print("pos run  level  sides  sigma  step  rotation  verdict")
    for c in rows:
        print(f"{c['pos']:3d} {c['run']:3d} {c['level']:+6.2f} {c['sides']:+6.2f} "
              f"{c['sigma']:6.3f} {c['step']:5.2f} {c['rotation']:6.2f}    "
              f"{c['verdict']:9s} {','.join(c['fails'])}  {c['note']}")

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "small-dike_quality.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["position", "run", "level_minus_igrf_uT", "left_minus_right_uT",
                    "far_sigma_uT", "max_step_uT", "max_rotation_uT", "verdict",
                    "fails", "note"])
        for c in rows:
            w.writerow([c["pos"], c["run"], f"{c['level']:.2f}", f"{c['sides']:.2f}",
                        f"{c['sigma']:.3f}", f"{c['step']:.2f}", f"{c['rotation']:.2f}",
                        c["verdict"], " ".join(c["fails"]), c["note"]])
    print(f"Saved {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
