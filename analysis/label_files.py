"""Copy the field runs into named_data/ with descriptive names.

The originals in magnetic_field_data/ are never modified, and
named_data/metadata.csv maps every labeled file back to its original.

Each name ends with the original recording time, e.g.
large-dike-position-01-run-1_12-36-58.csv.

Each copy gets a distance_cm column right after time. Every run covered
the same 90 cm at a constant speed, so time maps linearly onto distance:
distance_cm = 90 * (t - t_first) / (t_last - t_first).

Layout (folders and names sort in this order):
  large-dike/
    1_positions-01-to-20/      large-dike-position-01-run-1 ... position-20-run-3
    2_position-10-runs/        large-dike-position-10-run-04 ... run-15
    3_position-10-speeds/      large-dike-position-10-03sec-speed-run-1 ... 15sec-speed-run-3
    4_position-10-elevations/  large-dike-position-10-elevation-001.1cm ... 100cm
  small-dike/
    1_positions-01-to-20/      small-dike-position-01-run-1 ... position-20-run-3

Experiments (from the field notes). Both dikes are divided into 10 cm
intervals numbered 1-20, and each position was scanned 3 times.
  Small dike: run 1 is the first transect (positions 1->20); runs 2 and 3
    were taken back to back on the return pass (positions 20->1).
  Large dike:
    1. Positions 1-20, 3 identical passes each
    2. Position 10 repeated under the same conditions
    3. Position 10 at different scan speeds (seconds to cross the 90 cm)
    4. Position 10 at different phone heights above the rock

Usage (from the repo root):
    python analysis/label_files.py
"""
import csv
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "magnetic_field_data"
OUT_DIR = REPO / "named_data"
RUN_LENGTH_CM = 90

# ---- Small dike (11:09:51 to 12:02:52) ----
# Position -> runs 1, 2, 3
SMALL_POSITIONS = {
    1: ["11-09-51", "12-02-25", "12-02-52"],
    2: ["11-13-45", "12-01-33", "12-01-58"],
    3: ["11-15-24", "12-00-34", "12-00-58"],
    4: ["11-16-39", "11-59-36", "12-00-02"],
    5: ["11-17-28", "11-58-45", "11-59-12"],
    6: ["11-18-00", "11-57-53", "11-58-15"],
    7: ["11-18-40", "11-57-08", "11-57-29"],
    8: ["11-19-26", "11-56-14", "11-56-44"],
    9: ["11-20-19", "11-54-34", "11-54-57"],
    10: ["11-22-19", "11-53-49", "11-54-12"],
    11: ["11-22-43", "11-52-59", "11-53-24"],
    12: ["11-23-12", "11-52-09", "11-52-31"],
    13: ["11-23-47", "11-51-22", "11-51-43"],
    14: ["11-24-16", "11-50-28", "11-50-55"],
    15: ["11-24-52", "11-40-59", "11-41-43"],
    16: ["11-25-35", "11-39-32", "11-40-13"],
    17: ["11-26-35", "11-38-39", "11-39-05"],
    18: ["11-27-30", "11-35-40", "11-38-03"],
    19: ["11-28-02", "11-34-10", "11-34-44"],
    20: ["11-28-30", "11-32-47", "11-33-21"],
}

# ---- Large dike (12:36:58 to 13:26:39) ----
# 1. Position -> runs 1, 2, 3
LARGE_POSITIONS = {
    1: ["12-36-58", "12-37-53", "12-38-42"],
    2: ["12-39-14", "12-39-40", "12-40-07"],
    3: ["12-40-30", "12-40-53", "12-41-23"],
    4: ["12-41-52", "12-42-17", "12-42-37"],
    5: ["12-43-00", "12-43-24", "12-43-45"],
    6: ["12-44-16", "12-44-37", "12-45-03"],
    7: ["12-45-37", "12-45-59", "12-46-23"],
    8: ["12-46-57", "12-47-19", "12-47-42"],
    9: ["12-48-06", "12-48-31", "12-48-53"],
    10: ["12-49-22", "12-49-46", "12-50-09"],
    11: ["12-50-34", "12-50-55", "12-51-15"],
    12: ["12-51-39", "12-52-01", "12-52-25"],
    13: ["12-52-50", "12-53-09", "12-53-31"],
    14: ["12-53-53", "12-54-14", "12-54-35"],  # notes said 12-54-25 (typo)
    15: ["12-55-06", "12-55-27", "12-55-50"],
    16: ["12-56-18", "12-56-39", "12-56-58"],
    17: ["12-57-19", "12-57-39", "12-58-08"],
    18: ["12-58-35", "12-58-55", "12-59-17"],
    19: ["12-59-52", "13-00-11", "13-00-34"],
    20: ["13-00-58", "13-01-22", "13-01-51"],
}

# 2. Position 10, repeated under the same conditions
LARGE_REPEATS = ["13-18-20", "13-18-44", "13-19-10", "13-19-31", "13-19-53",
                 "13-20-14", "13-20-41", "13-21-04", "13-21-24", "13-21-45",
                 "13-22-06", "13-22-27"]

# 3. Position 10, seconds taken to cross the 90 cm
LARGE_SPEEDS_S = {
    3: ["13-26-03", "13-26-18", "13-26-39"],
    5: ["13-23-05", "13-23-29", "13-23-49"],
    15: ["13-24-35", "13-25-03", "13-25-41"],
}

# 4. Position 10, phone height above the rock (cm)
LARGE_HEIGHTS_CM = {
    "13-07-13": 1.1, "13-09-14": 2, "13-10-40": 3, "13-11-17": 6,
    "13-11-45": 9, "13-12-09": 12, "13-12-53": 15, "13-13-36": 30,
    "13-14-17": 50, "13-15-16": 70, "13-16-16": 100,
}

FRICTION = "phone had some friction on the rock (field notes)"
NOTES = {
    # small dike: first measurement at positions 1, 2 and 12
    "11-09-51": FRICTION,
    "11-13-45": FRICTION,
    "11-23-12": FRICTION,
    # large dike
    "12-54-35": "listed as 12-54-25 in field notes; confirmed as position 14 run 3",
    "13-20-14": "missing from field notes; confirmed as position 10 run 9",
}


def original(ts):
    return DATA_DIR / f"magnetometer_2026-09-19_{ts}.csv"


def build_rows():
    rows = []

    def add(dike, folder, name, ts, experiment, position=10, run="",
            height_cm="", scan_time_s=""):
        rows.append({
            "labeled_file": f"{dike}/{folder}/{dike}-{name}_{ts}.csv",
            "original_file": original(ts).name,
            "dike": dike,
            "experiment": experiment,
            "position": position,
            "run": run,
            "height_cm": height_cm,
            "scan_time_s": scan_time_s,
            "note": NOTES.get(ts, ""),
        })

    # Numbers are zero-padded so an alphabetical listing (ls, GitHub) shows
    # them in numeric order; the folders are numbered for the same reason.
    for dike, positions in [("large-dike", LARGE_POSITIONS),
                            ("small-dike", SMALL_POSITIONS)]:
        for pos, stamps in sorted(positions.items()):
            for run, ts in enumerate(stamps, 1):
                add(dike, "1_positions-01-to-20", f"position-{pos:02d}-run-{run}",
                    ts, "positions", position=pos, run=run)

    # Position 10 already has runs 1-3 from the positions experiment, so the
    # repeats continue the count (13-20-14 is run 9, as in the field notes).
    first = len(LARGE_POSITIONS[10]) + 1
    for run, ts in enumerate(LARGE_REPEATS, first):
        add("large-dike", "2_position-10-runs", f"position-10-run-{run:02d}", ts,
            "repeats", run=run)
    for secs, stamps in sorted(LARGE_SPEEDS_S.items()):
        for run, ts in enumerate(stamps, 1):
            add("large-dike", "3_position-10-speeds",
                f"position-10-{secs:02d}sec-speed-run-{run}", ts, "speeds",
                run=run, scan_time_s=secs)
    for ts, h in sorted(LARGE_HEIGHTS_CM.items(), key=lambda kv: kv[1]):
        whole, _, frac = f"{h:g}".partition(".")
        label = f"{int(whole):03d}" + (f".{frac}" if frac else "")
        add("large-dike", "4_position-10-elevations",
            f"position-10-elevation-{label}cm", ts, "heights", height_cm=h)
    return rows


def copy_with_distance(src, dest):
    """Copy a run, inserting distance_cm (0 to RUN_LENGTH_CM) after time."""
    lines = src.read_text().splitlines()
    comments = [ln for ln in lines if ln.startswith("#")]
    header, *data = [ln for ln in lines if ln and not ln.startswith("#")]
    rows = [ln.split(",") for ln in data]
    t = [float(r[0]) for r in rows]
    span = t[-1] - t[0]

    cols = header.split(",")
    out = comments + [",".join([cols[0], "distance_cm", *cols[1:]])]
    for r, ti in zip(rows, t):
        d = RUN_LENGTH_CM * (ti - t[0]) / span
        out.append(",".join([r[0], f"{d:.2f}", *r[1:]]))
    dest.write_text("\n".join(out) + "\n")


def main():
    rows = build_rows()

    missing = [r["original_file"] for r in rows
               if not (DATA_DIR / r["original_file"]).exists()]
    if missing:
        raise SystemExit(f"Missing original files: {missing}")
    labeled = [r["labeled_file"] for r in rows]
    if len(set(labeled)) != len(labeled):
        raise SystemExit("Two runs map to the same labeled name.")
    originals = [r["original_file"] for r in rows]
    if len(set(originals)) != len(originals):
        raise SystemExit("One original file is labeled twice.")

    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)  # rebuilt from scratch so it always matches the mapping
    for r in rows:
        dest = OUT_DIR / r["labeled_file"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        copy_with_distance(DATA_DIR / r["original_file"], dest)

    with open(OUT_DIR / "metadata.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Copied {len(rows)} files into {OUT_DIR.relative_to(REPO)}/")


if __name__ == "__main__":
    main()
