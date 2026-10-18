"""Copy the big-dike runs into labeled_data/ with descriptive names.

The originals in magnetic_field_data/ are never modified, and
labeled_data/metadata.csv maps every labeled file back to its original.

Layout (folders and names sort in this order):
  1_positions-01-to-20/      large-dike-position-01-run-1 ... position-20-run-3
  2_position-10-runs/        large-dike-position-10-run-04 ... run-15
  3_position-10-speeds/      large-dike-position-10-03sec-speed-run-1 ... 15sec-speed-run-3
  4_position-10-elevations/  large-dike-position-10-elevation-001.1cm ... 100cm

Experiments (from the field notes):
  1. Locations 1-20 along the dike, 10 cm apart, 3 identical passes each
  2. Location 10 at different phone heights above the rock
  3. Location 10 repeated under the same conditions
  4. Location 10 at different scan speeds (seconds to cross the 90 cm)

Usage (from the repo root):
    python analysis/label_files.py
"""
import csv
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "magnetic_field_data"
OUT_DIR = REPO / "labeled_data"
PREFIX = "large-dike-"  # every run in this mapping is from the large dike

# 1. Location -> three repeat passes
LOCATIONS = {
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

# 2. Location 10, phone height above the rock (cm)
HEIGHTS_CM = {
    "13-07-13": 1.1, "13-09-14": 2, "13-10-40": 3, "13-11-17": 6,
    "13-11-45": 9, "13-12-09": 12, "13-12-53": 15, "13-13-36": 30,
    "13-14-17": 50, "13-15-16": 70, "13-16-16": 100,
}

# 3. Location 10, repeated under the same conditions
REPEATS = ["13-18-20", "13-18-44", "13-19-10", "13-19-31", "13-19-53",
           "13-20-14", "13-20-41", "13-21-04", "13-21-24", "13-21-45",
           "13-22-06", "13-22-27"]

# 4. Location 10, seconds taken to cross the 90 cm
SPEEDS_S = {
    3: ["13-26-03", "13-26-18", "13-26-39"],
    5: ["13-23-05", "13-23-29", "13-23-49"],
    15: ["13-24-35", "13-25-03", "13-25-41"],
}

NOTES = {
    "12-54-35": "listed as 12-54-25 in field notes; confirmed as position 14 run 3",
    "13-20-14": "missing from field notes; confirmed as position 10 run 9",
}


def original(ts):
    return DATA_DIR / f"magnetometer_2026-09-19_{ts}.csv"


def build_rows():
    rows = []

    def add(folder, name, ts, experiment, location=10, run="",
            height_cm="", scan_time_s=""):
        rows.append({
            "labeled_file": f"{folder}/{PREFIX}{name}.csv",
            "original_file": original(ts).name,
            "experiment": experiment,
            "location": location,
            "run": run,
            "height_cm": height_cm,
            "scan_time_s": scan_time_s,
            "note": NOTES.get(ts, ""),
        })

    # Numbers are zero-padded so an alphabetical listing (ls, GitHub) shows
    # them in numeric order; the folders are numbered for the same reason.
    for loc, stamps in sorted(LOCATIONS.items()):
        for run, ts in enumerate(stamps, 1):
            add("1_positions-01-to-20", f"position-{loc:02d}-run-{run}", ts,
                "locations", location=loc, run=run)
    # Position 10 already has runs 1-3 from the locations experiment, so the
    # repeats continue the count (13-20-14 is run 9, as in the field notes).
    first = len(LOCATIONS[10]) + 1
    for run, ts in enumerate(REPEATS, first):
        add("2_position-10-runs", f"position-10-run-{run:02d}", ts, "repeats",
            run=run)
    for secs, stamps in sorted(SPEEDS_S.items()):
        for run, ts in enumerate(stamps, 1):
            add("3_position-10-speeds",
                f"position-10-{secs:02d}sec-speed-run-{run}", ts, "speeds",
                run=run, scan_time_s=secs)
    for ts, h in sorted(HEIGHTS_CM.items(), key=lambda kv: kv[1]):
        whole, _, frac = f"{h:g}".partition(".")
        label = f"{int(whole):03d}" + (f".{frac}" if frac else "")
        add("4_position-10-elevations", f"position-10-elevation-{label}cm", ts,
            "heights", height_cm=h)
    return rows


def main():
    rows = build_rows()

    missing = [r["original_file"] for r in rows
               if not (DATA_DIR / r["original_file"]).exists()]
    if missing:
        raise SystemExit(f"Missing original files: {missing}")
    labeled = [r["labeled_file"] for r in rows]
    if len(set(labeled)) != len(labeled):
        raise SystemExit("Two runs map to the same labeled name.")

    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)  # rebuilt from scratch so it always matches the mapping
    for r in rows:
        dest = OUT_DIR / r["labeled_file"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(DATA_DIR / r["original_file"], dest)

    with open(OUT_DIR / "metadata.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Copied {len(rows)} files into {OUT_DIR.relative_to(REPO)}/")


if __name__ == "__main__":
    main()
