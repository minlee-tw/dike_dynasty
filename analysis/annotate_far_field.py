"""Write the small-dike far-field results into named_data/.

Reads figures/small-dike_far_field_cutoffs.csv (from small_dike_far_field.py)
and adds:
  - to every small-dike run file, a `region` column: `left` or `right` for
    far-field samples on either side of the dike, `dike` for samples between
    the cutoffs. A run with no dike found is far field throughout, split at
    the dike-template match into `left` and `right`.
  - to named_data/metadata.csv, per-run columns for the small-dike runs:
    dike_found, far_left_below_cm, far_right_above_cm (cutoffs in cm along the
    run, as in distance_cm), far_mean_uT, far_sd_uT, far_se_uT and kept.
    Large-dike rows are left blank.
`kept` is "no" for the five runs dropped for offset baselines (EXCLUDED: far-
field value outside 2σ of the run average, P18 r2 once the other four are
removed). Running it again replaces the columns it wrote.

label_files.py rebuilds named_data/ from scratch, which removes these columns,
so run this again after it.

Usage (from the repo root):
    python analysis/label_files.py              # only if the mapping changed
    python analysis/small_dike_quality.py
    python analysis/small_dike_far_field.py
    python analysis/annotate_far_field.py
"""
import csv

from plot_btotal import FIG_DIR, REPO
from plot_positions_overlay import NAMED_DIR

EXCLUDED = {(1, 1), (18, 1), (18, 2), (15, 2), (20, 3)}  # (position, run)
META_COLUMNS = ["dike_found", "far_left_below_cm", "far_right_above_cm", "far_mean_uT",
                "far_sd_uT", "far_se_uT", "kept"]


def label_run(path, cut):
    """Rewrite one run file with a `region` column (replacing an old one)."""
    lines = path.read_text().splitlines()
    out, d_col, drop = [], None, None
    for ln in lines:
        if not ln or ln.startswith("#"):
            out.append(ln)
            continue
        cells = ln.split(",")
        if d_col is None:  # header
            drop = cells.index("region") if "region" in cells else None
            if drop is not None:
                del cells[drop]
            d_col = cells.index("distance_cm")
            out.append(",".join(cells + ["region"]))
            continue
        if drop is not None:
            del cells[drop]
        d = float(cells[d_col])
        if cut["dike_found"] == "yes":
            a, b = float(cut["far_left_below_cm"]), float(cut["far_right_above_cm"])
            region = "left" if d < a else "right" if d > b else "dike"
        else:
            region = "left" if d < float(cut["match_cm"]) else "right"
        out.append(",".join(cells + [region]))
    path.write_text("\n".join(out) + "\n")


def main():
    with open(FIG_DIR / "small-dike_far_field_cutoffs.csv", newline="") as fh:
        cuts = {(int(r["position"]), int(r["run"])): r for r in csv.DictReader(fh)}

    meta_path = NAMED_DIR / "metadata.csv"
    with open(meta_path, newline="") as fh:
        reader = csv.DictReader(fh)
        fields = [f for f in reader.fieldnames if f not in META_COLUMNS] + META_COLUMNS
        rows = list(reader)

    labelled = 0
    for row in rows:
        for col in META_COLUMNS:
            row[col] = ""
        if row["dike"] != "small-dike" or row["experiment"] != "positions":
            continue
        key = (int(row["position"]), int(row["run"]))
        cut = cuts[key]
        label_run(NAMED_DIR / row["labeled_file"], cut)
        labelled += 1
        row.update({col: cut[col] for col in META_COLUMNS[:-1]})
        row["kept"] = "no" if key in EXCLUDED else "yes"

    with open(meta_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Added region to {labelled} small-dike run files and far-field columns to "
          f"{meta_path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
