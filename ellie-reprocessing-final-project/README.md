# Ellie's Bz reprocessing: start here

This folder holds Ellie's part of the Dike Dynasty final project: both dikes reprocessed **using Bz only**. It is self-contained. The scripts, the processed data and the method all live here, and nothing here depends on the group's earlier Btotal analysis.

## The short version

- **Small dike:** each run's dike is found from Bz (LOESS slope plus a template check), the far field is cut out, and a single far-field Bz is computed.
  **Result: −42.15 ± 0.23 µT** (standard error, 52 of 60 runs). The dike is confirmed in **47 of 60** runs.
- **Large dike:** the phone's recalibration jumps are stitched out of Bz. Each run is then shifted so its far end (80–90 cm) sits at the small-dike far-field value, −42.15 µT. No separate far field was worked out for the large dike.
- **Every sample of every run** is saved as a CSV, so you can use the processed data without rerunning anything.

## Where everything is

```
ellie-reprocessing-final-project/
├── README.md                         ← this file
├── common.py                         shared helpers (loading runs, LOESS slope, robust SD)
├── small_dike_bz.py                  small dike: find the dike, far-field cutoffs, far-field Bz
├── large_dike_bz.py                  large dike: remove recalibration steps, anchor to −42.15 µT
├── plot_small_dike_bz_raw.py         plain plots of the small-dike Bz as recorded
├── data/
│   ├── small-dike/
│   │   ├── runs/                     60 files, one per run: every sample + its region
│   │   ├── small-dike_bz_runs_summary.csv      one row per run
│   │   └── small-dike_bz_far_field_value.csv   the final far-field number
│   └── large-dike/
│       ├── runs/                     60 files, one per run: raw, corrected and anchored Bz
│       └── large-dike_bz_runs_summary.csv      one row per run
└── figures/                          only created with --figures; not tracked by git
```

The original measurements are in `named_data/` at the top of the repo. They are **only read**, never changed. Run names match: `named_data/small-dike/1_positions-01-to-20/small-dike-position-05-run-1_11-17-28.csv` becomes `data/small-dike/runs/small-dike-position-05-run-1_11-17-28_bz.csv`.

## Using the data

Each per-run file starts with one `#` line summarising the run (cutoffs, or the shift applied), so skip it when reading:

```python
import pandas as pd
run = pd.read_csv("ellie-reprocessing-final-project/data/small-dike/runs/"
                  "small-dike-position-05-run-1_11-17-28_bz.csv", comment="#")
far_field = run[(run.region != "dike") & (run.repeated_sample == "no")]
```

### Small-dike run files (`data/small-dike/runs/`)

| Column | Meaning |
| --- | --- |
| `time`, `distance_cm`, `Bz` | as recorded (µT) |
| `repeated_sample` | `yes` where the logger wrote the previous sample again; left out of all calculations |
| `distance_from_template_centre_cm` | distance from where the dike template matched (its steepest drop) |
| `loess_slope_uT_per_cm` | LOESS slope of Bz (6 cm span); blank near the ends and on repeated samples |
| `steep` | `yes` where the slope is above the dike threshold (0.19 µT/cm) |
| `region` | `far_field_left`, `dike` or `far_field_right`. Runs with no confirmed dike are all far field, split at the template centre |

### Small-dike summary (`small-dike_bz_runs_summary.csv`)

There is one row per run, with these columns:

- **Template match:** `match` (strong, weak or none), `correlation`, `size_uT` (dip size) and `centre_cm`.
- **Dike:** `dike_confirmed`, plus the cutoffs `far_left_below_cm` and `far_right_above_cm`. Far field is below the first and above the second, in cm along the run.
- **Far field:** `far_mean_uT`, `far_sd_uT`, `far_se_uT` and `far_samples`.
- **`kept`:** `no` for the 8 outlier runs left out of the final value.

`small-dike_bz_far_field_value.csv` has the final number and how it was obtained.

### Large-dike run files (`data/large-dike/runs/`)

| Column | Meaning |
| --- | --- |
| `time`, `distance_cm`, `Bz_raw` | as recorded (µT) |
| `recalibration_step` | `yes` on the first sample after each removed recalibration jump |
| `Bz_step_corrected` | jumps removed; the level is arbitrary |
| `Bz_anchored` | **use this one**: corrected and shifted so the 80–90 cm median is −42.15 µT |
| `anchor_window` | `yes` for the 80–90 cm samples used to set the level |

`large-dike_bz_runs_summary.csv` gives, per run: the number of steps removed, the shift applied, the anchored range and `far_end_loess_slope_uT_per_cm`. That last column shows how flat the anchor stretch is; larger values make the level less certain.

## Rerunning everything

From the repo root, in this order (the large dike reads the small-dike result):

```bash
python ellie-reprocessing-final-project/small_dike_bz.py
```

```bash
python ellie-reprocessing-final-project/large_dike_bz.py
```

Add `--figures` to either one to also make the figures in `figures/`. Each takes under a minute.

## What we did, and why

### Small dike

1. **Slope (LOESS).** At each sample, a straight line is fitted to the samples within ±3 cm, weighting nearby samples most (tricube weights; local-linear LOESS, 6 cm span). Its slope is taken. LOESS was preferred to a moving average, and the 6 cm span has the same noise level as the 4 cm least-squares window the group used before.
2. **Steep stretches.** |slope| > **2.5σ** for at least **0.5 cm**. σ is the robust SD of the slope at least 30 cm from the dike, over all 60 runs: σ = 0.076 µT/cm, so the threshold is 0.19 µT/cm.
3. **Template.** The 29 runs whose Bz varies by more than 3 µT are lined up on their dip, refined by cross-correlation, and averaged. That gives the dike's typical Bz shape: a small rise, a sharp drop, a dip, and a slow recovery over about 30 cm. Its centre is the steepest point of the drop.
4. **Match.** The template slides along each run, and the best correlation marks where the dike is. Each match is rated:
   - **strong:** correlation ≥ 0.7 and dip ≥ 2 µT (23 runs)
   - **weak:** correlation ≥ 0.5 and dip ≥ 0.5 µT (34 runs)
   - **none:** anything less (3 runs)
5. **Confirmation: the slope finds the dike, the template confirms it.** A steep stretch counts only if it overlaps the template's footprint, −4 to +12.5 cm from its centre. The dike then runs from those stretches through the template centre. Nearby stretches less than 5 cm apart are joined, but only if they start within ±15 cm of the centre, which keeps unrelated dips near the run ends out. A **3 cm margin** is added on each side.
   - **One manual call:** P18 r2 is set to not confirmed (weak match, and its only steep stretch sits at the edge of the footprint).
   - **Not confirmed (13 runs):** P02 r1–r3, P11 r3, P14 r1, P15 r1–r3, P18 r2, P19 r2, P20 r1–r3. Positions 2 and 20 have no confirmed run.
6. **Far-field value per run.** The mean of every sample outside the dike, or of the whole run if no dike was confirmed. The SE comes from 4 cm block means, because neighbouring samples are correlated. Cutting a template-width window out of the unconfirmed runs instead changes their values by at most 0.2 µT.
7. **Outlier runs.** Runs more than **1 SD** from the mean of all 60 runs are excluded, in one pass: P01 r2, P01 r3, P18 r1, P18 r2, P18 r3, P19 r1, P19 r3 and P20 r3. Most are runs where the phone was tilted.

**Result:** far-field Bz = **−42.15 ± 0.23 µT** (SE, n = 52). Median −42.19 µT; SD across runs 1.64 µT.

For comparison, the NOAA World Magnetic Model (WMM2025) puts the vertical field at the site on 2026-09-19 at 46.86 µT downward, which a level, screen-up phone reads as **−46.86 µT**. Tilting the phone can only make Bz less negative, so −42.15 µT means "Bz as the phone was typically held", not the true vertical field. It fits an average tilt of about 25°.

### Large dike

1. **Recalibration steps.** The phone recalibrated about every 0.41 s without pausing the 100 Hz recording, so Bz jumps by hundreds to thousands of µT. A step is flagged when Bz changes by more than 20 µT from one sample to the next and by more than 10× the median of the previous 5 changes, plus 5 µT. Bz after the step is shifted to continue the straight-line trend of the 5 samples before it (one-sided linear extrapolation). There is no smoothing. Between 8 and 25 steps are removed per run.
   - **Testing:** on small-dike Bz with fake steps added, this recovered the true profile to about 1.5 µT RMS. Fancier methods weren't meaningfully better, and on the real data all methods agreed to within 1 µT.
2. **Anchoring.** Stitching keeps each run's shape but not its level, and the raw level before the first recalibration is unreliable (−171 to +305 µT). So each run is shifted so that the **median Bz over 80–90 cm equals −42.15 µT**, the small-dike far field.
   - **Why the far end:** the dike's steepest drop is at a median of 31 cm along the run, so the far end is about twice as far from it as the start. The far end is also the flatter end in 54 of 60 runs (about 1.6 against 5 µT/cm).
   - **Alternatives:** anchoring the start, or the mean of both ends, gave the same agreement between repeat runs (within about 1 µT) but would put the far, flattest end 70–145 µT off background.
   - **Why the median over 10 cm:** it uses about 100 samples and ignores spikes. Using 85–90 cm or a mean instead changes the anchor by about 1 µT typically.

## Things to keep in mind

- **Bz only.** No Btotal, and no labels or results from the group's Btotal analysis, are used anywhere here. Please don't mix the two when using these files.
- **Small-dike far field** depends on phone tilt. Repeat runs at one position can differ by a few µT. The SE covers random scatter between runs, not the tilt bias.
- **Weak runs.** At positions 1, 2 and 10–20 the dike is only about 0.5–1 µT in Bz, close to the noise. Keeping only confirmed runs favours the larger-looking dikes there, so report how many runs passed at each position.
- **The template centre** marks the steepest drop, not necessarily the dike's geometric centre. In a magnetization model fit, leave the dike centre as a free parameter.
- **Large-dike levels** are relative to the far end, taken as background. Some dike field remains at 80–90 cm, so absolute levels may be off by tens of µT where that end is still steep (P02, P03, P08, P09). Shapes and anomaly sizes aren't affected.
- **Small leftover kinks:** large-dike recalibrations that moved Bz by less than 20 µT weren't caught.
- **Large-dike runs can't be rejected for an odd far-field level,** because every run's far end is set to the same value.

## Relation to the rest of the repo

- `named_data/`: the original labelled runs. Read only.
- `analysis/`, `figures/`, `reports/`: the group's earlier work, mostly Btotal-based. It isn't used here.
- `corrections/large-dike/`: an earlier copy of the large-dike step-removal figures, from before this folder existed. The current script is `large_dike_bz.py` here.
