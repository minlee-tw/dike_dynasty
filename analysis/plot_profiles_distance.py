"""Plot each magnetometer run as field vs distance along the 90 cm profile.

Time is converted to distance assuming a constant scan speed over the whole
recording: x = L * (t - t_first) / (t_last - t_first). A 10 s run over 0.9 m
is therefore 9 cm per second.

--smooth W applies an overlapping (sliding) moving average: every point is
replaced by the mean of all samples within ±W/2 metres of it, so neighbouring
windows overlap almost completely. Near the ends of a run the window is cut
off by the edge, so it averages fewer samples there.

Usage (from the repo root):
    python analysis/plot_profiles_distance.py                  # Btotal, first file through 12-02-52
    python analysis/plot_profiles_distance.py --component Bz
    python analysis/plot_profiles_distance.py --overlay        # all runs on one graph
    python analysis/plot_profiles_distance.py --overlay --smooth 0.05
    python analysis/plot_profiles_distance.py --start 12-36-58 --end 13-26-39
"""
import argparse
import math

import matplotlib.pyplot as plt
import numpy as np

from plot_btotal import FIG_DIR, REPO, load, select_files, stamp


def to_distance(t, length):
    return length * (t - t.iloc[0]) / (t.iloc[-1] - t.iloc[0])


def moving_average(x, y, width):
    """Centred sliding mean of y over a window `width` wide in x (x sorted)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    csum = np.concatenate([[0.0], np.cumsum(y)])
    lo = np.searchsorted(x, x - width / 2, side="left")
    hi = np.searchsorted(x, x + width / 2, side="right")
    return (csum[hi] - csum[lo]) / (hi - lo)


def smooth_tag(args):
    return f"_smooth{args.smooth * 100:g}cm" if args.smooth else ""


def title(args, files):
    first, last = stamp(files[0]), stamp(files[-1])
    t = (f"{args.component} vs distance, 2026-09-19, {first} to {last} "
         f"({len(files)} files, {args.length:g} m runs)")
    if args.smooth:
        t += f"\n{args.smooth * 100:g} cm overlapping moving average"
    return t


def plot_overlay(files, args):
    colors = plt.cm.turbo(np.linspace(0, 1, len(files)))
    fig, ax = plt.subplots(figsize=(16, 10), layout="constrained")

    for color, f in zip(colors, files):
        df = load(f)
        x = to_distance(df["time"], args.length)
        y = df[args.component]
        if args.smooth:
            y = moving_average(x, y, args.smooth)
        ax.plot(x, y, lw=0.6 if args.smooth else 0.4, color=color, label=stamp(f))

    first, last = stamp(files[0]), stamp(files[-1])
    ax.set_xlim(0, args.length)
    ax.set_xlabel("Distance along profile (m)")
    ax.set_ylabel(f"{args.component} (µT)")
    ax.set_title(title(args, files))
    ax.grid(alpha=0.3)
    leg = ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), ncols=2,
                    fontsize=9, title="Run (start time)", frameon=False)
    for line in leg.get_lines():
        line.set_linewidth(3)

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / (f"{args.component.lower()}_distance_overlay{smooth_tag(args)}"
                     f"_{first}_to_{last}.png")
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(REPO)} ({len(files)} files)")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--start", default=None, help="first timestamp, HH-MM-SS")
    p.add_argument("--end", default="12-02-52", help="last timestamp, HH-MM-SS")
    p.add_argument("--length", type=float, default=0.9, help="profile length (m)")
    p.add_argument("--component", default="Btotal",
                   choices=["Btotal", "Bx", "By", "Bz"])
    p.add_argument("--ncols", type=int, default=6)
    p.add_argument("--overlay", action="store_true",
                   help="draw every run on one set of axes with a legend")
    p.add_argument("--smooth", type=float, default=None,
                   help="moving-average window width in metres, e.g. 0.05")
    args = p.parse_args()

    files = select_files(args.start, args.end)
    if not files:
        raise SystemExit("No files in that time range.")
    if args.overlay:
        plot_overlay(files, args)
        return

    ncols = args.ncols
    nrows = math.ceil(len(files) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.2 * nrows),
                             squeeze=False, layout="constrained")

    for ax, f in zip(axes.flat, files):
        df = load(f)
        x = to_distance(df["time"], args.length)
        y = df[args.component]
        if args.smooth:
            ax.plot(x, y, lw=0.5, color="0.75")  # raw data underneath
            ax.plot(x, moving_average(x, y, args.smooth), lw=1.0, color="tab:blue")
        else:
            ax.plot(x, y, lw=0.8, color="tab:blue")
        ax.set_title(f"{stamp(f)}  ({df['time'].iloc[-1]:.1f} s)", fontsize=9)
        ax.set_xlim(0, args.length)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    for ax in axes.flat[len(files):]:
        ax.axis("off")

    fig.supxlabel("Distance along profile (m)")
    fig.supylabel(f"{args.component} (µT)")
    first, last = stamp(files[0]), stamp(files[-1])
    fig.suptitle(title(args, files))

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / (f"{args.component.lower()}_distance{smooth_tag(args)}"
                     f"_{first}_to_{last}.png")
    fig.savefig(out, dpi=130)
    print(f"Saved {out.relative_to(REPO)} ({len(files)} files)")


if __name__ == "__main__":
    main()
