#!/usr/bin/env python3
# ============================================================
#  Visualiser for the 2D grain-growth output (grain_growth_2d.cpp)
#
#  Reads:
#    <outdir>/phi_XXXXX.csv    -> phi = sum_i eta_i^2   (boundaries dip below 1)
#    <outdir>/grain_XXXXX.csv  -> argmax_i eta_i        (integer grain map)
#    <outdir>/eta_profiles.dat -> final mid-row line scan (the "X" crossing)
#
#  Produces (in <outdir>/):
#    grains.gif        coloured grain map, evolving in time
#    phi.gif           boundary-energy field, evolving in time
#    grains_first_last.png   initial vs final microstructure side by side
#    profiles.png      the eta_i(x) line scan across the mid row
#
#  Usage:   python plot_grain_growth.py fig2_out
#  Needs:   numpy, matplotlib   (gif export uses matplotlib's Pillow writer)
# ============================================================
import sys, os, glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")               # no display needed, just write files
import matplotlib.pyplot as plt
from matplotlib import animation


def load_series(outdir, prefix):
    """Load all <prefix>_XXXXX.csv frames, sorted by step number."""
    files = glob.glob(os.path.join(outdir, prefix + "_*.csv"))
    def step_of(fn):
        m = re.search(prefix + r"_(\d+)\.csv$", os.path.basename(fn))
        return int(m.group(1)) if m else -1
    files.sort(key=step_of)
    steps  = [step_of(f) for f in files]
    frames = [np.loadtxt(f) for f in files]
    return steps, frames


def make_gif(steps, frames, out_png, title, cmap, vmin=None, vmax=None, fps=8):
    if not frames:
        print(f"  (no frames for {title}, skipping)")
        return
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(frames[0], cmap=cmap, vmin=vmin, vmax=vmax,
                   origin="lower", interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    ttl = ax.set_title(f"{title}   step {steps[0]}")

    def update(f):
        im.set_data(frames[f])
        ttl.set_text(f"{title}   step {steps[f]}")
        return im, ttl

    anim = animation.FuncAnimation(fig, update, frames=len(frames),
                                   interval=1000 / fps, blit=False)
    anim.save(out_png, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)
    print(f"  wrote {out_png}  ({len(frames)} frames)")


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "fig2_out"
    if not os.path.isdir(outdir):
        sys.exit(f"no such directory: {outdir}  (run the C++ program first)")

    # ---- grain map: colour each orientation. tab20 gives distinct grain colours. ----
    gsteps, gframes = load_series(outdir, "grain")
    if gframes:
        nop = int(max(f.max() for f in gframes)) + 1
        make_gif(gsteps, gframes, os.path.join(outdir, "grains.gif"),
                 "grains", cmap=plt.get_cmap("tab20", max(nop, 2)),
                 vmin=0, vmax=nop - 1)

        # initial vs final side by side
        fig, axs = plt.subplots(1, 2, figsize=(9, 4.6))
        for ax, idx, lab in [(axs[0], 0, f"initial (Voronoi, step {gsteps[0]})"),
                             (axs[1], -1, f"final (step {gsteps[-1]})")]:
            ax.imshow(gframes[idx], cmap=plt.get_cmap("tab20", max(nop, 2)),
                      vmin=0, vmax=nop - 1, origin="lower", interpolation="nearest")
            ax.set_title(lab); ax.set_xticks([]); ax.set_yticks([])
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "grains_first_last.png"), dpi=140)
        plt.close(fig)
        print(f"  wrote {os.path.join(outdir, 'grains_first_last.png')}")

    # ---- phi field: grain interiors ~1, boundaries darker. ----
    psteps, pframes = load_series(outdir, "phi")
    make_gif(psteps, pframes, os.path.join(outdir, "phi.gif"),
             "phi = sum eta_i^2", cmap="magma", vmin=0.0, vmax=1.05)

    # ---- the "X" line scan: eta_i(x) across the mid row ----
    prof = os.path.join(outdir, "eta_profiles.dat")
    if os.path.exists(prof):
        data = np.loadtxt(prof)
        x = data[:, 0]
        fig, ax = plt.subplots(figsize=(8, 4))
        for p in range(1, data.shape[1]):
            col = data[:, p]
            if col.max() > 0.5:                 # only draw orientations present on this row
                ax.plot(x, col, lw=1.3, label=f"eta_{p-1}")
        ax.set_xlabel("x"); ax.set_ylabel("eta_i")
        ax.set_title("mid-row line scan (each crossing = a grain boundary)")
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "profiles.png"), dpi=140)
        plt.close(fig)
        print(f"  wrote {os.path.join(outdir, 'profiles.png')}")

    print("done.")


if __name__ == "__main__":
    main()