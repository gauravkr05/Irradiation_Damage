#!/usr/bin/env python3
# ============================================================================
#  Plotter for fig7.cpp  (Fig. 7 of Jiang et al. 2022)
#  Reads fig7_out/{h_*.csv, sg_*.csv} and composites the paper's colouring:
#     red = bubble (h high), white = grain boundary (grain sum low),
#     blue = grain interior. Writes fig7_out/fig5poly_evolution.png (+ bubbles-only).
#  Usage:  python fig7_plot.py
# ============================================================================
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

OUTDIR = "fig5poly_out"
PANELS = [("F20", "2.0e21 f/cm$^3$"), ("F35", "3.5e21 f/cm$^3$"),
          ("F40", "4.0e21 f/cm$^3$"), ("F45", "4.5e21 f/cm$^3$")]

BLUE, WHITE, RED = (np.array(to_rgb("#1f3b8c")),
                    np.array(to_rgb("#f5f5f5")),
                    np.array(to_rgb("#c62828")))

def smoothstep(x, a, b):
    t = np.clip((x - a) / (b - a), 0.0, 1.0)
    return t * t * (3 - 2 * t)

def composite(h, sg):
    wb = smoothstep(h, 0.4, 0.7)
    notb = 1.0 - wb
    gbness = smoothstep(0.6 - sg, 0.0, 0.5)
    gb = notb * gbness
    grain = notb * (1.0 - gbness)
    return np.clip(grain[..., None]*BLUE + gb[..., None]*WHITE + wb[..., None]*RED, 0, 1)

def main():
    if not os.path.isdir(OUTDIR):
        raise SystemExit(f"{OUTDIR} not found. Run ./fig5poly first.")

    fig, axs = plt.subplots(2, 2, figsize=(8.6, 8.8)); axs = axs.ravel()
    for ax, (tag, label), let in zip(axs, PANELS, "abcd"):
        hf, sf = f"{OUTDIR}/h_{tag}.csv", f"{OUTDIR}/sg_{tag}.csv"
        if not (os.path.exists(hf) and os.path.exists(sf)):
            ax.set_title(f"({let}) {label}\n[missing]"); ax.set_xticks([]); ax.set_yticks([]); continue
        ax.imshow(composite(np.loadtxt(hf), np.loadtxt(sf)), origin="lower", interpolation="nearest")
        ax.set_title(f"({let}) {label}", fontsize=11); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Fig. 5 (polycrystalline): recrystallization + bubbles vs fission density", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f"{OUTDIR}/fig5poly_evolution.png", dpi=150); plt.close(fig)
    print("wrote", f"{OUTDIR}/fig5poly_evolution.png")

    fig, axs = plt.subplots(2, 2, figsize=(8.6, 8.8)); axs = axs.ravel()
    for ax, (tag, label), let in zip(axs, PANELS, "abcd"):
        hf = f"{OUTDIR}/h_{tag}.csv"
        if not os.path.exists(hf): ax.set_xticks([]); ax.set_yticks([]); continue
        ax.imshow(np.loadtxt(hf), origin="lower", cmap="inferno", vmin=0, vmax=1, interpolation="nearest")
        ax.set_title(f"({let}) {label}", fontsize=11); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Fig. 5 (polycrystalline): bubble field h", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f"{OUTDIR}/fig5poly_bubbles_only.png", dpi=150); plt.close(fig)
    print("wrote", f"{OUTDIR}/fig5poly_bubbles_only.png")

if __name__ == "__main__":
    main()
