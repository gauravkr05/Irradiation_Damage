#!/usr/bin/env python3
# All six snapshots of the shrinking circular grain in ONE image.
#   blue = matrix (orientation 1),  red = circular grain (orientation 2)
# reads f5_step_<n>/eta_step_<n>.dat written by fig5_circle.cpp
# panels have no ticks or frames; only the step number is written under each
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

STEPS = [200, 2000, 4000, 4600, 4800, 5000]
NROW, NCOL = 2, 3
DT = 0.1

cmap = LinearSegmentedColormap.from_list(
    "orient", ["#1b4f9c", "#6f8fc4", "#efe8e0", "#cf7a6b", "#b32a1e"])

print("working directory :", os.getcwd())
missing = [s for s in STEPS if not os.path.isfile(f"f5_step_{s}/eta_step_{s}.dat")]
if missing:
    print("missing data for steps:", missing)
    print("run  ./fig5_circle  in THIS directory first")
    raise SystemExit(1)

fig, axes = plt.subplots(NROW, NCOL, figsize=(3 * NCOL, 3.15 * NROW))

for ax, s in zip(axes.ravel(), STEPS):
    raw = np.loadtxt(f"f5_step_{s}/eta_step_{s}.dat")
    nx  = int(raw[:, 0].max()) + 1
    ny  = int(raw[:, 1].max()) + 1
    e1  = raw[:, 2].reshape(nx, ny)
    e2  = raw[:, 3].reshape(nx, ny)

    c = e2**2 / (e1**2 + e2**2 + 1e-12)      # 0 = matrix, 1 = grain

    ax.imshow(c.T, cmap=cmap, origin="lower", vmin=0.0, vmax=1.0,
              interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f"step {s}   (t = {s * DT:g})", fontsize=10, pad=6)

    print(f"step {s:5d} : grain area fraction = {(c > 0.5).mean():.4f}")

fig.subplots_adjust(wspace=0.06, hspace=0.16)
fig.savefig("fig5_evolution.png", dpi=200, bbox_inches="tight")
print("wrote", os.path.abspath("fig5_evolution.png"))