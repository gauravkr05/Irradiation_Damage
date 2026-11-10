#!/usr/bin/env python3
# Shrinking of a circular grain, coloured by orientation.
#   blue  = matrix, orientation 1 (eta_1 = 1)
#   red   = circular grain, orientation 2 (eta_2 = 1)
# reads f5_step_<n>/eta_step_<n>.dat written by fig5_circle.cpp
# bare images: no axes, no ticks, no labels, one pixel per lattice point
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

STEPS = [200, 2000, 4000, 4600, 4800, 5000]

# blue -> pale boundary -> red
cmap = LinearSegmentedColormap.from_list(
    "orient", ["#1b4f9c", "#6f8fc4", "#efe8e0", "#cf7a6b", "#b32a1e"])

for s in STEPS:
    raw = np.loadtxt(f"f5_step_{s}/eta_step_{s}.dat")
    nx  = int(raw[:, 0].max()) + 1
    ny  = int(raw[:, 1].max()) + 1
    e1  = raw[:, 2].reshape(nx, ny)
    e2  = raw[:, 3].reshape(nx, ny)

    # orientation indicator: 0 where orientation 1 wins, 1 where orientation 2 wins
    c = e2**2 / (e1**2 + e2**2 + 1e-12)

    plt.imsave(f"fig5_step_{s}.png", c.T, cmap=cmap,
               origin="lower", vmin=0.0, vmax=1.0)
    print(f"wrote fig5_step_{s}.png   grain fraction = {(c > 0.5).mean():.4f}")
