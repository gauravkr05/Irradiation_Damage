"""
plot_morphology.py  -  top views (bubble thickness along z, i.e. what the
coverage metric sees) and x-z cross-sections through the GB, at identical
simulation times for all cases. Snapshot times are read from snapshots.csv.

    python3 plot_morphology.py [root=cases] [times=0,9,21,39]
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NX, NY, NZ = 128, 128, 64
CASES = ["1.2", "1.5", "2.0", "3.0"]
args = dict(a.split("=") for a in sys.argv[1:])
root = args.get("root", "cases")
want = [float(x) for x in args.get("times", "0,9,21,39").split(",")]
props = json.load(open("results/case_properties.json"))


def load(b, t):
    s = pd.read_csv(os.path.join(root, f"bGB_{b}", "snapshots.csv"))
    r = s.iloc[(s.time_dimensionless - t).abs().argmin()]
    if abs(r.time_dimensionless - t) > 1e-6:
        return None
    return np.fromfile(os.path.join(root, f"bGB_{b}", r.c_file), dtype=np.float32).reshape(NZ, NY, NX)


for kind in ("top", "side"):
    fig, axs = plt.subplots(len(CASES), len(want), figsize=(3.2 * len(want), 3.0 * len(CASES) if kind == "top" else 1.9 * len(CASES)))
    for i, b in enumerate(CASES):
        for k, t in enumerate(want):
            ax = axs[i, k]
            c = load(b, t)
            ax.set_xticks([]); ax.set_yticks([])
            if c is None:
                ax.text(0.5, 0.5, "n/a", ha="center"); continue
            if kind == "top":
                thick = (c >= 0.5).sum(axis=0)
                ax.imshow(thick, origin="lower", cmap="magma", vmin=0, vmax=30)
                ax.contour(thick >= 3, levels=[0.5], colors="c", linewidths=0.5)
            else:
                ax.imshow(c[:, NY // 2, :], origin="lower", cmap="viridis", vmin=0, vmax=1, aspect="equal")
                ax.axhline(NZ / 2 - 0.5, color="r", lw=0.4, ls="--")
            if i == 0:
                ax.set_title(f"t = {t:g}")
            if k == 0:
                ax.set_ylabel(f"bGB={b}\nθ≈{props[b]['theta_meas']:.0f}°", fontsize=9)
    ttl = ("Top view: number of cells with c ≥ 0.5 in each column (cyan: covered, ≥3 cells)"
           if kind == "top" else "x-z section at y = 64 (red dashed: GB plane)")
    fig.suptitle(ttl, fontsize=10)
    fig.tight_layout()
    fig.savefig(f"results/morphology_{kind}view.png", dpi=150)
print("done")
