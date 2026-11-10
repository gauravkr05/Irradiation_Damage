"""Panels of a single bubble growing in a supersaturated matrix, plus radius vs time.

Run the solver first:
    g++ -O3 -fopenmp -std=c++17 single_bubble.cpp -o sb
    ./sb 160 160 240000 24000
    python3 plot_single_bubble.py
"""
import glob
import numpy as np
import matplotlib.pyplot as plt

NPANELS = 3

meta = {}
for line in open("meta.txt"):
    k, v = line.split()
    meta[k] = float(v)
tscale = meta["tscale_s_per_tsim"]     # seconds of real time per unit of sim time
dt = meta["dt"]


def read_snap(path):
    with open(path, "rb") as fh:
        nx, ny = np.fromfile(fh, dtype=np.int32, count=2)
        return np.fromfile(fh, dtype=np.float32, count=nx * ny).reshape(ny, nx)


files = sorted(glob.glob("snap_*.bin"))
pick = [files[i] for i in np.linspace(0, len(files) - 1, NPANELS).astype(int)]

fig, axs = plt.subplots(1, NPANELS, figsize=(4.0 * NPANELS, 4.4))
for n, (ax, f) in enumerate(zip(np.atleast_1d(axs), pick)):
    c = read_snap(f)
    im = ax.imshow(c, cmap="RdBu_r", vmin=0.0, vmax=1.0, interpolation="bilinear")
    step = int(f.split("_")[1].split(".")[0])
    hours = step * dt * tscale / 3600.0
    ax.set_title(f"({chr(97 + n)})  t = {hours:.2f} h", fontsize=12, loc="left")
    ax.set_xticks([]); ax.set_yticks([])

cb = fig.colorbar(im, ax=np.atleast_1d(axs).tolist(), fraction=0.025, pad=0.02)
cb.set_label("c")
plt.savefig("single_bubble_panels.png", dpi=200, bbox_inches="tight")

# radius against time
d = np.loadtxt("bubble_radius.dat", comments="#")
fig2, ax2 = plt.subplots(figsize=(5.6, 4.4))
ax2.plot(d[:, 2], d[:, 4] * meta["dx_nm"], lw=1.8, color="crimson")
ax2.set_xlabel("time (hrs)")
ax2.set_ylabel("bubble radius (nm)")
ax2.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("single_bubble_radius.png", dpi=200)
plt.show()
