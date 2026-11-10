"""Render the bubble morphology on the boundary plane from the snapshots.
Run the solver with DUMP set, e.g.  DUMP=50000 ./pf2 256 256 600000 500
"""
import glob
import numpy as np
import matplotlib.pyplot as plt

files = sorted(glob.glob("snap_*.bin"))
if len(files) > 6:                      # keep six evenly spaced frames
    files = [files[i] for i in np.linspace(0, len(files) - 1, 6).astype(int)]

fig, axs = plt.subplots(1, len(files), figsize=(3.0 * len(files), 3.2))
for ax, f in zip(np.atleast_1d(axs), files):
    with open(f, "rb") as fh:
        nx, ny = np.fromfile(fh, dtype=np.int32, count=2)
        c = np.fromfile(fh, dtype=np.float32, count=nx * ny).reshape(ny, nx)
    bubble = c > 0.5
    ax.imshow(bubble, cmap="Greens", interpolation="nearest")
    ax.set_title(f"coverage = {bubble.mean():.2f}", fontsize=10)
    ax.axis("off")

plt.tight_layout()
plt.savefig("morphology.png", dpi=150)
plt.show()
