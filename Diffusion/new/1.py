import numpy as np
import matplotlib.pyplot as plt

save_times = [0, 200, 500, 1000]

# middle row slice from 2D: y = 50 (j = ny/2)
slice_y = 50.0

fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

# --- left plot: 1D explicit ---
ax = axes[0]
for tm in save_times:
    data = np.loadtxt("time_%d/c_time_%d.dat" % (tm, tm), skiprows=1)
    ax.plot(data[:, 0], data[:, 1], label="time = %d" % tm)
ax.set_xlabel("x")
ax.set_ylabel("c")
ax.set_title("1D Explicit")
ax.set_ylim(-0.05, 1.05)
ax.legend()
ax.grid(True, alpha=0.3)

# --- right plot: middle-row slice of 2D explicit ---
ax = axes[1]
for tm in save_times:
    data = np.loadtxt("time2d_%d/c_time_%d.dat" % (tm, tm), skiprows=1)
    # extract slice at y = slice_y
    mask = np.isclose(data[:, 1], slice_y)
    row  = data[mask]
    ax.plot(row[:, 0], row[:, 2], label="time = %d" % tm)
ax.set_xlabel("x")
ax.set_title("2D Explicit  (slice at y = %g)" % slice_y)
ax.set_ylim(-0.05, 1.05)
ax.legend()
ax.grid(True, alpha=0.3)

plt.suptitle("1D vs 2D Explicit Diffusion  —  middle row slice", fontsize=12)
plt.tight_layout()
plt.savefig("compare_1d_vs_2d.png", dpi=150, bbox_inches="tight")
plt.show()