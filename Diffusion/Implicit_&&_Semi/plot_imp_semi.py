import numpy as np
import matplotlib.pyplot as plt
import os

save_times = [0, 200, 500, 1000]
schemes = [
    ("implicit",      "Implicit"),
    ("semi_implicit", "Semi-Implicit "),
]

fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

for ax, (folder_prefix, title) in zip(axes, schemes):
    for tm in save_times:
        fname = "%s_time_%d/c_time_%d.dat" % (folder_prefix, tm, tm)
        data = np.loadtxt(fname, skiprows=1)
        ax.plot(data[:, 0], data[:, 1], label="time = %d" % tm)
    ax.set_xlabel("x")
    ax.set_ylabel("c")
    ax.set_title(title)
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("implicit_vs_semi.png", dpi=150, bbox_inches="tight")
plt.show()
