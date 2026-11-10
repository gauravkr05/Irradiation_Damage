import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("bubble_data.dat", comments="#")
cov, ratio = d[:, 2], d[:, 3]

fig, ax = plt.subplots(figsize=(6.2, 5.0))
ax.plot(cov, ratio, "-", lw=1.8, color="navy")
ax.set_xlabel("Coverage")
ax.set_ylabel(r"$A_{largest\ bubble}\ /\ A_{total\ bubble}$")
ax.set_xlim(0, 1.0)
ax.set_ylim(0, 1.05)
ax.axhline(1.0, color="0.5", lw=0.8, ls="--")
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("bubble_ratio_vs_coverage.png", dpi=200)
plt.show()
