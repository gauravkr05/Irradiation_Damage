import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("bubble_ratio.dat")

coverage = data[:, 0]
ratio = data[:, 1]

plt.figure(figsize=(7,5))

plt.plot(
    coverage,
    ratio,
    'o-',
    markersize=3
)

plt.xlabel("Bubble coverage")
plt.ylabel("Largest bubble volume / Total bubble volume")

plt.xlim(0, 1)
plt.ylim(0, 1)

plt.grid(True)

plt.tight_layout()

plt.savefig(
    "largest_bubble_vs_coverage.png",
    dpi=300
)

plt.show()