import numpy as np
import matplotlib.pyplot as plt

cases = [
    (0.25, "dt=0.25  r=0.25  stable"),
    (0.40, "dt=0.40  r=0.40  stable"),
    (0.50, "dt=0.50  r=0.50  borderline"),
    (0.60, "dt=0.60  r=0.60  UNSTABLE"),
    (1.00, "dt=1.00  r=1.00  UNSTABLE"),
]

plt.figure(figsize=(9, 5))
for dt, label in cases:
    data = np.loadtxt("stability_dt%g/c.dat" % dt, skiprows=1)
    plt.plot(data[:, 0], data[:, 1], label=label)

plt.xlabel("x")
plt.ylabel("c")
plt.title("Explicit Scheme Stability Check  (physical time = 10)")
plt.ylim(-20, 20)
plt.axhline(0, color="black", linewidth=0.5)
plt.axhline(1, color="black", linewidth=0.5, linestyle="--")
plt.legend(fontsize=8)
plt.grid(True, alpha=0.3)
plt.savefig("stability_check.png", dpi=150, bbox_inches="tight")
plt.show()
