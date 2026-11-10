"""
plot_bubble_growth.py

Reads bubble_growth.csv produced by single_bubble_growth.cpp.

It produces:
  1. projected_area_vs_time.png
  2. radius_vs_sqrt_time.png

The second plot compares the simulation with the analytical
sqrt(t) scaling used in the paper:

    r - r0 = [(ci-ceq)/(cb-ceq)] * sqrt(D*t)

For the demonstration:
    ci  = 0.12
    ceq = 0
    cb  = 1
    D   = 1 in dimensionless units.

Run:
    python plot_bubble_growth.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


CSV_FILE = "bubble_growth.csv"

# Same values as the C++ simulation
ci = 0.12
ceq = 0.0
cb = 1.0
D = 1.0

df = pd.read_csv(CSV_FILE)

t = df["time_dimensionless"].to_numpy()
area = df["projected_area_grid2"].to_numpy()
radius = df["equivalent_radius_grid"].to_numpy()

# Initial equivalent radius measured from the first numerical point.
r0 = radius[0]

# ---------------------------------------------------------
# Plot 1: projected bubble area versus time
# ---------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(t, area, marker="o", markersize=3)
plt.xlabel("Dimensionless time")
plt.ylabel("Projected bubble area (grid$^2$)")
plt.title("Single bubble growth: projected area")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("projected_area_vs_time.png", dpi=300)
plt.show()

# ---------------------------------------------------------
# Plot 2: equivalent radius versus sqrt(time)
# ---------------------------------------------------------
sqrt_t = np.sqrt(t)

# Analytical relation from the paper:
# r - r0 = [(ci-ceq)/(cb-ceq)] sqrt(D*t)
analytical_radius = (
    r0
    + (ci - ceq) / (cb - ceq) * np.sqrt(D * t)
)

plt.figure(figsize=(7, 5))
plt.plot(
    sqrt_t,
    radius,
    marker="o",
    markersize=3,
    label="Phase-field simulation",
)
plt.plot(
    sqrt_t,
    analytical_radius,
    linestyle="--",
    label=r"Analytical: $r-r_0=\frac{c_i-c_{eq}}{c_b-c_{eq}}\sqrt{D t}$",
)
plt.xlabel(r"$\sqrt{t}$")
plt.ylabel("Equivalent bubble radius (grid units)")
plt.title("Single bubble growth: comparison with sqrt(t) law")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("radius_vs_sqrt_time.png", dpi=300)
plt.show()

print("Saved:")
print("  projected_area_vs_time.png")
print("  radius_vs_sqrt_time.png")
