import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Read bubble_data.dat
# ============================================================

filename = "bubble_data.dat"

data = np.loadtxt(
    filename,
    comments="#"
)


# ============================================================
# Columns
#
# 0 = step
# 1 = time
# 2 = coverage
# 3 = Alargest / Atotal
# 4 = Atotal
# 5 = Nbubbles
# ============================================================

step = data[:, 0]

time = data[:, 1]

coverage = data[:, 2]

largest_ratio = data[:, 3]

Atotal = data[:, 4]

Nbubbles = data[:, 5]


# ============================================================
# Plot 1
#
# Largest bubble area / total bubble area
# vs
# GB coverage
# ============================================================

plt.figure(figsize=(7, 5))

plt.plot(
    coverage,
    largest_ratio,
    marker="o",
    markersize=3,
    linewidth=1.5
)

plt.xlabel(
    "Bubble coverage"
)

plt.ylabel(
    r"$A_{\mathrm{largest}} / A_{\mathrm{total}}$"
)

plt.title(
    "Largest Bubble Area vs Bubble Coverage"
)

plt.grid(True, alpha=0.3)

plt.tight_layout()

plt.savefig(
    "largest_area_vs_coverage.png",
    dpi=300
)

plt.show()


# ============================================================
# Plot 2
#
# Coverage vs time
# ============================================================

plt.figure(figsize=(7, 5))

plt.plot(
    time,
    coverage,
    linewidth=1.5
)

plt.xlabel(
    "Simulation time"
)

plt.ylabel(
    "Bubble coverage"
)

plt.title(
    "Bubble Coverage vs Simulation Time"
)

plt.grid(True, alpha=0.3)

plt.tight_layout()

plt.savefig(
    "coverage_vs_time.png",
    dpi=300
)

plt.show()


# ============================================================
# Plot 3
#
# Largest bubble fraction vs time
# ============================================================

plt.figure(figsize=(7, 5))

plt.plot(
    time,
    largest_ratio,
    linewidth=1.5
)

plt.xlabel(
    "Simulation time"
)

plt.ylabel(
    r"$A_{\mathrm{largest}} / A_{\mathrm{total}}$"
)

plt.title(
    "Largest Bubble Area Fraction vs Simulation Time"
)

plt.grid(True, alpha=0.3)

plt.tight_layout()

plt.savefig(
    "largest_area_fraction_vs_time.png",
    dpi=300
)

plt.show()


# ============================================================
# Plot 4
#
# Number of bubbles vs coverage
# ============================================================

plt.figure(figsize=(7, 5))

plt.plot(
    coverage,
    Nbubbles,
    marker="o",
    markersize=3,
    linewidth=1.5
)

plt.xlabel(
    "Bubble coverage"
)

plt.ylabel(
    "Number of bubbles"
)

plt.title(
    "Number of Bubbles vs Bubble Coverage"
)

plt.grid(True, alpha=0.3)

plt.tight_layout()

plt.savefig(
    "number_of_bubbles_vs_coverage.png",
    dpi=300
)

plt.show()