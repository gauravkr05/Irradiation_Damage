#!/usr/bin/env python3
# ============================================================================
#  Area vs time for the shrinking circular grain.
#  Reads f5_area.dat (columns: step  time  area) and plots the measured area.
#
#  Usage:  python plot_area.py            (defaults to f5_area.dat)
#          python plot_area.py f5_area.dat
# ============================================================================
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    fn = sys.argv[1] if len(sys.argv) > 1 else "f5_area.dat"
    if not os.path.exists(fn):
        raise SystemExit(f"{fn} not found. Compile and run the Fig.5 C++ first.")

    data = np.loadtxt(fn, skiprows=1)          # step time area
    t, area = data[:, 1], data[:, 2]

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(t, area, "o", ms=3, color="#1f77b4")

    ax.set_xlabel("time  t")
    ax.set_ylabel(r"grain area  $A$")
    ax.set_title("Curvature-driven shrinking: area vs time")
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig("area_vs_time.png", dpi=150)
    plt.close(fig)
    print("wrote area_vs_time.png")

if __name__ == "__main__":
    main()