"""
measure_contact_angles.py

Reads the 2-D lens equilibration runs (angle/circle_bGB*, angle/flat_bGB*),
measures the dihedral angle through the bubble at every snapshot, and writes
results/case_properties.json with, per case:
    gGB               1-D equilibrium GB excess energy (lattice, GB between cells)
    theta_young       Young prediction 2 acos(gGB / 2 gPG)
    theta_meas        measured angle: mean of the final circle-start and
                      flat-start values (both methods A and B)
    theta_meas_spread half the spread between the two starts and methods
Also writes results/contact_angle_history.csv and a figure of the
equilibrated lenses.
"""
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, "angle")
sys.path.insert(0, "analysis")
from measure_angle import measure, NX, NZ, ZGB  # noqa: E402

CASES = ["1.2", "1.5", "2.0", "3.0"]
os.makedirs("results", exist_ok=True)
T = json.load(open("analysis/interface_table.json"))
gPG = T["gPG_lattice"]

hist = []
props = {"gamma_PG": gPG, "gamma_PG_continuum": T["gPG_cont"]}
finals = {}
for b in CASES:
    vals = []
    for start in ("circle", "flat"):
        d = f"angle/{start}_bGB{b}"
        rows = list(csv.DictReader(open(os.path.join(d, "snapshots.csv"))))
        for r in rows:
            c = np.fromfile(os.path.join(d, r["c_file"]), dtype=np.float32).reshape(NZ, NX)
            thA, thB, h, a, _ = measure(c)
            t = float(r["time_dimensionless"])
            hist.append(dict(bGB=b, start=start, time=t, h=h, a=a, theta_A=thA, theta_B=thB))
        finals[(b, start)] = c
        vals += [thA, thB]
    g = T[b]["gGB_mid"]
    props[b] = dict(
        gGB=g,
        gGB_continuum=T[b]["gGB_cont"],
        gb_width_cells=T[b]["width"],
        theta_young=float(np.degrees(2 * np.arccos(g / (2 * gPG)))),
        theta_meas=float(np.mean(vals)),
        theta_meas_spread=float((max(vals) - min(vals)) / 2),
        theta_final_values=[float(v) for v in vals],
    )

with open("results/contact_angle_history.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(hist[0].keys()))
    w.writeheader()
    w.writerows(hist)
json.dump(props, open("results/case_properties.json", "w"), indent=1)

for b in CASES:
    p = props[b]
    print(f"bGB={b}: gGB={p['gGB']:.3f}  Young={p['theta_young']:.1f}  measured={p['theta_meas']:.1f} "
          f"+/- {p['theta_meas_spread']:.1f}  (circle A,B / flat A,B = "
          + ", ".join(f"{v:.1f}" for v in p["theta_final_values"]) + ")")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 4, figsize=(15, 5.5))
for i, start in enumerate(("circle", "flat")):
    for jx, b in enumerate(CASES):
        ax = axs[i, jx]
        c = finals[(b, start)]
        ax.imshow(c[16:80, 30:130], origin="lower", cmap="viridis", vmin=0, vmax=1,
                  extent=[30, 130, 16, 80])
        ax.contour(np.arange(NX), np.arange(NZ), c, levels=[0.5], colors="w", linewidths=0.8)
        ax.axhline(ZGB, color="r", lw=0.6, ls="--")
        ax.set_xlim(30, 130); ax.set_ylim(16, 80)
        ax.set_title(f"bGB={b}, {start} start\n$\\theta$ = {props[b]['theta_meas']:.0f}° (meas.)", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("Equilibrated 2-D lenses at t = 400 (c field, c = 0.5 contour; red = GB plane)")
fig.tight_layout()
fig.savefig("results/contact_angle_lenses.png", dpi=150)
