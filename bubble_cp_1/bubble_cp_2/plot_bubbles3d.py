"""
plot_bubbles3d.py

Plots the output of phasefield_gb_bubbles3d.cpp.

Creates:
    largest_area_fraction_vs_coverage.png   A_largest / A_total  vs  GB coverage
    porosity_coverage_vs_time.png           porosity and coverage vs time
    morphology_step_N.png                   3-D green bubbles on the GB (c = 0.5)
    morphology_panel.png                    several snapshots side by side
    topview_panel.png                       projected bubbles seen from above

Needs: numpy pandas matplotlib scikit-image
    pip install numpy pandas matplotlib scikit-image

Run it with the folder that holds bubble_metrics.csv and the .raw files:
    python plot_bubbles3d.py /path/to/output_folder
or cd into that folder and run:
    python plot_bubbles3d.py
All PNGs are written into that same folder.
"""

import glob
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from skimage.measure import marching_cubes

NX, NY, NZ = 128, 128, 64
DX_NM = 10.0
Z_GB = 31.5

GREEN = "#2fbf3a"
GB_BLUE = "#bcd7f5"


# ------------------------------------------------------------
# 1. Largest projected area / total projected area vs coverage
# ------------------------------------------------------------

folder = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
os.chdir(folder)
print("Reading results from:", folder)

if not os.path.exists("bubble_metrics.csv"):
    raise SystemExit(
        "bubble_metrics.csv not found in " + folder + "\n"
        "Run the simulation in this folder first, or pass the output folder:\n"
        "    python plot_bubbles3d.py /path/to/output_folder")

data = pd.read_csv("bubble_metrics.csv")

needed = ["step", "time_dimensionless", "coverage", "largest_area_fraction", "porosity"]
missing = [col for col in needed if col not in data.columns]
if missing:
    raise SystemExit(
        "bubble_metrics.csv in " + folder + " is missing columns " + str(missing) + ".\n"
        "This file was written by the OLD simulation. Compile and run the corrected\n"
        "phasefield_gb_bubbles3d.cpp in this folder (it overwrites the CSV and .raw files),\n"
        "then run this script again.")

fig, ax = plt.subplots(figsize=(7, 5.5))
ax.step(data["coverage"], data["largest_area_fraction"],
        where="post", color="navy", linewidth=2)
ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.05)
ax.set_xlabel("Coverage")
ax.set_ylabel(r"$A_{largest\ bubble}\ /\ A_{total\ bubble}$")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("largest_area_fraction_vs_coverage.png", dpi=300)
plt.close(fig)


# ------------------------------------------------------------
# 2. Porosity and coverage vs time
# ------------------------------------------------------------

fig, ax1 = plt.subplots(figsize=(7, 5))
ax1.plot(data["time_dimensionless"], data["porosity"], color="black", linewidth=2,
         label="porosity")
ax1.set_xlabel("time (dimensionless)")
ax1.set_ylabel("porosity")
ax1.set_ylim(0, max(0.2, data["porosity"].max() * 1.1))
ax2 = ax1.twinx()
ax2.plot(data["time_dimensionless"], data["coverage"], color="gray", linewidth=2,
         linestyle="--", label=r"coverage $X_{GB}^C$")
ax2.set_ylabel(r"GB coverage $X_{GB}^C$")
ax2.set_ylim(0, 1.0)
lines = ax1.get_lines() + ax2.get_lines()
ax1.legend(lines, [l.get_label() for l in lines], loc="upper left")
ax1.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("porosity_coverage_vs_time.png", dpi=300)
plt.close(fig)


# ------------------------------------------------------------
# 3-D morphology helpers
# ------------------------------------------------------------

def read_field(filename):
    field = np.fromfile(filename, dtype=np.float32)
    if field.size != NX * NY * NZ:
        raise ValueError("Wrong RAW file size: " + filename)
    # C++ index = (z*NY + y)*NX + x  ->  numpy shape (z, y, x)
    return field.reshape((NZ, NY, NX))


def step_of(filename):
    return int(os.path.basename(filename).split("_")[-1].split(".")[0])


def draw_morphology(ax, field, title):
    # c = 0.5 isosurface. marching_cubes returns (z, y, x) coordinates.
    verts, faces, _, _ = marching_cubes(field, level=0.5,
                                        spacing=(DX_NM, DX_NM, DX_NM))
    z = verts[:, 0]
    y = verts[:, 1]
    x = verts[:, 2]

    # GB plane
    X, Y = np.meshgrid([0, NX * DX_NM], [0, NY * DX_NM])
    ax.plot_surface(X, Y, np.full_like(X, Z_GB * DX_NM, dtype=float),
                    color=GB_BLUE, alpha=0.55, shade=False, zorder=0)

    ax.plot_trisurf(x, y, z, triangles=faces, color=GREEN,
                    linewidth=0.0, antialiased=True, shade=True,
                    lightsource=LightSource(azdeg=315, altdeg=50))

    half = 16 * DX_NM
    ax.set_xlim(0, NX * DX_NM)
    ax.set_ylim(0, NY * DX_NM)
    ax.set_zlim(Z_GB * DX_NM - half, Z_GB * DX_NM + half)
    ax.set_box_aspect((NX * DX_NM, NY * DX_NM, 2 * half))
    ax.view_init(elev=40, azim=-60)
    ax.set_axis_off()
    ax.set_title(title)


files = sorted(glob.glob("morphology_step_*.raw"), key=step_of)
if not files:
    raise SystemExit("No morphology_step_*.raw files found.")

cov_by_step = dict(zip(data["step"], data["coverage"]))


def label(fn):
    s = step_of(fn)
    cov = cov_by_step.get(s)
    return f"step {s}" + (f",  coverage {cov:.2f}" if cov is not None else "")


# ------------------------------------------------------------
# 3. One 3-D image per snapshot
# ------------------------------------------------------------

for fn in files:
    print("Plotting", fn)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    draw_morphology(ax, read_field(fn), label(fn))
    fig.tight_layout()
    fig.savefig(f"morphology_step_{step_of(fn)}.png", dpi=200)
    plt.close(fig)


# ------------------------------------------------------------
# 4. Panel of 6 snapshots (evenly spread through the run)
# ------------------------------------------------------------

pick = [files[int(round(i))] for i in np.linspace(0, len(files) - 1, min(6, len(files)))]
fig = plt.figure(figsize=(15, 9))
for k, fn in enumerate(pick):
    ax = fig.add_subplot(2, 3, k + 1, projection="3d")
    draw_morphology(ax, read_field(fn), label(fn))
fig.tight_layout()
fig.savefig("morphology_panel.png", dpi=200)
plt.close(fig)


# ------------------------------------------------------------
# 5. Top view: projected bubble area on the GB (what the
#    coverage and A_largest/A_total are measured from)
# ------------------------------------------------------------

fig, axes = plt.subplots(2, 3, figsize=(13, 9))
for ax, fn in zip(axes.flat, pick):
    proj = (read_field(fn) >= 0.5).any(axis=0)
    ax.imshow(proj, origin="lower", cmap=matplotlib.colors.ListedColormap([GB_BLUE, GREEN]),
              extent=[0, NX * DX_NM, 0, NY * DX_NM])
    ax.set_title(label(fn))
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
for ax in list(axes.flat)[len(pick):]:
    ax.axis("off")
fig.tight_layout()
fig.savefig("topview_panel.png", dpi=200)
plt.close(fig)

print("Done.")