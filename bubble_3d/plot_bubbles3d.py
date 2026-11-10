
"""
plot_bubbles3d.py

Plots the results of:

    phasefield_gb_bubbles3d.cpp

It creates:

1. largest_area_fraction_vs_coverage.png

       A_largest / A_total
                   vs
              GB coverage

2. morphology_step_XXXXXX.png

       3-D green bubble morphology

The morphology is extracted from:

       c = 0.5

which is the bubble/matrix phase-field interface.

Required Python packages:

    numpy
    pandas
    matplotlib
    scikit-image

Install scikit-image if needed:

    pip install scikit-image

Run:

    python plot_bubbles3d.py
"""

import glob
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Grid
# ------------------------------------------------------------

NX = 128
NY = 128
NZ = 64

# One grid spacing = 10 nm.
DX_NM = 10.0

# ------------------------------------------------------------
# Read metrics
# ------------------------------------------------------------

METRICS_FILE = "bubble_metrics.csv"

data = pd.read_csv(METRICS_FILE)


# ------------------------------------------------------------
# Plot 1
#
# Largest projected bubble area / total projected bubble area
# versus grain-boundary coverage.
# ------------------------------------------------------------

plt.figure(figsize=(7, 5))

plt.plot(
    data["coverage"],
    data["largest_area_fraction"],
    marker="o",
    markersize=4,
    linewidth=1.5
)

plt.xlabel(
    "Grain-boundary coverage"
)

plt.ylabel(
    "Largest projected area / total projected area"
)

plt.title(
    "Largest bubble fraction vs grain-boundary coverage"
)

plt.grid(
    True,
    alpha=0.3
)

plt.xlim(
    0.0,
    max(
        1.0,
        data["coverage"].max() * 1.05
    )
)

plt.ylim(
    0.0,
    1.05
)

plt.tight_layout()

plt.savefig(
    "largest_area_fraction_vs_coverage.png",
    dpi=300
)

plt.show()


# ------------------------------------------------------------
# Morphology plotting
# ------------------------------------------------------------

try:
    from skimage.measure import marching_cubes
except ImportError:
    print(
        "\nscikit-image is required for the 3-D morphology plot."
    )

    print(
        "Install it with:"
    )

    print(
        "    pip install scikit-image"
    )

    raise SystemExit


# ------------------------------------------------------------
# Read one RAW file
# ------------------------------------------------------------

def read_morphology(filename):

    number_of_values = NX * NY * NZ

    field = np.fromfile(
        filename,
        dtype=np.float32
    )

    if len(field) != number_of_values:
        raise ValueError(
            "Wrong RAW file size: "
            + filename
        )

    # C++ writes:
    #
    # index = (z*NY+y)*NX+x
    #
    # therefore:
    #
    # shape = (z,y,x)
    #
    field = field.reshape(
        (NZ, NY, NX)
    )

    return field


# ------------------------------------------------------------
# Plot one 3-D morphology
# ------------------------------------------------------------

def plot_morphology(filename):

    field = read_morphology(filename)

    # --------------------------------------------------------
    # Extract c = 0.5 isosurface.
    # --------------------------------------------------------

    vertices, faces, normals, values = marching_cubes(
        field,
        level=0.5,
        spacing=(DX_NM, DX_NM, DX_NM)
    )

    # marching_cubes coordinates are:
    #
    # vertices[:,0] = z
    # vertices[:,1] = y
    # vertices[:,2] = x
    #
    z = vertices[:, 0]
    y = vertices[:, 1]
    x = vertices[:, 2]

    fig = plt.figure(
        figsize=(8, 7)
    )

    ax = fig.add_subplot(
        111,
        projection="3d"
    )

    ax.plot_trisurf(
        x,
        y,
        z,
        triangles=faces,
        linewidth=0.0,
        antialiased=True,
        alpha=0.90
    )

    # --------------------------------------------------------
    # Show the GB plane.
    #
    # z = 32*10 nm = 320 nm approximately.
    # --------------------------------------------------------

    x_plane = np.array(
        [0, NX * DX_NM]
    )

    y_plane = np.array(
        [0, NY * DX_NM]
    )

    X, Y = np.meshgrid(
        x_plane,
        y_plane
    )

    Z = np.ones_like(X) * (
        31.5 * DX_NM
    )

    ax.plot_surface(
        X,
        Y,
        Z,
        alpha=0.10
    )

    ax.set_xlabel(
        "x (nm)"
    )

    ax.set_ylabel(
        "y (nm)"
    )

    ax.set_zlabel(
        "z (nm)"
    )

    # Keep the physical aspect ratio.
    ax.set_box_aspect(
        (
            NX,
            NY,
            NZ
        )
    )

    step_text = os.path.basename(
        filename
    )

    step_text = step_text.replace(
        "morphology_step_",
        ""
    )

    step_text = step_text.replace(
        ".raw",
        ""
    )

    ax.set_title(
        "3-D bubble morphology, step = "
        + step_text
        + ", c = 0.5"
    )

    plt.tight_layout()

    output_name = (
        "morphology_"
        + step_text
        + ".png"
    )

    plt.savefig(
        output_name,
        dpi=300
    )

    plt.show()


# ------------------------------------------------------------
# Plot all saved morphology snapshots.
# ------------------------------------------------------------

files = sorted(
    glob.glob(
        "morphology_step_*.raw"
    )
)

if len(files) == 0:

    print(
        "No morphology_step_*.raw files found."
    )

else:

    print(
        "\nFound",
        len(files),
        "morphology files."
    )

    for filename in files:

        print(
            "Plotting:",
            filename
        )

        plot_morphology(
            filename
        )
