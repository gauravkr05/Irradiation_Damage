#!/usr/bin/env python3
# ============================================================
#  plot_percolation.py
#
#  Plots the output of gb_bubble_percolation.cpp.
#
#  Main figure (the sketch) :
#      V_largest_bubble / V_total_bubble   vs   time
#  flat and low while the bubbles are separate (it sits near 1/N),
#  then a steep jump to 1 when they coalesce into one percolated
#  tunnel network.  One curve per run; the arrow marks the jump.
#
#  Also produced, since the runs already write them :
#      - Omega_V (vented GB coverage)  vs  t      [paper Fig. 6a]
#      - Omega_V                       vs  Omega_C [paper Fig. 6b]
#      - porosity and Omega_C          vs  t      [paper Fig. 4b]
#      - bubble density                vs  t      [paper Fig. 5a]
#      - bubble density vs average bubble area + White's Eq. 15 [Fig. 5b]
#      - projected bubble maps                    [paper Fig. 4a,c-e]
#
#  usage :
#      python3 plot_percolation.py                      # every run_* folder
#      python3 plot_percolation.py run_theta40_rho20    # only these
#      python3 plot_percolation.py --x hours            # x axis in hours
#      python3 plot_percolation.py --validate           # Fick's law check
# ============================================================
import os
import sys
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- house style ----------------
plt.rcParams.update({
    "font.size":        11,
    "axes.labelsize":   12,
    "axes.titlesize":   12,
    "legend.fontsize":  10,
    "axes.grid":        True,
    "grid.alpha":       0.25,
    "grid.linestyle":   ":",
    "figure.dpi":       120,
    "savefig.dpi":      180,
    "savefig.bbox":     "tight",
})

COLS = ["t", "t_h", "porosity", "n_bub", "vmax_frac",
        "omega_c", "omega_v", "avg_area", "rho_b"]

T0_S = 625.0          # paper's quoted time scale, seconds


# ============================================================
#  loading
# ============================================================
def load_run(folder):
    """Read <folder>/metrics.dat into a dict of arrays + metadata."""
    path = os.path.join(folder, "metrics.dat")
    if not os.path.isfile(path):
        return None

    header = ""
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                header += line
            else:
                break

    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.size == 0 or data.shape[0] < 2:
        return None

    run = {c: data[:, i] for i, c in enumerate(COLS[:data.shape[1]])}
    run["folder"] = folder
    run["header"] = header

    # theta / rho0 out of the folder name, e.g. run_theta40_rho20
    base = os.path.basename(folder.rstrip("/"))
    theta, rho = None, None
    for part in base.split("_"):
        if part.startswith("theta"):
            theta = float(part[5:])
        elif part.startswith("rho"):
            rho = float(part[3:])
    run["theta"] = theta
    run["rho0"] = rho

    bits = []
    if theta is not None:
        bits.append(r"$\theta=%g^\circ$" % theta)
    if rho is not None:
        bits.append(r"$\rho_b^o=%g\,\mu m^{-2}$" % rho)
    run["label"] = ", ".join(bits) if bits else base
    return run


def xaxis(run, mode):
    if mode == "hours":
        return run["t_h"], "time  (h,  $t_o$ = %g s)" % T0_S
    return run["t"], "time  (reduced units,  $t_o = \\tilde{l}^2/D_g$)"


def crossing(x, y, level):
    """First x where y rises through `level`, linearly interpolated."""
    for i in range(1, len(y)):
        if y[i - 1] < level <= y[i]:
            if y[i] == y[i - 1]:
                return x[i]
            f = (level - y[i - 1]) / (y[i] - y[i - 1])
            return x[i - 1] + f * (x[i] - x[i - 1])
    return None


# ============================================================
#  Figure 1 : the sketch -- V_largest / V_total vs time
# ============================================================
def fig_largest_fraction(runs, xmode, out="fig_largest_bubble_fraction.png"):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    cmap = plt.cm.viridis(np.linspace(0.05, 0.85, max(len(runs), 2)))

    for r, col in zip(runs, cmap):
        x, xlabel = xaxis(r, xmode)
        y = r["vmax_frac"]
        ax.plot(x, y, "-", lw=2.0, color=col, label=r["label"])

        # mark the jump : first time the largest bubble holds half the gas
        t50 = crossing(x, y, 0.5)
        if t50 is not None:
            ax.plot([t50], [0.5], "o", ms=5, color=col, zorder=5)
            ax.annotate("", xy=(t50, 0.0), xytext=(t50, 0.5),
                        arrowprops=dict(arrowstyle="->", color=col,
                                        lw=1.2, alpha=0.75))

    ax.axhline(1.0, color="0.4", ls="--", lw=1.0)
    ax.text(0.015, 1.0, " fully percolated", transform=ax.get_yaxis_transform(),
            va="bottom", ha="left", fontsize=9, color="0.35")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$V_{\mathrm{largest\ bubble}}\,/\,V_{\mathrm{total\ bubble}}$")
    ax.set_title("Bubble coalescence and percolation on a grain boundary")
    ax.set_ylim(0, 1.08)
    ax.set_xlim(left=0)
    ax.legend(loc="lower right", frameon=True, framealpha=0.9)
    fig.savefig(out)
    plt.close(fig)
    print("wrote", out)


# ============================================================
#  Figure 2 : vented coverage  (paper Fig. 6)
# ============================================================
def fig_vented(runs, xmode, out="fig_vented_coverage.png"):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    cmap = plt.cm.viridis(np.linspace(0.05, 0.85, max(len(runs), 2)))

    for r, col in zip(runs, cmap):
        x, xlabel = xaxis(r, xmode)
        axes[0].plot(x, r["omega_v"], "-", lw=2.0, color=col, label=r["label"])
        axes[1].plot(r["omega_c"], r["omega_v"], "-", lw=2.0, color=col)

    axes[0].set_xlabel(xlabel)
    axes[0].set_ylabel(r"vented GB coverage  $\Omega^V_{GB}$")
    axes[0].set_title("(a)  vs time", loc="left")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_xlim(left=0)
    axes[0].legend(loc="lower right")

    axes[1].axvspan(0.4, 0.6, color="0.85", zorder=0)
    axes[1].text(0.5, 0.05, "0.4 - 0.6", ha="center", fontsize=9, color="0.4")
    axes[1].set_xlabel(r"GB coverage  $\Omega^C_{GB}$")
    axes[1].set_ylabel(r"vented GB coverage  $\Omega^V_{GB}$")
    axes[1].set_title("(b)  vs coverage", loc="left")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_xlim(0, 1.0)

    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print("wrote", out)


# ============================================================
#  Figure 3 : porosity and GB coverage  (paper Fig. 4b)
# ============================================================
def fig_porosity(runs, xmode, out="fig_porosity_coverage.png"):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax2 = ax.twinx()
    ax2.grid(False)
    cmap = plt.cm.viridis(np.linspace(0.05, 0.85, max(len(runs), 2)))

    for r, col in zip(runs, cmap):
        x, xlabel = xaxis(r, xmode)
        ax.plot(x, r["porosity"], "-", lw=1.4, color="0.6", zorder=1)
        ax2.plot(x, r["omega_c"], "-", lw=2.0, color=col,
                 label=r["label"], zorder=2)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("porosity  (grey, all $\\theta$ collapse)", color="0.4")
    ax.tick_params(axis="y", colors="0.4")
    ax2.set_ylabel(r"GB coverage  $\Omega^C_{GB}$  (colour)")
    ax2.set_ylim(0, 1.05)
    ax.set_xlim(left=0)
    ax2.legend(loc="lower right")
    ax.set_title("Porosity is set by the gas budget; coverage by the contact angle")
    fig.savefig(out)
    plt.close(fig)
    print("wrote", out)


# ============================================================
#  Figure 4 : bubble density  (paper Fig. 5)
# ============================================================
def fig_density(runs, xmode, out="fig_bubble_density.png"):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    cmap = plt.cm.viridis(np.linspace(0.05, 0.85, max(len(runs), 2)))

    for r, col in zip(runs, cmap):
        x, xlabel = xaxis(r, xmode)
        axes[0].plot(x, r["rho_b"], "-o", lw=1.8, ms=3, color=col, label=r["label"])
        m = r["avg_area"] > 0
        axes[1].plot(r["avg_area"][m], r["rho_b"][m], "o", ms=4, color=col,
                     label=r["label"])

    axes[0].set_xlabel(xlabel)
    axes[0].set_ylabel(r"bubble density  $\rho_b$  ($\mu m^{-2}$)")
    axes[0].set_title("(a)  coalescence in time", loc="left")
    axes[0].set_xlim(left=0)
    axes[0].legend()

    # White's Eq. 15 :  rho_b = rho_b0 / (1 + 2 rho_b0 (Ab - Ab0))
    Ab = np.logspace(-3, 0, 300)
    for rho0_w, Ab0, ls in ((40.0, 0.006, "--"),):
        pred = rho0_w / (1.0 + 2.0 * rho0_w * (Ab - Ab0))
        ok = pred > 0
        axes[1].plot(Ab[ok], pred[ok], ls, color="k", lw=1.3,
                     label=r"White Eq. 15 ($\rho_b^o$=40, $A_b^o$=0.006)")

    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"average bubble area  $A_b$  ($\mu m^2$)")
    axes[1].set_ylabel(r"bubble density  $\rho_b$  ($\mu m^{-2}$)")
    axes[1].set_title("(b)  universal curve", loc="left")
    axes[1].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print("wrote", out)


# ============================================================
#  Figure 5 : projected bubble maps  (paper Fig. 4a, c-e)
# ============================================================
def fig_maps(runs, out="fig_bubble_maps.png"):
    grids = []
    for r in runs:
        files = sorted(glob.glob(os.path.join(r["folder"], "proj_*.dat")))
        if files:
            grids.append((r, files))
    if not grids:
        return

    ncol = max(len(f) for _, f in grids)
    nrow = len(grids)
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(2.15 * ncol, 2.3 * nrow),
                             squeeze=False)

    for i, (r, files) in enumerate(grids):
        for j in range(ncol):
            ax = axes[i][j]
            ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
            if j >= len(files):
                ax.axis("off")
                continue

            with open(files[j]) as f:
                t = float(f.readline().split("=")[1])
            g = np.loadtxt(files[j], comments="#")
            ax.imshow(g, cmap="bone_r", vmin=0, vmax=1, interpolation="nearest")

            # the venting circle : 90 % of the side length
            n = g.shape[0]
            th = np.linspace(0, 2 * np.pi, 200)
            ax.plot(0.5 * n + 0.45 * n * np.cos(th),
                    0.5 * n + 0.45 * n * np.sin(th),
                    "-", color="tab:red", lw=0.9, alpha=0.8)

            ax.set_title("t = %.0f" % t, fontsize=9)
            if j == 0:
                ax.set_ylabel(r["label"], fontsize=9)

    fig.suptitle("Bubbles projected on the grain boundary plane "
                 "(red = venting circle)", y=1.0)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print("wrote", out)


# ============================================================
#  Fick's-law validation  (paper Fig. 3 inset)
# ============================================================
def fig_validate(path="validate/fick.dat", out="fig_fick_validation.png"):
    if not os.path.isfile(path):
        print("no", path, "- run ./gbbub --validate first")
        return
    d = np.loadtxt(path, comments="#")
    fig, ax = plt.subplots(figsize=(6.6, 4.8))

    # The sharp initial profile relaxes to its equilibrium tanh width before
    # Fickian growth takes over, which costs a fixed offset in r and a fixed
    # lag in t.  Both are properties of the diffuse interface, not of beta_S,
    # so they are taken out by referencing every curve to its first sample.
    t1 = d[:, 1].min()
    for bS in np.unique(d[:, 0]):
        m = d[:, 0] == bS
        r = d[m, 2] - d[m, 2][0]
        ax.plot(d[m, 1] - t1, r, "o", ms=4, label=r"$\beta_S$ = %g" % bS)

    t = np.linspace(0, d[:, 1].max() - t1, 300)
    ax.plot(t, 0.1 * np.sqrt(t), "k--", lw=1.4,
            label=r"Eq. 14 : $\frac{c_g^i-c_g^{eq}}{c_g^b-c_g^{eq}}\sqrt{D_g t}$")

    ax.set_xlabel(r"time since interface relaxed  (reduced units)")
    ax.set_ylabel(r"$r - r_o$  (cells)")
    ax.set_title("Fickian growth recovered for every $\\beta_S$  (Eqs. 13, 14)")
    ax.legend(loc="lower right")
    fig.savefig(out)
    plt.close(fig)
    print("wrote", out)


# ============================================================
def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]

    xmode = "reduced"
    if "--x" in sys.argv:
        i = sys.argv.index("--x")
        if i + 1 < len(sys.argv):
            xmode = sys.argv[i + 1]
            if xmode in args:
                args.remove(xmode)

    if "--validate" in flags:
        fig_validate()
        if not args:
            return

    folders = args if args else sorted(glob.glob("run_*"))
    runs = [r for r in (load_run(f) for f in folders) if r is not None]

    if not runs:
        print("no run_*/metrics.dat found. Run the simulation first, e.g.")
        print("    ./gbbub 40 20 60 1200 60 0.20")
        return

    # order the legend by contact angle, then density
    runs.sort(key=lambda r: (r["theta"] if r["theta"] is not None else 0,
                             r["rho0"] if r["rho0"] is not None else 0))

    print("loaded %d run(s):" % len(runs))
    for r in runs:
        t50 = crossing(r["t"], r["vmax_frac"], 0.5)
        v50 = crossing(r["t"], r["omega_v"], 0.5)
        print("  %-24s  t_end=%7.1f   t(Vmax/Vtot=0.5)=%s   t(OmegaV=0.5)=%s"
              % (os.path.basename(r["folder"]), r["t"][-1],
                 ("%.1f" % t50) if t50 else "  -  ",
                 ("%.1f" % v50) if v50 else "  -  "))

    fig_largest_fraction(runs, xmode)
    fig_vented(runs, xmode)
    fig_porosity(runs, xmode)
    fig_density(runs, xmode)
    fig_maps(runs)
    fig_validate()


if __name__ == "__main__":
    main()