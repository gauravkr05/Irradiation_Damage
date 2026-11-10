#!/usr/bin/env python3
"""
analyze_growth.py

Measure single-bubble growth rate from one or more runs of single_bubble_gb.

Usage:
    python3 analyze_growth.py run1 run2 run3
    python3 analyze_growth.py run_* --tskip 2 --tmax 50

Each argument is an output folder containing growth.csv and params.txt.

What it reports for every run (inside the fit window tskip <= t <= tmax):

    dR/dt      mean growth velocity of the equivalent radius
               (cells per unit time, and nm/s using t_phys_s)
    dV/dt      mean volume growth rate (cells^3 per unit time)
    k          parabolic rate constant from  R^2 = R^2(tskip) + k (t - tskip)
               (diffusion-controlled growth at constant supersaturation
                gives a straight line here)
    m          exponent in  dV ~ dt^m
               m ~ 1.5 : bulk-diffusion controlled (R^2 ~ t)
               m ~ 1   : constant flux (interface / GB-supply limited)
               m < 1   : matrix depleting, growth slowing down
               m < 0 or dV < 0 : bubble below critical radius, shrinking
    gas drift  max |relative change of total gas|. Must stay ~1e-12 or
               smaller; if not, conservation is broken (dt too big or
               clampC=1).

It also writes growth_analysis.png with R(t), V(t), dR/dt, c_far(t),
contact-angle estimate and gas drift, and a c-field slice for each run.
"""

import argparse
import glob
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_params(folder):
    p = {}
    with open(os.path.join(folder, "params.txt")) as f:
        for line in f:
            if "=" in line:
                k, v = line.strip().split("=", 1)
                try:
                    p[k] = float(v)
                except ValueError:
                    p[k] = v
    return p


def label_of(folder, p):
    return (f"{os.path.basename(folder.rstrip('/'))}: c0={p['c0']:g}, R0={p['R0']:g}, "
            f"Mb={p['Mb']:g}, th={int(p['theta'])}")


def smooth_derivative(t, y, w=5):
    """central-difference derivative, lightly smoothed with a moving mean"""
    dy = np.gradient(y, t)
    if len(dy) > w:
        k = np.ones(w) / w
        dy = np.convolve(dy, k, mode="same")
    return dy


def analyse(folder, tskip, tmax):
    d = np.genfromtxt(os.path.join(folder, "growth.csv"), delimiter=",", names=True)
    p = read_params(folder)

    t = d["t"]
    tp = d["t_phys_s"]
    R = d["R_eq"]
    V = d["V"]

    sel = (t >= tskip) & (t <= (tmax if tmax > 0 else t.max()))
    res = {"folder": folder, "p": p, "d": d, "sel": sel}

    if sel.sum() < 3:
        res["ok"] = False
        return res
    res["ok"] = True

    ts, tps, Rs, Vs = t[sel], tp[sel], R[sel], V[sel]

    # linear fits
    dRdt = np.polyfit(ts, Rs, 1)[0]
    dVdt = np.polyfit(ts, Vs, 1)[0]
    k = np.polyfit(ts - ts[0], Rs ** 2 - Rs[0] ** 2, 1)[0]

    # physical units (t_phys_s is proportional to t)
    sec_per_t = (tps[-1] - tps[0]) / (ts[-1] - ts[0]) if ts[-1] > ts[0] else np.nan
    dx_nm = 10.0
    dRdt_nm_s = dRdt * dx_nm / sec_per_t

    # growth exponent of the volume increment
    dV = Vs - Vs[0]
    dt_ = ts - ts[0]
    m = np.nan
    good = (dV > 0) & (dt_ > 0)
    if good.sum() >= 3:
        m = np.polyfit(np.log(dt_[good]), np.log(dV[good]), 1)[0]

    # mass balance: largest volume the bubble can reach if all excess
    # gas in the box (above c_eq ~ 0 in the matrix) went into it
    Ntot = p["NX"] * p["NY"] * p["NZ"]
    Vmax = p["c0"] * Ntot
    Rmax = (3 * Vmax / (4 * np.pi)) ** (1 / 3)

    res.update(dict(dRdt=dRdt, dVdt=dVdt, k=k, m=m, dRdt_nm_s=dRdt_nm_s,
                    drift=np.nanmax(np.abs(d["gas_drift_rel"])),
                    Rmax=Rmax, R_end=R[-1], theta_end=d["theta_est_deg"][-1],
                    cfar_end=d["c_far"][-1]))
    return res


def plot_slice(ax, folder, p):
    raws = sorted(glob.glob(os.path.join(folder, "c_*.raw")))
    if not raws:
        ax.text(0.5, 0.5, "no raw file", ha="center", va="center")
        return
    NX, NY, NZ = int(p["NX"]), int(p["NY"]), int(p["NZ"])
    c = np.fromfile(raws[-1], dtype=np.float32).reshape(NZ, NY, NX)
    im = ax.imshow(c[:, NY // 2, :], origin="lower", cmap="viridis", vmin=0, vmax=1)
    ax.contour(c[:, NY // 2, :], levels=[0.5], colors="w", linewidths=0.8)
    ax.axhline(NZ / 2 - 0.5, color="r", lw=0.6, ls="--")
    ax.set_title(f"{os.path.basename(folder.rstrip('/'))}: c(x,z) at y=NY/2\n"
                 f"{os.path.basename(raws[-1])}  (red = GB)", fontsize=8)
    ax.set_xlabel("x (cells)")
    ax.set_ylabel("z (cells)")
    plt.colorbar(im, ax=ax, fraction=0.046)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--tskip", type=float, default=2.0,
                    help="ignore t < tskip (initial interface relaxation)")
    ap.add_argument("--tmax", type=float, default=-1.0,
                    help="ignore t > tmax (e.g. once the matrix is depleted)")
    ap.add_argument("--physical", action="store_true",
                    help="plot against physical time (s) instead of model time")
    ap.add_argument("--out", default="growth_analysis.png")
    a = ap.parse_args()

    results = [analyse(f, a.tskip, a.tmax) for f in a.runs
               if os.path.exists(os.path.join(f, "growth.csv"))]
    if not results:
        sys.exit("No growth.csv found in the given folders.")

    # ---------------- table ----------------
    hdr = (f"{'run':<22}{'c0':>6}{'R0':>6}{'Mb':>7}{'th':>4}"
           f"{'dR/dt':>10}{'nm/s':>11}{'dV/dt':>10}{'k(R^2)':>10}{'m':>7}"
           f"{'R_end':>8}{'R_max':>8}{'th_est':>8}{'drift':>10}")
    print(f"\nFit window: t >= {a.tskip}" + (f", t <= {a.tmax}" if a.tmax > 0 else ""))
    print(hdr)
    print("-" * len(hdr))
    for r in results:
        p = r["p"]
        name = os.path.basename(r["folder"].rstrip("/"))[:21]
        if not r["ok"]:
            print(f"{name:<22} not enough points in fit window")
            continue
        print(f"{name:<22}{p['c0']:>6.3g}{p['R0']:>6.3g}{p['Mb']:>7.3g}{int(p['theta']):>4d}"
              f"{r['dRdt']:>10.4f}{r['dRdt_nm_s']:>11.3e}{r['dVdt']:>10.2f}{r['k']:>10.3f}"
              f"{r['m']:>7.2f}{r['R_end']:>8.2f}{r['Rmax']:>8.2f}{r['theta_end']:>8.1f}"
              f"{r['drift']:>10.1e}")
        if r["dRdt"] < 0:
            print(f"{'':<22}-> shrinking: R0 below critical radius for this c0"
                  " (or initN=0 left n=1 inside the bubble)")
        if r["drift"] > 1e-8:
            print(f"{'':<22}-> gas not conserved: lower dt or set clampC=0")
        if r["R_end"] > 0.7 * r["Rmax"]:
            print(f"{'':<22}-> close to mass-balance limit: box too small or c0 too low")
    print("\nUnits: R in cells (1 cell = 10 nm), t in model time units;"
          " nm/s uses the bulk-diffusivity time scaling in t_phys_s.\n")

    # ---------------- figure ----------------
    n = len(results)
    fig = plt.figure(figsize=(15, 9 + 3.2 * ((n + 2) // 3)))
    gs = fig.add_gridspec(2 + (n + 2) // 3, 3)
    axR = fig.add_subplot(gs[0, 0])
    axV = fig.add_subplot(gs[0, 1])
    axv = fig.add_subplot(gs[0, 2])
    axc = fig.add_subplot(gs[1, 0])
    axt = fig.add_subplot(gs[1, 1])
    axg = fig.add_subplot(gs[1, 2])

    for r in results:
        d, p = r["d"], r["p"]
        t = d["t_phys_s"] if a.physical else d["t"]
        lab = label_of(r["folder"], p)
        axR.plot(t, d["R_eq"], label=lab)
        axV.plot(t, d["V"])
        axv.plot(t, smooth_derivative(d["t"], d["R_eq"]))
        axc.plot(t, d["c_far"])
        axt.plot(t, d["theta_est_deg"])
        axg.semilogy(t, np.abs(d["gas_drift_rel"]) + 1e-18)

    xl = "t (s)" if a.physical else "t (model units)"
    for ax, yl, title in [
        (axR, "R_eq (cells)", "Equivalent radius  (3V/4pi)^(1/3)"),
        (axV, "V (cells^3, c >= 0.5)", "Bubble volume"),
        (axv, "dR_eq/dt (cells / model time)", "Growth velocity"),
        (axc, "c far from bubble", "Matrix concentration (depletion)"),
        (axt, "theta estimate (deg)", "Lens angle 2 atan(h/a)"),
        (axg, "|total gas drift| (relative)", "Gas conservation check"),
    ]:
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
    if not a.physical:
        for ax in (axR, axV, axv):
            ax.axvspan(0, a.tskip, color="grey", alpha=0.15)
    axR.legend(fontsize=7)

    for i, r in enumerate(results):
        ax = fig.add_subplot(gs[2 + i // 3, i % 3])
        plot_slice(ax, r["folder"], r["p"])

    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print(f"Figure written to {a.out}")


if __name__ == "__main__":
    main()
