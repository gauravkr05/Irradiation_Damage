#!/usr/bin/env python3
"""Figures for the UGP report, written to ugp_report/figures/.

Only reads existing simulation output (no new simulations are run):
  f1_single      gb_energy_study/results/contact_angle_lenses.png (cropped),
                 gb_energy_study/preliminary/crit2/*/bubble_metrics.csv
  f2_a1          gb_energy_study/percolation_study/runs/*/bubble_metrics.csv
  f3_morphology  gb_energy_study/percolation_study/runs/N33_*/thick/*.npy
  f4_percolation percolation_study runs + results/{events_summary,gaps}.csv
  f5_flatgb      Flat_Gb/{gb_position_all_runs,velocity_summary}.csv
  f6_coupled     Flat_Gb_bubbles/runs/{coupled_df0.1,gb_only_df0.1,...}
  f7_shift       Flat_Gb_bubbles/results/robustness_summary.json
  f8_poly        Flat_Gb_bubbles/runs/sphere_*/sphere.csv, GB2/fig2_out/*.csv

Run:  python3 ugp_report/make_figures.py
"""
import csv
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PERC = os.path.join(ROOT, "gb_energy_study", "percolation_study")
FGB = os.path.join(ROOT, "Flat_Gb_bubbles")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- style
INK, INK2, MUTED, GRID, BAND = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#dcdbd6"
BLUE, ORANGE, AQUA, YELLOW, VIOLET, RED = ("#2a78d6", "#eb6834", "#1baf7a",
                                           "#eda100", "#6250d6", "#e34948")
# contact-angle cases: fixed order and colour in every figure
ANGLES = [("1.2", 148, BLUE), ("1.5", 128, ORANGE), ("2.0", 120, AQUA), ("3.0", 104, YELLOW)]

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": INK2, "axes.linewidth": 0.6, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "lines.linewidth": 1.3, "legend.frameon": False,
    "savefig.dpi": 300, "pdf.fonttype": 42, "axes.titlelocation": "left",
})


def read_csv(path):
    """CSV -> dict of float arrays (no pandas on this machine)."""
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    head, body = rows[0], [r for r in rows[1:] if r]
    out = {}
    for j, h in enumerate(head):
        col = []
        for r in body:
            try:
                col.append(float(r[j]))
            except (ValueError, IndexError):
                col.append(np.nan)
        out[h] = np.array(col)
    return out


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote", name)


def bare(ax, edge=True):
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(edge); s.set_color(MUTED); s.set_linewidth(0.5)


def ramp(n, lo=0.35, hi=0.95, name="Blues"):
    cm = plt.get_cmap(name)
    return [cm(lo + (hi - lo) * i / max(n - 1, 1)) for i in range(n)]


def events():
    return {r["case"]: r for r in read_rows(os.path.join(PERC, "results", "events_summary.csv"))}


# ================================================ f1: a single bubble on the GB
def fig_single():
    im = Image.open(os.path.join(ROOT, "gb_energy_study", "results",
                                 "contact_angle_lenses.png")).convert("RGB")
    frames = [(28, 551), (585, 1108), (1142, 1665), (1699, 2222)]   # top-row panels
    y0, y1 = 165, 380
    fig = plt.figure(figsize=(6.5, 2.75))
    outer = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.05], wspace=0.28)
    left = outer[0].subgridspec(2, 2, wspace=0.06, hspace=0.32)
    for k, ((x0, x1), (b, th, _)) in enumerate(zip(frames, ANGLES)):
        ax = fig.add_subplot(left[k // 2, k % 2])
        w = x1 - x0
        ax.imshow(np.asarray(im.crop((x0 + int(0.2 * w), y0, x1 - int(0.2 * w), y1))))
        bare(ax, edge=False)
        ax.set_title(("(a)  " if k == 0 else "") + f"$\\theta$ = {th}°", fontsize=7.5)

    ax = fig.add_subplot(outer[1])
    crit = os.path.join(ROOT, "gb_energy_study", "preliminary", "crit2")
    radii = [2.5, 3.0, 3.4, 3.8, 4.5, 5.0, 5.5]
    cols = dict(zip(radii, ramp(len(radii))))
    box = 48 * 48 * 64
    for r in radii:
        m = read_csv(os.path.join(crit, f"bGB1.2_r{r}", "bubble_metrics.csv"))
        ax.plot(m["time_dimensionless"], m["porosity"] * box, color=cols[r], lw=1.4,
                label=f"$r_0$ = {r}")
    m = read_csv(os.path.join(crit, "bGB3.0_r4.5", "bubble_metrics.csv"))
    ax.plot(m["time_dimensionless"], m["porosity"] * box, "--", color=cols[4.5], lw=1.4,
            label="$r_0$ = 4.5, $\\theta$ = 104°")
    ax.text(10.2, 0, "dissolve", fontsize=6.5, color=INK2, va="center")
    ax.set_xlim(0, 12); ax.set_ylim(-80, 2750)
    ax.set_xlabel("time"); ax.set_ylabel("bubble volume (cells)")
    ax.set_title("(b)  one bubble ($\\theta$ = 148°): starting size")
    ax.legend(loc="upper left", fontsize=6.3, ncol=2, handlelength=1.6, columnspacing=0.8)
    save(fig, "f1_single.pdf")


# =============================================== f2: A1/Atot vs coverage
def fig_a1():
    ev = events()
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.55), sharey=True)
    for ax, N, lab in zip(axs, (33, 25), "ab"):
        for b, th, col in ANGLES:
            case = f"N{N}_bGB{b}"
            m = read_csv(os.path.join(PERC, "runs", case, "bubble_metrics.csv"))
            cov, a1 = m["coverage"], m["largest_area_fraction"]
            stop = int(np.argmax(cov)) + 1          # until coverage peaks (it falls later at high theta)
            ax.plot(cov[:stop], a1[:stop], color=col, lw=1.3, label=f"$\\theta$ = {th}°")
            e = ev[case]
            cf = float(e["cov_first_merger"])
            i = int(np.argmin(abs(m["time_dimensionless"] - float(e["t_first_merger"]))))
            ax.plot(cf, a1[i], "o", mfc="white", mec=col, mew=1.2, ms=5, zorder=5)
            if e["t_single"]:
                ax.plot(float(e["cov_single"]), 1.0, "o", color=col, mec="white", mew=0.6, ms=6, zorder=5)
        ax.set_xlim(0.08, 0.76); ax.set_ylim(0, 1.04)
        ax.set_xlabel("coverage (fraction of GB covered by bubbles)")
        ax.set_title(f"({lab})  N = {N} bubbles")
    axs[0].set_ylabel("$A_1/A_{tot}$ (largest bubble / all bubbles)")
    h, l = axs[0].get_legend_handles_labels()
    h += [plt.Line2D([], [], ls="", marker="o", mfc="white", mec=INK2, mew=1.2, ms=5),
          plt.Line2D([], [], ls="", marker="o", color=INK2, ms=6)]
    l += ["percolation starts", "percolation ends"]
    axs[0].legend(h, l, loc="upper left", fontsize=6.5, handlelength=1.6)
    fig.tight_layout(w_pad=1.2)
    save(fig, "f2_a1_coverage.pdf")


# ============================================ f3: top-view morphology maps
def fig_morphology():
    ev = events()
    times = [5, 12, 20]
    fig, axs = plt.subplots(4, 4, figsize=(5.4, 5.55))
    im = None
    for i, (b, th, col) in enumerate(ANGLES):
        case = f"N33_bGB{b}"
        m = read_csv(os.path.join(PERC, "runs", case, "bubble_metrics.csv"))
        tsing = ev[case]["t_single"]
        t_end = math.ceil(float(tsing)) if tsing else 150
        for j, t in enumerate(times + [t_end]):
            ax = axs[i, j]
            a = np.load(os.path.join(PERC, "runs", case, "thick", f"t{t:07.2f}.npy")).astype(float)
            im = ax.imshow(np.ma.masked_less(a, 1), cmap="Blues", vmin=0, vmax=24,
                           origin="lower", interpolation="nearest")
            ax.contour(a, levels=[2.5], colors=INK, linewidths=0.4, origin="lower")
            bare(ax)
            n = int(m["number_of_projected_bubbles"][np.argmin(abs(m["time_dimensionless"] - t))])
            ax.set_title(f"t = {t}, {n} bubble" + ("" if n == 1 else "s"), fontsize=7, loc="center", pad=2)
            if j == 0:
                ax.set_ylabel(f"$\\theta$ = {th}°", fontsize=8.5, color=INK)
    axs[0, 3].annotate("end state", xy=(0.5, 1.2), xycoords="axes fraction",
                       ha="center", fontsize=7.5, color=INK2)
    fig.subplots_adjust(wspace=0.06, hspace=0.22, right=0.88)
    cax = fig.add_axes([0.9, 0.3, 0.018, 0.4])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("bubble thickness (cells)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    cb.outline.set_linewidth(0.4)
    save(fig, "f3_morphology.pdf")


# =========================================== f4: percolation measurements
def fig_percolation():
    ev = events()
    gaps = read_csv(os.path.join(PERC, "results", "gaps.csv"))
    gcase = np.array([r["case"] for r in read_rows(os.path.join(PERC, "results", "gaps.csv"))])
    fig, axs = plt.subplots(2, 3, figsize=(6.5, 4.1))
    letters = iter("abcdef")
    for r, N in enumerate((33, 25)):
        for b, th, col in ANGLES:
            case = f"N{N}_bGB{b}"
            m = read_csv(os.path.join(PERC, "runs", case, "bubble_metrics.csv"))
            t, cov, n = m["time_dimensionless"], m["coverage"], m["number_of_projected_bubbles"]
            e = ev[case]
            lab = f"$\\theta$ = {th}°"
            axs[r, 0].plot(cov, n, color=col, lw=1.1, label=lab)
            axs[r, 1].plot(t, cov, color=col, lw=1.3, label=lab)
            cf, tf = float(e["cov_first_merger"]), float(e["t_first_merger"])
            axs[r, 0].plot(cf, N, "o", mfc="white", mec=col, mew=1.2, ms=5, zorder=5)
            axs[r, 1].plot(tf, cf, "o", mfc="white", mec=col, mew=1.2, ms=5, zorder=5)
            if e["t_single"]:
                cs, ts = float(e["cov_single"]), float(e["t_single"])
                axs[r, 0].plot(cs, 1, "o", color=col, mec="white", mew=0.6, ms=6, zorder=5)
                axs[r, 1].plot(ts, cs, "o", color=col, mec="white", mew=0.6, ms=6, zorder=5)
            k = gcase == case
            axs[r, 2].plot(gaps["t"][k], gaps["min_gap"][k], color=col, lw=0.9, label=lab)
        axs[r, 0].set_xlabel("coverage"); axs[r, 0].set_ylabel("number of bubbles")
        axs[r, 0].set_ylim(0, N + 2)
        axs[r, 1].set_xlabel("time"); axs[r, 1].set_ylabel("coverage"); axs[r, 1].set_xlim(0, 152)
        axs[r, 2].set_xlabel("time"); axs[r, 2].set_ylabel("smallest gap (cells)"); axs[r, 2].set_xlim(0, 152)
        axs[r, 0].set_title(f"({next(letters)})  N = {N}: count vs coverage")
        axs[r, 1].set_title(f"({next(letters)})  N = {N}: coverage vs time")
        axs[r, 2].set_title(f"({next(letters)})  N = {N}: smallest gap")
    h, l = axs[0, 1].get_legend_handles_labels()
    h += [plt.Line2D([], [], ls="", marker="o", mfc="white", mec=INK2, mew=1.2, ms=5),
          plt.Line2D([], [], ls="", marker="o", color=INK2, ms=6)]
    l += ["percolation starts", "percolation ends"]
    fig.legend(h, l, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.035),
               handlelength=1.6, columnspacing=1.2)
    fig.tight_layout(h_pad=1.2, w_pad=1.0)
    save(fig, "f4_percolation.pdf")


# =============================================== f5: flat GB, no bubbles
def fig_flatgb():
    d = os.path.join(ROOT, "Flat_Gb")
    pos = read_csv(os.path.join(d, "gb_position_all_runs.csv"))
    case = np.array([r["case"] for r in read_rows(os.path.join(d, "gb_position_all_runs.csv"))])
    vs = read_rows(os.path.join(d, "velocity_summary.csv"))
    dfs = [("df0", 0.0), ("df0.025", 0.025), ("df0.0375", 0.0375),
           ("df0.05", 0.05), ("df0.075", 0.075), ("df0.1", 0.1)]
    cols = dict(zip([c for c, _ in dfs], ramp(len(dfs), 0.3, 0.95)))
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.45))
    ax = axs[0]
    for c, df in dfs:
        k = case == c
        ax.plot(pos["time"][k], pos["gb_z_mean"][k] - 20, color=cols[c], label=f"$\\Delta f$ = {df:g}")
    ax.set_xlim(0, 252); ax.set_ylim(-0.5, 31)
    ax.set_xlabel("time"); ax.set_ylabel("GB displacement (cells)")
    ax.set_title("(a)  GB position, flat GB without bubbles")
    ax.legend(loc="upper left", fontsize=6.5, ncol=3, handlelength=1.5, columnspacing=0.8)

    ax = axs[1]
    x = np.linspace(0, 0.105, 50)
    ax.plot(x, 1.865 * x, "--", color=INK2, lw=1.0, label="sharp-interface theory")
    ax.axvspan(0, 0.031, color=GRID, lw=0)
    ax.text(0.003, 0.172, "pinned\nby the grid", fontsize=6.5, color=INK2, va="top")
    for r in vs:
        if r["case"] in cols:
            ax.plot(float(r["df_bias"]), float(r["v_fit"]), "o", ms=6, color=cols[r["case"]],
                    mec="white", mew=0.6, zorder=5)
    ax.plot([], [], "o", color=INK2, ms=6, label="simulation")
    ax.set_xlim(0, 0.105); ax.set_ylim(0, 0.2)
    ax.set_xlabel("driving force $\\Delta f$"); ax.set_ylabel("GB velocity (cells / time)")
    ax.set_title("(b)  velocity vs driving force")
    ax.legend(loc="lower right", fontsize=6.5)
    fig.tight_layout(w_pad=1.5)
    save(fig, "f5_flatgb.pdf")


# ===================================== f6: bubbles on a migrating flat GB
def read_map(fn, NX=128, NY=128, NZ=64):
    b = open(fn, "rb").read()
    o = NX * NY + 4 * NX * NY + NX * NY + NX * NY   # skip thick (u8), gbz (f32), zb, zt (u8)
    sc = np.frombuffer(b, np.float32, NZ * NX, o).reshape(NZ, NX); o += 4 * NZ * NX
    sn = np.frombuffer(b, np.float32, NZ * NX, o).reshape(NZ, NX)
    return sc, sn


def fig_coupled():
    fig = plt.figure(figsize=(6.5, 3.75))
    outer = GridSpec(2, 1, figure=fig, height_ratios=[0.62, 1.6], hspace=0.42)
    top = outer[0].subgridspec(1, 4, wspace=0.08)
    bot = outer[1].subgridspec(1, 2, wspace=0.3)
    for k, t in enumerate((0, 10, 25, 100)):
        ax = fig.add_subplot(top[0, k])
        sc, sn = read_map(os.path.join(FGB, "runs", "coupled_df0.1", "maps", f"t{t:07.2f}.bin"))
        ax.imshow(sc, origin="lower", cmap="Greys", vmin=0, vmax=1,
                  extent=(0, 128, 0, 64), interpolation="bilinear")
        ax.contour(np.arange(128) + 0.5, np.arange(64) + 0.5,
                   np.ma.masked_where(sc > 0.5, sn), levels=[0], colors=RED, linewidths=1.0)
        bare(ax)
        ax.set_title(("(a)  " if k == 0 else "") + f"t = {t}", fontsize=7.5)

    ax = fig.add_subplot(bot[0, 0])
    free = read_csv(os.path.join(FGB, "runs", "gb_only_df0.1", "gb_position.csv"))
    cpl = read_csv(os.path.join(FGB, "runs", "coupled_df0.1", "gb_bubble.csv"))
    ax.fill_between(cpl["time"], cpl["bub_zbot_mean"], cpl["bub_ztop_mean"],
                    color=BAND, lw=0, label="bubble layer")
    ok = free["gb_z_mean"] <= 63
    ax.plot(free["time"][ok], free["gb_z_mean"][ok], "--", color=RED, label="GB without bubbles")
    g = cpl["cols_real_gb"] > 0
    ax.plot(cpl["time"][g], cpl["gb_real_mean"][g], "-", color=RED, label="GB with 33 bubbles")
    ax.plot(cpl["time"], 64 - cpl["vol_II"] / 128**2, "-", color=INK2, lw=1.0,
            label="edge of grain II (after GB is gone)")
    ax.axvline(100, color=MUTED, lw=0.7, ls=":")
    ax.text(104, 3, "coverage = 1: GB replaced\nby a gas sheet", fontsize=6.5, color=INK2, va="bottom")
    ax.set_xlim(0, 352); ax.set_ylim(0, 72)
    ax.set_xlabel("time"); ax.set_ylabel("height z (cells)")
    ax.set_title("(b)  GB height, $\\Delta f$ = 0.1")
    ax.legend(loc="upper left", fontsize=6.3, handlelength=1.8)

    ax = fig.add_subplot(bot[0, 1])
    for fixed, moving, col, lab in (("bubbles_only", "coupled_df0.1", AQUA, "$\\theta$ = 120°, N = 33"),
                                    ("theta128_bubbles_only", "theta128_coupled_df0.1", ORANGE, "$\\theta$ = 128°, N = 33"),
                                    ("N22_bubbles_only", "N22_coupled_df0.1", VIOLET, "$\\theta$ = 120°, N = 22")):
        for run, ls, w in ((fixed, "-", 1.4), (moving, "--", 1.1)):
            m = read_csv(os.path.join(FGB, "runs", run, "bubble_metrics.csv"))
            ax.plot(m["coverage"], m["largest_area_fraction"], ls, color=col, lw=w,
                    label=lab if ls == "-" else None)
    h, l = ax.get_legend_handles_labels()
    h += [plt.Line2D([], [], color=INK2, ls="-"), plt.Line2D([], [], color=INK2, ls="--")]
    l += ["GB fixed", "GB moving"]
    ax.legend(h, l, loc="upper left", fontsize=6.3, handlelength=1.8)
    ax.set_xlim(0.3, 0.85); ax.set_ylim(0, 1.04)
    ax.set_xlabel("coverage"); ax.set_ylabel("$A_1/A_{tot}$")
    ax.set_title("(c)  GB fixed vs GB moving")
    save(fig, "f6_coupled.pdf")


# ========================================== f7: path shift vs GB motion
def fig_shift():
    js = json.load(open(os.path.join(FGB, "results", "robustness_summary.json")))
    fig, ax = plt.subplots(figsize=(4.0, 2.45))
    nf = js["noise_floor"]["path_shift"]["max_abs"]
    ax.axhspan(-nf, nf, color=GRID, lw=0)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.text(2.95, nf + 0.001, "run-to-run noise", fontsize=6.3, color=INK2, ha="right", va="bottom")
    lab = {"ref": "reference", "theta128": "$\\theta$ 128°", "N22": "N 22",
           "theta104": "$\\theta$ 104°, $\\Delta f$ 0.3", "df05": "$\\Delta f$ 0.5", "L10": "$L_{GB}$ 10"}
    off = {"ref": (6, 4), "theta128": (6, -9), "N22": (6, 4), "theta104": (6, 3),
           "df05": (6, 2), "L10": (-6, 5)}
    for k in ("ref", "theta128", "N22", "theta104", "df05", "L10"):
        v = js[k]
        x, ps = v["gb_displacement_in_window"], v["path_shift"]
        ax.plot([x, x], [ps["p25"], ps["p75"]], color=INK2, lw=1.2, alpha=0.6)
        ax.plot(x, ps["median"], "o", ms=6, mfc="white" if k == "df05" else BLUE,
                mec=BLUE, mew=1.3, zorder=5)
        dx, dy = off[k]
        ax.annotate(lab[k], (x, ps["median"]), xytext=(dx, dy), textcoords="offset points",
                    fontsize=6.5, color=INK, ha="right" if dx < 0 else "left")
    ax.set_xlim(0, 3.0); ax.set_ylim(-0.035, 0.04)
    ax.set_xlabel("how far the GB moved during percolation (cells)")
    ax.set_ylabel("shift of the curve (coverage)")
    save(fig, "f7_shift.pdf")


# ================================== f8: closed GB and polycrystal outlook
def fig_poly():
    fig = plt.figure(figsize=(6.5, 2.4))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.35, 0.8, 0.8], wspace=0.3)
    ax = fig.add_subplot(gs[0])
    for df, col in (("0", BLUE), ("0.1", RED)):
        nb = read_csv(os.path.join(FGB, "runs", f"sphere_nobub_df{df}", "sphere.csv"))
        wb = read_csv(os.path.join(FGB, "runs", f"sphere_bub_df{df}", "sphere.csv"))
        k = nb["R_vol"] > 0.5
        ax.plot(nb["time"][k], nb["R_vol"][k], "--", color=col, lw=1.1)
        k = wb["R_vol"] > 0.5
        ax.plot(wb["time"][k], wb["R_vol"][k], "-", color=col, label=f"$\\Delta f$ = {df}")
        i = np.argmax(wb["coverage"] >= 0.9999)
        ax.plot(wb["time"][i], wb["R_vol"][i], "o", mfc="white", mec=col, mew=1.2, ms=5, zorder=5)
    h, l = ax.get_legend_handles_labels()
    h += [plt.Line2D([], [], color=INK2, ls="-"), plt.Line2D([], [], color=INK2, ls="--"),
          plt.Line2D([], [], ls="", marker="o", mfc="white", mec=INK2, mew=1.2, ms=5)]
    l += ["with 25 bubbles", "no bubbles", "gas shell closes"]
    ax.legend(h, l, loc="upper right", fontsize=6.3, handlelength=1.8)
    ax.set_xlim(0, 240); ax.set_ylim(0, 34)
    ax.set_xlabel("time"); ax.set_ylabel("grain radius (cells)")
    ax.set_title("(a)  shrinking round grain")

    d = os.path.join(ROOT, "GB2", "fig2_out")
    for k, (s, lab) in enumerate((("00000", "(b)  2-D polycrystal, start"),
                                  ("05000", "(c)  after 5000 steps"))):
        ax = fig.add_subplot(gs[k + 1])
        phi = np.loadtxt(os.path.join(d, f"phi_{s}.csv"))
        g = np.loadtxt(os.path.join(d, f"grain_{s}.csv"))
        rgb = plt.get_cmap("tab20")((g % 20) / 20.0)[..., :3]
        rgb = 0.55 * rgb + 0.45                      # pastel grains
        shade = np.clip((phi - 0.55) / 0.45, 0, 1)[..., None]
        ax.imshow(rgb * shade + (1 - shade) * 0.1, origin="lower", interpolation="bilinear")
        bare(ax)
        ax.set_title(lab, fontsize=7.5)
    save(fig, "f8_poly.pdf")


if __name__ == "__main__":
    fig_single()
    fig_a1()
    fig_morphology()
    fig_percolation()
    fig_flatgb()
    fig_coupled()
    fig_shift()
    fig_poly()
