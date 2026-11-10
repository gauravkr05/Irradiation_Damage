"""Figures for the percolation study (reads results/events_summary.csv + runs/*/bubble_metrics.csv)."""
import os, csv, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analyze import load_metrics, ROOT, OUT, THETA

COL = {"1.2": "#2a78d6", "1.5": "#eb6834", "2.0": "#1baf7a", "3.0": "#eda100"}
MK = {"1.2": "o", "1.5": "s", "2.0": "^", "3.0": "D"}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6, "lines.linewidth": 2})
S = {r["case"]: r for r in csv.DictReader(open(f"{OUT}/events_summary.csv"))}
f_ = lambda v: float(v) if v not in ("", "None", None) else None
lab = lambda b: f"θ = {THETA[b]:.0f}° (bGB {b})"


def series(N):
    for b in ("1.2", "1.5", "2.0", "3.0"):
        d = f"{ROOT}/runs/N{N}_bGB{b}"
        if os.path.exists(f"{d}/bubble_metrics.csv"): yield b, load_metrics(d), S.get(f"N{N}_bGB{b}", {})


# Fig 1: bubble count vs coverage, and coverage vs time
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
for j, N in enumerate((33, 25)):
    for b, m, s in series(N):
        t, n, cov = m["time_dimensionless"], m["number_of_projected_bubbles"], m["coverage"]
        ax[0, j].plot(cov, n, color=COL[b], label=lab(b))
        ax[1, j].plot(t, cov, color=COL[b], label=lab(b))
        for key, fill in (("first_merger", "white"), ("single", COL[b])):
            if s and f_(s[f"cov_{key}"]) is not None:
                ax[0, j].plot(f_(s[f"cov_{key}"]), N if key == "first_merger" else 1, MK[b], ms=8,
                              mfc=fill, mec=COL[b], mew=2, zorder=5)
                ax[1, j].plot(f_(s[f"t_{key}"]), f_(s[f"cov_{key}"]), MK[b], ms=8, mfc=fill, mec=COL[b], mew=2, zorder=5)
    ax[0, j].set(title=f"N = {N}: projected bubble count vs coverage", xlabel="coverage", ylabel="number of projected bubbles")
    ax[1, j].set(title=f"N = {N}: coverage vs time", xlabel="time (dimensionless)", ylabel="coverage")
    ax[0, j].legend(frameon=False, fontsize=8)
fig.text(0.5, 0.005, "open marker = first merger (percolation begin); filled marker = single bubble left (percolation end)",
         ha="center", fontsize=9, color="#52514e")
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(f"{OUT}/fig1_count_and_coverage.png", dpi=150); plt.close(fig)

# Fig 2: coverage and time at each event vs contact angle
EV = [("first_merger", "first merger"), ("half", "A1/Atot > 0.5"), ("single", "single bubble left")]
fig, ax = plt.subplots(2, 3, figsize=(13, 7.5))
for k, (key, name) in enumerate(EV):
    for N, ls, mk in ((33, "-", "o"), (25, "--", "s")):
        th, cv, cp, tt = [], [], [], []
        for b in ("1.2", "1.5", "2.0", "3.0"):
            s = S.get(f"N{N}_bGB{b}")
            if s and f_(s[f"cov_{key}"]) is not None:
                th.append(THETA[b]); cv.append(f_(s[f"cov_{key}"])); cp.append(f_(s[f"covprev_{key}"])); tt.append(f_(s[f"t_{key}"]))
        if th:
            ax[0, k].errorbar(th, cv, yerr=[np.array(cv) - np.array(cp), 0 * np.array(cv)], fmt=mk + ls,
                              color="#2a78d6" if N == 33 else "#eb6834", ms=8, capsize=3, label=f"N = {N}")
            ax[1, k].plot(th, tt, mk + ls, color="#2a78d6" if N == 33 else "#eb6834", ms=8, label=f"N = {N}")
    ax[0, k].set(title=f"coverage at {name}", xlabel="measured contact angle θ (°)", ylabel="coverage")
    ax[1, k].set(title=f"time of {name}", xlabel="measured contact angle θ (°)", ylabel="time (dimensionless)")
    ax[0, k].legend(frameon=False); ax[1, k].legend(frameon=False)
fig.text(0.5, 0.005, "error bar = coverage at the previous 0.1 output (event lies inside this bracket); θ itself carries ±6°",
         ha="center", fontsize=9, color="#52514e")
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(f"{OUT}/fig2_events_vs_angle.png", dpi=150); plt.close(fig)

# Fig 3: last thickness map of every run (what the 'single bubble' looks like)
fig, ax = plt.subplots(2, 4, figsize=(14, 7.2))
for j, N in enumerate((33, 25)):
    for i, b in enumerate(("1.2", "1.5", "2.0", "3.0")):
        fs = sorted(glob.glob(f"{ROOT}/runs/N{N}_bGB{b}/thick/t*.npy"))
        a = ax[j, i]; a.set_xticks([]); a.set_yticks([]); a.grid(False)
        if not fs: continue
        s = S.get(f"N{N}_bGB{b}", {}); ts = f_(s.get("t_single"))
        pick = [f for f in fs if ts is not None and float(os.path.basename(f)[1:-4]) >= ts]
        f = pick[0] if pick else fs[-1]
        th = np.load(f)
        a.imshow(th, cmap="Blues", origin="lower", vmin=0, vmax=max(1, th.max()))
        a.contour(th >= 3, levels=[0.5], colors="#0b0b0b", linewidths=0.6)
        a.set_title(f"N={N}, θ={THETA[b]:.0f}°, t={float(os.path.basename(f)[1:-4]):.0f}", fontsize=9)
fig.suptitle("Column thickness (cells with c ≥ 0.5) at the first map with a single bubble (or the last map); black = covered boundary")
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_final_maps.png", dpi=130); plt.close(fig)
print("figures written to", OUT)

# Fig 4: late-stage behaviour - min gap between islands, coverage and porosity vs time
G = list(csv.DictReader(open(f"{OUT}/gaps.csv")))
fig, ax = plt.subplots(2, 3, figsize=(15, 7.5))
for j, N in enumerate((33, 25)):
    for b, m, s in series(N):
        g = [(float(r["t"]), float(r["min_gap"])) for r in G if r["case"] == f"N{N}_bGB{b}" and r["min_gap"]]
        if g: ax[j, 0].plot(*zip(*g), color=COL[b], lw=1.5, label=lab(b))
        ax[j, 1].plot(m["time_dimensionless"], m["coverage"], color=COL[b], label=lab(b))
        ax[j, 2].plot(m["time_dimensionless"], m["porosity"], color=COL[b], label=lab(b))
        if s and f_(s["t_single"]) is not None:
            for a in ax[j, 1:]: a.axvline(f_(s["t_single"]), color=COL[b], lw=0.8, ls=":")
    ax[j, 0].set(title=f"N = {N}: smallest gap between two islands", xlabel="time", ylabel="empty cells between islands")
    ax[j, 1].set(title=f"N = {N}: coverage (dotted = single bubble)", xlabel="time", ylabel="coverage")
    ax[j, 2].set(title=f"N = {N}: porosity (gas volume fraction)", xlabel="time", ylabel="porosity")
    ax[j, 0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{OUT}/fig4_late_stage.png", dpi=140); plt.close(fig)
print("fig4 written")
