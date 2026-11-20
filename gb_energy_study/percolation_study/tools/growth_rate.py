"""Is coverage vs time a straight line (constant growth rate), and is the rate / acceleration the
same for two contact angles?  Reads only runs/*/bubble_metrics.csv (coverage every 0.1).

Two angles compared: theta = 104 deg (bGB 3.0) and 148 deg (bGB 1.2); all other parameters and the
initial condition are identical.  Fit window = independent-growth stage:
  start t = 1.0  (excludes the first-unit relaxation of the sharp initial interfaces)
  end   = last output before the first merger (after that, coverage jumps when necks fill).
Measured per run:
  straight-line fit C = a + b t        -> rate b, R^2, max |residual|
  quadratic fit C = a + b t + c t^2   -> mean acceleration 2c
  local rate dC/dt = slope of a linear fit over a sliding 1.0-wide window (11 outputs)
  local acceleration = slope of the local-rate curve over a sliding 2.0-wide window
Outputs: results/growth_rate/{growth_rate.png, growth_fits.csv, local_rate_*.csv}
"""
import os, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f"{ROOT}/results/growth_rate"; os.makedirs(OUT, exist_ok=True)
CASES = [("3.0", 103.9, "#eda100"), ("1.2", 148.0, "#2a78d6")]   # same colours as the other figures
T0 = 1.0


def load(N, b):
    R = list(csv.DictReader(open(f"{ROOT}/runs/N{N}_bGB{b}/bubble_metrics.csv")))
    g = lambda k: np.array([float(r[k]) for r in R])
    return g("time_dimensionless"), g("coverage"), g("number_of_projected_bubbles").astype(int)


def local_slope(t, y, width):
    h = int(round(width / 0.1 / 2)); s = np.full(len(t), np.nan)
    for i in range(h, len(t) - h):
        s[i] = np.polyfit(t[i - h:i + h + 1], y[i - h:i + h + 1], 1)[0]
    return s


def fits(t, C, lo, hi):
    m = (t >= lo - 1e-9) & (t <= hi + 1e-9); tt, cc = t[m], C[m]
    b, a = np.polyfit(tt, cc, 1); res = cc - (a + b * tt)
    r2 = 1 - (res ** 2).sum() / ((cc - cc.mean()) ** 2).sum()
    q = np.polyfit(tt, cc, 2)
    return dict(t_lo=lo, t_hi=hi, n_pts=int(m.sum()), rate=b, intercept=a, R2=r2,
                max_abs_resid=np.abs(res).max(), accel_quadfit=2 * q[0]), (tt, cc, a + b * tt, res)


rows = []; plotdata = {}
for N in (33, 25):
    for b, th, col in CASES:
        t, C, n = load(N, b)
        i_m = int(np.argmax(n < n[0])); t_m = t[i_m]; hi = t[i_m - 1]
        rate = local_slope(t, C, 1.0); acc = local_slope(t, rate, 2.0)
        f, pd = fits(t, C, T0, hi)
        k_lo = np.nanargmin(np.abs(t - (T0 + 0.5))); k_hi = np.nanargmin(np.abs(t - (hi - 0.5)))
        f.update(N=N, bGB=b, theta=th, t_first_merger=t_m, rate_start=rate[k_lo], rate_end=rate[k_hi],
                 rate_change_pct=100 * (rate[k_hi] - rate[k_lo]) / rate[k_lo])
        # common window for both angles (both still before any merger)
        fc, _ = fits(t, C, T0, 7.9)
        f.update(rate_common_1_7p9=fc["rate"], accel_common_1_7p9=fc["accel_quadfit"], R2_common=fc["R2"])
        rows.append(f); plotdata[(N, b)] = (t, C, n, rate, acc, pd, f)
        with open(f"{OUT}/local_rate_N{N}_bGB{b}.csv", "w") as fh:
            w = csv.writer(fh); w.writerow(["time", "coverage", "bubbles", "dC_dt_1.0window", "d2C_dt2_2.0window"])
            for z in zip(t, C, n, rate, acc): w.writerow([f"{z[0]:.1f}", f"{z[1]:.6f}", z[2], f"{z[3]:.6g}", f"{z[4]:.6g}"])

keys = ["N", "bGB", "theta", "t_first_merger", "t_lo", "t_hi", "n_pts", "rate", "intercept", "R2", "max_abs_resid",
        "accel_quadfit", "rate_start", "rate_end", "rate_change_pct", "rate_common_1_7p9", "accel_common_1_7p9", "R2_common"]
with open(f"{OUT}/growth_fits.csv", "w") as fh:
    w = csv.DictWriter(fh, keys); w.writeheader()
    for r in rows: w.writerow({k: (f"{r[k]:.6g}" if isinstance(r[k], float) else r[k]) for k in keys})
for r in rows:
    print({k: (round(r[k], 6) if isinstance(r[k], float) else r[k]) for k in keys})

# window robustness: same time window for both angles, three start times
with open(f"{OUT}/window_robustness.csv", "w") as fh:
    w = csv.writer(fh); w.writerow(["N", "window", "rate_104", "rate_148", "ratio", "accel_104", "accel_148", "R2_104", "R2_148"])
    for N in (33, 25):
        for lo, hi in ((1, 7.8), (2, 7.8), (3, 7.8)):
            fa, _ = fits(*load(N, "3.0")[:2], lo, hi); fb, _ = fits(*load(N, "1.2")[:2], lo, hi)
            w.writerow([N, f"{lo}-{hi}", f"{fa['rate']:.5f}", f"{fb['rate']:.5f}", f"{fa['rate']/fb['rate']:.3f}",
                        f"{fa['accel_quadfit']:.5f}", f"{fb['accel_quadfit']:.5f}", f"{fa['R2']:.5f}", f"{fb['R2']:.5f}"])
            print("window", N, lo, hi, round(fa["rate"] / fb["rate"], 3), round(fa["accel_quadfit"], 5), round(fb["accel_quadfit"], 5))

# rate at equal coverage (before first merger, t >= 1.5 so the 1.0-wide window is past the relaxation)
eq = {}
with open(f"{OUT}/rate_at_equal_coverage.csv", "w") as fh:
    w = csv.writer(fh); w.writerow(["N", "coverage", "t_104", "rate_104", "t_148", "rate_148", "rate_ratio"])
    for N in (33, 25):
        cur = {}
        for b, th, col in CASES:
            t, C, n, rate, acc, pd, f = plotdata[(N, b)]
            m = (t >= 1.5) & (t < f["t_first_merger"] - 0.45)
            cur[b] = (C[m], rate[m], t[m])
        grid = np.arange(0.18, min(cur["3.0"][0].max(), cur["1.2"][0].max()) - 1e-9, 0.02)
        ratios = []
        for cv in grid:
            ra = np.interp(cv, cur["3.0"][0], cur["3.0"][1]); rb = np.interp(cv, cur["1.2"][0], cur["1.2"][1])
            ta = np.interp(cv, cur["3.0"][0], cur["3.0"][2]); tb = np.interp(cv, cur["1.2"][0], cur["1.2"][2])
            ratios.append((cv, ra / rb)); w.writerow([N, f"{cv:.2f}", f"{ta:.2f}", f"{ra:.5f}", f"{tb:.2f}", f"{rb:.5f}", f"{ra/rb:.3f}"])
        k = np.mean([r for cv, r in ratios if 0.199 <= cv <= 0.301])   # mean ratio over coverage 0.20-0.30
        eq[N] = (cur, k); print("equal coverage", N, [(round(a, 2), round(b, 2)) for a, b in ratios], "mean 0.20-0.30:", round(k, 2))

# ---------------- figure (N = 33 on top, N = 25 below) ----------------
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e4e3df", "grid.linewidth": 0.6})
fig, ax = plt.subplots(2, 4, figsize=(19, 8.8))
for j, N in enumerate((33, 25)):
    for b, th, col in CASES:
        t, C, n, rate, acc, (tt, cc, line, res), f = plotdata[(N, b)]
        lbl = f"θ = {th:.0f}°"
        m = t <= f["t_first_merger"] + 6
        ax[j, 0].plot(t[m], C[m], color=col, lw=2.2, label=f"{lbl}: simulation")
        ax[j, 0].plot(tt, line, color="#0b0b0b", lw=1.2, ls="--",
                      label=f"{lbl}: line, slope {f['rate']:.4f}/time, R² {f['R2']:.4f}")
        ax[j, 0].axvline(f["t_first_merger"], color=col, lw=0.8, ls=":")
        ax[j, 1].plot(tt, res * 1e3, color=col, lw=2, label=lbl)
        ax[j, 2].plot(t[m], rate[m], color=col, lw=2, label=f"{lbl}: local dC/dt")
        ax[j, 2].hlines(f["rate"], tt[0], tt[-1], color=col, lw=1.2, ls="--")
        ax[j, 2].axvline(f["t_first_merger"], color=col, lw=0.8, ls=":")
    cur, k = eq[N]
    ax[j, 3].plot(cur["3.0"][0], cur["3.0"][1], color=CASES[0][2], lw=2, label="θ = 104°")
    ax[j, 3].plot(cur["1.2"][0], cur["1.2"][1], color=CASES[1][2], lw=2, label="θ = 148°")
    ax[j, 3].plot(cur["1.2"][0], k * cur["1.2"][1], color=CASES[1][2], lw=1.3, ls="--",
                  label=f"θ = 148° × {k:.2f} (mean ratio, coverage 0.20–0.30)")
    ax[j, 0].set(title=f"N = {N}: coverage vs time", xlabel="time (dimensionless)", ylabel="coverage")
    ax[j, 0].legend(frameon=False, fontsize=7.5, loc="upper left")
    ax[j, 1].axhline(0, color="#52514e", lw=0.8)
    ax[j, 1].set(title=f"N = {N}: simulation − straight line", xlabel="time (dimensionless)",
                 ylabel="residual (×10⁻³ coverage)")
    ax[j, 1].legend(frameon=False, fontsize=8)
    ax[j, 2].set(title=f"N = {N}: growth rate vs time", xlabel="time (dimensionless)",
                 ylabel="dC/dt (per time unit)", ylim=(0, None))
    ax[j, 2].legend(frameon=False, fontsize=8, loc="upper right")
    ax[j, 3].set(title=f"N = {N}: growth rate vs coverage (before 1st merger)", xlabel="coverage",
                 ylabel="dC/dt (per time unit)", ylim=(0, None))
    ax[j, 3].legend(frameon=False, fontsize=7.5, loc="lower left")
fig.text(0.5, 0.004, "Dashed black = straight-line fit from t = 1 (initial interface relaxation excluded) to the last output before the first merger "
         "(dotted vertical line). A straight line would give a flat residual and a flat dC/dt.", ha="center", fontsize=9, color="#52514e")
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(f"{OUT}/growth_rate.png", dpi=140); plt.close(fig)
print("written", OUT)
