"""
plot_gb_energy_study.py

Combines the four GB-energy cases (and the optional Ms = Mgb = 0 runs) into
one CSV and draws figures A-F.

Usage (from the study root):
    python3 plot_gb_energy_study.py

Expects:
    cases/bGB_<b>/bubble_metrics.csv          main runs
    cases_mob0/bGB_<b>/bubble_metrics.csv     mobility test (optional)
    results/case_properties.json              measured gamma_GB and theta
                                              (written by measure_contact_angles.py)
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

CASES = ["1.2", "1.5", "2.0", "3.0"]
COLORS = {"1.2": "#1f4e9c", "1.5": "#2a9d8f", "2.0": "#e09f3e", "3.0": "#c0392b"}
STYLES = {"1.2": "-", "1.5": "--", "2.0": "-.", "3.0": ":"}
OUT = "results"
os.makedirs(OUT, exist_ok=True)

props = json.load(open(os.path.join(OUT, "case_properties.json")))


def label(b):
    p = props[b]
    return (f"bGB={b}: $\\gamma_{{GB}}$={p['gGB']:.2f}, "
            f"$\\theta_{{meas}}$={p['theta_meas']:.0f}°")


def load(root):
    frames = []
    for b in CASES:
        fn = os.path.join(root, f"bGB_{b}", "bubble_metrics.csv")
        if not os.path.exists(fn):
            continue
        m = pd.read_csv(fn)
        m.insert(0, "bGB", float(b))
        m.insert(1, "gamma_GB", props[b]["gGB"])
        m.insert(2, "theta_measured_deg", props[b]["theta_meas"])
        frames.append(m)
    return pd.concat(frames, ignore_index=True) if frames else None


main = load("cases")
main.insert(0, "mobility_set", "baseline")
mob0 = load("cases_mob0")
combined = main
if mob0 is not None:
    mob0.insert(0, "mobility_set", "Ms=Mgb=0")
    combined = pd.concat([main, mob0], ignore_index=True)
combined.to_csv(os.path.join(OUT, "combined_metrics.csv"), index=False)


def per_case(ax, df, y, x="time_dimensionless", ls="-", lab=True, **kw):
    for b in CASES:
        d = df[df.bGB == float(b)]
        if len(d):
            ax.plot(d[x], d[y], STYLES[b] if ls == "-" else ls, color=COLORS[b], label=label(b) if lab else None, **kw)


# ---------------- A: main plot ----------------
fig, ax = plt.subplots(figsize=(7.5, 6))
per_case(ax, main, "largest_area_fraction", x="coverage", lw=2)
ax.axhline(1.0, color="grey", ls="--", lw=0.8)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.05)
ax.set_xlabel("Coverage")
ax.set_ylabel(r"$A_{largest}/A_{total}$")
ax.set_title("A. Largest connected bubble fraction vs coverage")
ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "A_largest_fraction_vs_coverage.png"), dpi=200)

# ---------------- B-E: time series ----------------
fig, axs = plt.subplots(2, 2, figsize=(12, 9))
per_case(axs[0, 0], main, "coverage", lw=1.8)
axs[0, 0].set_ylabel("Coverage"); axs[0, 0].set_title("B. Coverage vs time")
per_case(axs[0, 1], main, "largest_area_fraction", lw=1.8, lab=False)
axs[0, 1].set_ylabel(r"$A_{largest}/A_{total}$"); axs[0, 1].set_title("C. Largest connected fraction vs time")
per_case(axs[1, 0], main, "number_of_projected_bubbles", lw=1.8, lab=False)
axs[1, 0].set_ylabel("Connected projected bubbles"); axs[1, 0].set_title("D. Number of bubbles vs time")
per_case(axs[1, 1], main, "mass_drift", lw=1.8, lab=False)
axs[1, 1].set_ylabel("(gas - gas$_0$)/gas$_0$"); axs[1, 1].set_title("E. Relative gas-mass drift")
axs[1, 1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
for a in axs.flat:
    a.set_xlabel("Time (dimensionless)"); a.grid(alpha=0.3)
axs[0, 0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "BCDE_time_series.png"), dpi=200)

# ---------------- F: contact angle vs GB energy ----------------
fig, ax = plt.subplots(figsize=(6.5, 5))
g = [props[b]["gGB"] for b in CASES]
tm = [props[b]["theta_meas"] for b in CASES]
te = [props[b]["theta_meas_spread"] for b in CASES]
ty = [props[b]["theta_young"] for b in CASES]
ax.errorbar(g, tm, yerr=te, fmt="o", ms=7, capsize=4, color="k", label="measured (2-D equilibrated lens)")
import numpy as np
gg = np.linspace(0.3, 1.5, 200)
gpg = props["gamma_PG"]
ax.plot(gg, np.degrees(2 * np.arccos(gg / (2 * gpg))), "-", color="grey",
        label=r"Young: $\gamma_{GB}=2\gamma_{PG}\cos(\theta/2)$, $\gamma_{PG}$=%.3f" % gpg)
ax.plot(g, ty, "x", color="grey")
for b, x, y in zip(CASES, g, tm):
    ax.annotate(f"bGB={b}", (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8)
ax.set_xlabel(r"$\gamma_{GB}$ (1-D equilibrium excess energy, model units)")
ax.set_ylabel(r"Dihedral angle through bubble $\theta$ (deg)")
ax.set_title("F. Contact angle vs grain-boundary energy")
ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "F_angle_vs_gb_energy.png"), dpi=200)

# ---------------- mobility sensitivity ----------------
if mob0 is not None:
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.5))
    per_case(axs[0], main, "largest_area_fraction", x="coverage", lw=2)
    per_case(axs[0], mob0.iloc[::5], "largest_area_fraction", x="coverage", ls="o", lab=False, ms=2.5)
    axs[0].set_xlabel("Coverage"); axs[0].set_ylabel(r"$A_{largest}/A_{total}$")
    axs[0].set_xlim(0, 1); axs[0].set_ylim(0, 1.05)
    axs[0].set_title("Lines: baseline mobility; dots: Ms = Mgb = 0")
    per_case(axs[1], main, "coverage", lw=2, lab=False)
    per_case(axs[1], mob0.iloc[::5], "coverage", ls="o", lab=False, ms=2.5)
    axs[1].set_xlabel("Time (dimensionless)"); axs[1].set_ylabel("Coverage")
    for a in axs:
        a.grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "G_mobility_sensitivity.png"), dpi=200)

print("wrote", OUT)
