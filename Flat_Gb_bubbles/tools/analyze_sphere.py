"""analyze_sphere.py - shrinking spherical grain (R0 = 32, 96^3, theta = 120 deg), with and without bubbles.

runs/sphere_nobub_df{0,0.025,0.05,0.1} : no bubbles          (c = 0)
runs/sphere_bub_df{0,0.025,0.05,0.1}   : 25 bubbles on the GB (cInit = 0.18), same IC for all dfBias
Writes results/sphere_fig1_radius_vs_time.png, results/sphere_fig2_coverage.png,
results/sphere_fig3_slices.png, results/sphere_summary.json
"""
import csv, glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
NX = NY = NZ = 96
R0 = 32.0
GAMMA = 1.0724          # GB energy for bGB = 2 (Flat_Gb control run)
L, KAPPA = 1.0, 2.0     # LGB, jGB
DFS = ["0", "0.025", "0.05", "0.1"]
COL = {"0": "#1f77b4", "0.025": "#2ca02c", "0.05": "#ff7f0e", "0.1": "#d62728"}


def load(fn):
    if not os.path.exists(fn):
        return None
    rows = list(csv.DictReader(open(fn)))
    if not rows:
        return None
    return {k: np.array([float(r[k]) if r[k] not in ("", "nan") else np.nan for r in rows]) for k in rows[0]}


def sharp(df, tmax, r0=R0):
    """sharp-interface Allen-Cahn: dR/dt = -L*kappa*(2/R + df/gamma)"""
    t, r, ts, rs = 0.0, r0, [0.0], [r0]
    dt = 0.01
    while r > 0.3 and t < tmax:
        k1 = -L * KAPPA * (2 / r + df / GAMMA)
        rm = r + 0.5 * dt * k1
        if rm <= 0: break
        k2 = -L * KAPPA * (2 / rm + df / GAMMA)
        r += dt * k2; t += dt
        ts.append(t); rs.append(r)
    return np.array(ts), np.array(rs)


def first(t, cond):
    i = np.flatnonzero(cond)
    return float(t[i[0]]) if len(i) else None


nob = {d: load(R(f"runs/sphere_nobub_df{d}/sphere.csv")) for d in DFS}
bub = {d: load(R(f"runs/sphere_bub_df{d}/sphere.csv")) for d in DFS}
S = {}
for name, D in (("nobub", nob), ("bub", bub)):
    for d, x in D.items():
        if x is None: continue
        t = x["time"]
        tv = first(t, x["V_inner"] == 0)
        ts, rs = sharp(float(d), 1e4)
        e = {"t_end": float(t[-1]), "t_vanish": tv, "t_vanish_sharp": float(ts[-1]) if rs[-1] < 1 else None,
             "R_gb_final": float(x["r_gb_mean"][-1]), "R_vol_final": float(x["R_vol"][-1]),
             "t_R_gb_below_16": first(t, x["r_gb_mean"] < 16), "t_R_vol_below_16": first(t, x["R_vol"] < 16),
             "max_abs_mass_drift": float(np.nanmax(abs(x["mass_drift"]))),
             "free_energy_increases": int((np.diff(x["free_energy"]) > 1e-9 * abs(x["free_energy"][1:])).sum()),
             "nonfinite": int(x["nonfinite"].sum()), "max_centroid_offset": float(np.nanmax(x["centroid_offset"]))}
        if name == "bub":
            e.update({"coverage_max": float(np.nanmax(x["coverage"])), "t_cov_0.8": first(t, x["coverage"] >= 0.8),
                      "t_cov_0.9": first(t, x["coverage"] >= 0.9), "t_single_patch": first(t, (x["n_patches"] == 1) & (t > 0)),
                      "t_cov_1": first(t, x["coverage"] >= 0.9999), "gas_cells_final": float(x["gas_cells"][-1])})
        S[f"{name}_df{d}"] = e
json.dump(S, open(R("results/sphere_summary.json"), "w"), indent=1)
for k, v in S.items():
    print(k, {a: (round(b, 3) if isinstance(b, float) else b) for a, b in v.items()})

# ---------------- Fig 1: radius vs time
have_bub = any(v is not None for v in bub.values())
fig, axs = plt.subplots(1, 2 if have_bub else 1, figsize=(13 if have_bub else 7.5, 5.2), squeeze=False)
ax = axs[0, 0]
for d in DFS:
    x = nob[d]
    if x is None: continue
    ax.plot(x["time"], x["r_gb_mean"], color=COL[d], lw=2, label=f"$\\Delta f$ = {d}: simulation (GB radius, mean over rays)")
    ts, rs = sharp(float(d), x["time"][-1] + 20)
    ax.plot(ts, rs, color=COL[d], lw=1, ls="--", alpha=0.8)
ax.plot([], [], color="k", lw=1, ls="--", label="sharp-interface: $dR/dt = -L\\kappa\\,(2/R + \\Delta f/\\gamma_{GB})$")
ax.set_xlabel("time (dimensionless)"); ax.set_ylabel("radius of the inner grain (cells)")
ax.set_ylim(0, 34); ax.set_xlim(left=0)
ax.set_title("No bubbles: inner spherical grain shrinking ($R_0$ = 32, $\\theta$ = 120°)", fontsize=10.5)
ax.legend(fontsize=8); ax.grid(alpha=0.3)
if have_bub:
    ax = axs[0, 1]
    for d in DFS:
        x = nob[d]
        if x is not None:
            ax.plot(x["time"], x["R_vol"], color=COL[d], lw=1.2, ls="--", alpha=0.9)
        x = bub[d]
        if x is None: continue
        tv = S[f"bub_df{d}"]["t_vanish"]; tn = S[f"nobub_df{d}"]["t_vanish"] if nob[d] is not None else float("nan")
        ax.plot(x["time"], x["R_vol"], color=COL[d], lw=2.2, label=f"$\\Delta f$ = {d}: with 25 bubbles, vanishes t = {tv:.0f} (no bubbles: {tn:.0f})")
        ax.plot(x["time"], x["r_gb_mean"], color=COL[d], lw=1, ls=":", alpha=0.9)
        tc = S[f"bub_df{d}"]["t_cov_1"]
        if tc is not None:
            ax.plot([tc], [x["R_vol"][np.argmin(abs(x["time"] - tc))]], "o", color=COL[d], ms=6, mec="k", mew=0.6)
    ax.plot([], [], color="k", lw=1.2, ls="--", label="no bubbles (same measure)")
    ax.plot([], [], color="k", lw=1, ls=":", label="with bubbles: GB radius on the rays where GB necks remain")
    ax.plot([], [], "o", color="w", mec="k", label="gas shell closes (coverage = 1, no GB left)")
    ax.set_xlabel("time (dimensionless)"); ax.set_ylabel("equivalent radius of the inner grain  $(3V_{inner}/4\\pi)^{1/3}$ (cells)")
    ax.set_ylim(0, 34); ax.set_xlim(left=0)
    ax.set_title("With 25 bubbles on the GB (cInit = 0.18) vs without", fontsize=10.5)
    ax.legend(fontsize=7.5); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(R("results/sphere_fig1_radius_vs_time.png"), dpi=170)

# ---------------- Fig 0: the no-bubble result on its own (first deliverable), plus a CSV of the curves
fig, ax = plt.subplots(figsize=(7.5, 5.2))
tab = {}
for d in DFS:
    x = nob[d]
    if x is None: continue
    tv = S[f"nobub_df{d}"]["t_vanish"]; tsh = S[f"nobub_df{d}"]["t_vanish_sharp"]
    ax.plot(x["time"], x["r_gb_mean"], color=COL[d], lw=2, label=f"$\\Delta f$ = {d}: simulation, grain vanishes at t = {tv:.0f} (sharp: {tsh:.0f})")
    ts, rs = sharp(float(d), x["time"][-1] + 20)
    ax.plot(ts, rs, color=COL[d], lw=1, ls="--", alpha=0.8)
    tab[d] = (x["time"], x["r_gb_mean"], x["R_vol"])
ax.plot([], [], color="k", lw=1, ls="--", label="sharp-interface: $dR/dt = -L\\kappa\\,(2/R + \\Delta f/\\gamma_{GB})$")
ax.set_xlabel("time (dimensionless)"); ax.set_ylabel("radius of the inner grain (cells)")
ax.set_ylim(0, 34); ax.set_xlim(0, 145)
ax.set_title("No bubbles: inner spherical grain ($R_0$ = 32, $\\theta$ = 120°) shrinking inside the outer grain", fontsize=10)
ax.legend(fontsize=8.5); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(R("results/sphere_fig0_nobubbles_radius_vs_time.png"), dpi=170)
with open(R("results/sphere_nobubbles_radius_vs_time.csv"), "w") as f:
    f.write("time," + ",".join(f"R_gb_df{d},R_vol_df{d}" for d in tab) + "\n")
    n = max(len(v[0]) for v in tab.values()); t0 = max(tab.values(), key=lambda v: len(v[0]))[0]
    for i in range(n):
        f.write(f"{t0[i]:.1f}," + ",".join((f"{v[1][i]:.3f},{v[2][i]:.3f}" if i < len(v[0]) else ",") for v in tab.values()) + "\n")

# ---------------- Fig 2: coverage (bubble runs)
if have_bub:
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for d in DFS:
        x = bub[d]
        if x is None: continue
        axs[0].plot(x["time"], x["coverage"], color=COL[d], lw=2, label=f"$\\Delta f$ = {d}")
        axs[1].plot(x["coverage"], x["largest_fraction"], color=COL[d], lw=1.8, label=f"$\\Delta f$ = {d}")
    axs[0].set_xlabel("time"); axs[0].set_ylabel("fraction of the spherical GB covered by gas (solid angle)")
    axs[0].set_ylim(0, 1.02); axs[0].grid(alpha=0.3); axs[0].legend(fontsize=9)
    axs[1].set_xlabel("coverage of the spherical GB"); axs[1].set_ylabel("largest patch / total covered  $A_1/A_{tot}$")
    axs[1].set_ylim(0, 1.04); axs[1].grid(alpha=0.3); axs[1].legend(fontsize=9)
    fig.suptitle("Bubbles on the shrinking spherical GB", fontsize=11)
    fig.tight_layout()
    fig.savefig(R("results/sphere_fig2_coverage.png"), dpi=170)


# ---------------- Fig 3: xz slices
def read_map(fn):
    b = open(fn, "rb").read(); o = 0
    sc = np.frombuffer(b, np.float32, NZ * NX, o).reshape(NZ, NX); o += 4 * NZ * NX
    sn = np.frombuffer(b, np.float32, NZ * NX, o).reshape(NZ, NX)
    return sc, sn


runs = [(f"$\\Delta f$ = {d}, no bubbles", R(f"runs/sphere_nobub_df{d}")) for d in ("0", "0.1")]
if have_bub:
    runs += [(f"$\\Delta f$ = {d}, 25 bubbles", R(f"runs/sphere_bub_df{d}")) for d in ("0", "0.1")]
runs = [(n, p) for n, p in runs if glob.glob(p + "/maps/t*.bin")]
if runs:
    want = [0, 10, 25, 50, 75, 100, 125]
    fig, axs = plt.subplots(len(runs), len(want), figsize=(1.9 * len(want), 2.0 * len(runs) + 0.4), squeeze=False)
    for r, (name, p) in enumerate(runs):
        maps = sorted(glob.glob(p + "/maps/t*.bin")); tm = np.array([float(os.path.basename(f)[1:-4]) for f in maps])
        for k, w in enumerate(want):
            a = axs[r, k]; i = int(np.argmin(abs(tm - w)))
            if abs(tm[i] - w) > 0.6: a.axis("off"); continue
            sc, sn = read_map(maps[i])
            a.imshow(sc, origin="lower", cmap="Greys", vmin=0, vmax=1, extent=(0, NX, 0, NZ))
            a.contour(np.arange(NX) + 0.5, np.arange(NZ) + 0.5, np.ma.masked_where(sc > 0.5, sn), levels=[0], colors="r", linewidths=1.2)
            a.set_xticks([]); a.set_yticks([])
            if r == 0: a.set_title(f"t = {tm[i]:.0f}", fontsize=9)
            if k == 0: a.set_ylabel(name, fontsize=8)
    fig.suptitle("xz slice through the centre: black = gas, red = GB (n1 = n2)", fontsize=10)
    fig.tight_layout()
    fig.savefig(R("results/sphere_fig3_slices.png"), dpi=150)
print("done")
