"""analyze.py - figures and numbers for the coupled GB-migration + bubble-growth study.

Runs (all bGB = 2.0, theta = 120 deg, grid 128x128x64 unless noted, GB and bubble centres at z = 20):
  runs/bubbles_only     : 33 bubbles, dfBias = 0      (GB cannot migrate: independent bubble run)
  runs/gb_only_df0.1    : no bubbles, dfBias = 0.1    (independent GB run; planar, 32x32x64 is exact)
  runs/coupled_df0.1    : 33 bubbles, dfBias = 0.1    (both together)
Writes results/fig1_gb_position.png, results/fig2_largest_fraction_vs_coverage.png,
results/fig3_xz_slices.png and results/summary.json.
"""
import csv, glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
NX, NY, NZ = 128, 128, 64
Z_TOP = NZ - 1          # last cell centre (no-flux wall)


def load(fn):
    rows = list(csv.DictReader(open(fn)))
    return {k: np.array([float(r[k]) if r[k] not in ("", "nan") else np.nan for r in rows]) for k in rows[0]}


bo, bo_c = load(R("runs/bubbles_only/bubble_metrics.csv")), load(R("runs/bubbles_only/gb_bubble.csv"))
cp, cp_c = load(R("runs/coupled_df0.1/bubble_metrics.csv")), load(R("runs/coupled_df0.1/gb_bubble.csv"))
go = load(R("runs/gb_only_df0.1/gb_position.csv"))
t_bo, t_cp, t_go = bo["time_dimensionless"], cp["time_dimensionless"], go["time"]
assert np.allclose(t_cp, cp_c["time"]) and np.allclose(t_bo, bo_c["time"])


def first(t, cond):
    i = np.flatnonzero(cond)
    return float(t[i[0]]) if len(i) else None


S = {}
# ---- independent runs: the timing check
for lev in (0.7, 0.8, 0.85, 0.9):
    S[f"bubbles_only_t_cov_{lev}"] = first(t_bo, bo["coverage"] >= lev)
    S[f"coupled_t_cov_{lev}"] = first(t_cp, cp["coverage"] >= lev)
S["bubbles_only_max_cov"] = float(bo["coverage"].max())
S["bubbles_only_t_max_cov"] = float(t_bo[bo["coverage"].argmax()])
S["bubbles_only_t_single"] = first(t_bo, bo["number_of_projected_bubbles"] == 1)
S["bubbles_only_tend"] = float(t_bo[-1])
gone = go["columns_no_crossing"] == go["columns_with_gb"] + go["columns_no_crossing"]
for zl in (40, 50, 55, 60):
    S[f"gb_only_t_reach_z{zl}"] = first(t_go, go["gb_z_mean"] >= zl)
S["gb_only_t_gb_vanished"] = first(t_go, go["columns_with_gb"] == 0)
S["gb_only_tend"] = float(t_go[-1])
m = (t_go >= 10) & (go["gb_z_mean"] <= 50)
S["gb_only_v_fit_t10_to_z50"] = float(np.polyfit(t_go[m], go["gb_z_mean"][m], 1)[0])
# ---- coupled run
S["coupled_tend"] = float(t_cp[-1])
S["coupled_max_cov"] = float(cp["coverage"].max())
S["coupled_t_single"] = first(t_cp, cp["number_of_projected_bubbles"] == 1)
for zl in (25, 30, 40, 50, 55, 60):
    S[f"coupled_t_gb_reach_z{zl}"] = first(t_cp, cp_c["gb_real_mean"] >= zl)
S["coupled_gb_final"] = float(cp_c["gb_real_mean"][-1])
S["coupled_bubble_z_final"] = float(cp_c["bub_z_mean"][-1])
for k in ("gas", "free_energy"):
    pass
S["max_abs_mass_drift"] = {"bubbles_only": float(np.nanmax(abs(bo["mass_drift"]))),
                           "coupled": float(np.nanmax(abs(cp["mass_drift"])))}
S["free_energy_increases"] = {
    "bubbles_only": int((np.diff(bo["free_energy"]) > 1e-9 * abs(bo["free_energy"][1:])).sum()),
    "coupled": int((np.diff(cp["free_energy"]) > 1e-9 * abs(cp["free_energy"][1:])).sum()),
    "gb_only": int((np.diff(go["free_energy"]) > 1e-9 * abs(go["free_energy"][1:])).sum())}
S["nonfinite"] = int(bo["nonfinite"].sum() + cp["nonfinite"].sum() + go["nonfinite"].sum())
json.dump(S, open(R("results/summary.json"), "w"), indent=1)
for k, v in S.items():
    print(f"{k:34s} {v}")

# ---------------------------------------------------------------- Fig 1: GB position vs time
zII = NZ - cp_c["vol_II"] / (NX * NY)          # lower edge of grain II (= GB, or top of the gas sheet)
S["coupled_t_gb_vanished"] = first(t_cp, cp_c["cols_real_gb"] == 0)
S["coupled_t_cov_1"] = first(t_cp, cp["coverage"] >= 0.9999)
S["coupled_zII_final"] = float(zII[-1])
mt = t_cp >= (S["coupled_t_gb_vanished"] or 1e9) + 10
if mt.sum() > 20:
    S["coupled_sheet_front_velocity_after_gb_gone"] = float(np.polyfit(t_cp[mt], zII[mt], 1)[0])
    S["coupled_sheet_centroid_velocity_after_gb_gone"] = float(np.polyfit(t_cp[mt], cp_c["bub_z_mean"][mt], 1)[0])
m1 = (t_cp >= 10) & (cp["coverage"] < 0.94)
S["coupled_gb_velocity_while_cov_below_0.94"] = float(np.polyfit(t_cp[m1], cp_c["gb_real_mean"][m1], 1)[0])
m2 = (t_go >= 10) & (t_go <= t_cp[m1][-1])
S["gb_only_velocity_same_window"] = float(np.polyfit(t_go[m2], go["gb_z_mean"][m2], 1)[0])
json.dump(S, open(R("results/summary.json"), "w"), indent=1)
for k in [k for k in S if k.startswith("coupled_t_gb_van") or "velocity" in k or k.endswith("cov_1") or "zII" in k]:
    print(f"{k:34s} {S[k]}")

fig, ax = plt.subplots(figsize=(8.6, 5.4))
ax.fill_between(t_cp, cp_c["bub_zbot_mean"], cp_c["bub_ztop_mean"], color="#bbbbbb", alpha=0.45, lw=0,
                label="GB + bubbles: bubble layer (mean bottom to mean top)")
ax.plot(t_go, go["gb_z_mean"], color="#1f77b4", lw=2.2, label="GB only (no bubbles), $\\Delta f$ = 0.1")
ax.plot(t_cp, cp_c["gb_real_mean"], color="#d62728", lw=2.2, label="GB + 33 bubbles, $\\Delta f$ = 0.1: GB height (where a GB exists)")
ax.plot(t_cp, zII, color="#d62728", lw=1.3, ls="--", label="GB + bubbles: lower edge of grain II (= top of gas sheet once GB is gone)")
ax.plot(t_cp, cp_c["bub_z_mean"], color="#555555", lw=1.3, ls=":", label="GB + bubbles: bubble centroid height")
ax.plot(t_bo, bo_c["gb_real_mean"], color="#2ca02c", lw=1.5, label="bubbles only ($\\Delta f$ = 0): GB height")
ax.axhline(Z_TOP, color="k", lw=0.8)
ax.text(395, Z_TOP - 2.4, "top wall (z = 63)", fontsize=9, ha="right")
for lev, ls in ((0.8, "-."), (0.9, "--")):
    tc = S[f"bubbles_only_t_cov_{lev}"]
    if tc is not None:
        ax.axvline(tc, color="#9467bd", ls=ls, lw=1)
        pass
ax.text(S["bubbles_only_t_cov_0.9"] + 1.5, 4, "bubbles only ($\\Delta f$ = 0):\n"
        f"coverage 0.8 at t = {S['bubbles_only_t_cov_0.8']:.0f} (dash-dot)\ncoverage 0.9 at t = {S['bubbles_only_t_cov_0.9']:.0f} (dashed)",
        color="#9467bd", fontsize=8)
tv = S["coupled_t_gb_vanished"]
if tv is not None:
    ax.axvline(tv, color="#d62728", ls=":", lw=1)
    ax.text(tv + 1.5, 14, f"GB + bubbles: coverage = 1,\nGB fully replaced by\ngas sheet at t = {tv:.0f}", color="#d62728", fontsize=8)
ax.set_xlabel("time (dimensionless)")
ax.set_ylabel("height z (cells)")
ax.set_ylim(0, NZ + 1)
ax.set_xlim(0, max(t_go[-1], t_cp[-1]))
ax.legend(loc="upper left", fontsize=7.6, framealpha=0.92)
ax.set_title("Grain-boundary position: free GB vs GB carrying 33 growing bubbles ($\\theta$ = 120°)", fontsize=10.5)
fig.tight_layout()
fig.savefig(R("results/fig1_gb_position.png"), dpi=170)

# ---------------------------------------------------------------- Fig 2: A1/Atot vs coverage
fig, ax = plt.subplots(figsize=(7.2, 5.2))
ax.plot(bo["coverage"], bo["largest_area_fraction"], color="#2ca02c", lw=1.8, label="bubbles only ($\\Delta f$ = 0, GB fixed)")
ax.plot(cp["coverage"], cp["largest_area_fraction"], color="#d62728", lw=1.8, label="GB + bubbles ($\\Delta f$ = 0.1, GB free to move)")
for d, t, col in ((bo, t_bo, "#2ca02c"), (cp, t_cp, "#d62728")):
    for tm in (5, 10, 20):
        i = np.searchsorted(t, tm)
        if i < len(t) and abs(t[i] - tm) < 0.06:
            ax.plot(d["coverage"][i], d["largest_area_fraction"][i], "o", ms=4, color=col)
            ax.annotate(f"t={tm}", (d["coverage"][i], d["largest_area_fraction"][i]), fontsize=7.5, color=col,
                        xytext=(4, -11 if col == "#2ca02c" else 5), textcoords="offset points")
ax.set_xlabel("projected coverage  $A_{tot}/A_{GB}$")
ax.set_ylabel("largest bubble area / total bubble area  $A_1/A_{tot}$")
ax.set_ylim(0, 1.04)
ax.legend(loc="upper left", fontsize=9)
ax.set_title("Coalescence path: $A_1/A_{tot}$ vs coverage (33 bubbles, $\\theta$ = 120°)", fontsize=10.5)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(R("results/fig2_largest_fraction_vs_coverage.png"), dpi=170)


# ---------------------------------------------------------------- Fig 3: xz slices (coupled)
def read_map(fn):
    b = open(fn, "rb").read()
    o = 0
    thick = np.frombuffer(b, np.uint8, NX * NY, o).reshape(NY, NX); o += NX * NY
    gbz = np.frombuffer(b, np.float32, NX * NY, o).reshape(NY, NX); o += 4 * NX * NY
    zb = np.frombuffer(b, np.uint8, NX * NY, o).reshape(NY, NX); o += NX * NY
    zt = np.frombuffer(b, np.uint8, NX * NY, o).reshape(NY, NX); o += NX * NY
    sc = np.frombuffer(b, np.float32, NZ * NX, o).reshape(NZ, NX); o += 4 * NZ * NX
    sn = np.frombuffer(b, np.float32, NZ * NX, o).reshape(NZ, NX)
    return thick, gbz, zb, zt, sc, sn


maps = sorted(glob.glob(R("runs/coupled_df0.1/maps/t*.bin")))
if maps:
    tm = np.array([float(os.path.basename(f)[1:-4]) for f in maps])
    want = [0, 10, 25, 50, 100, 150, 200, tm[-1]]
    pick = sorted({int(np.argmin(abs(tm - w))) for w in want})
    fig, axs = plt.subplots(len(pick), 2, figsize=(11, 1.45 * len(pick) + 0.6), squeeze=False)
    for r, i in enumerate(pick):
        thick, gbz, zb, zt, sc, sn = read_map(maps[i])
        a = axs[r, 0]
        a.imshow(sc, origin="lower", cmap="Greys", vmin=0, vmax=1, aspect="auto", extent=(0, NX, 0, NZ))
        a.contour(np.arange(NX) + 0.5, np.arange(NZ) + 0.5, np.ma.masked_where(sc > 0.5, sn), levels=[0], colors="r", linewidths=1.2)
        a.set_ylabel(f"t={tm[i]:.0f}\nz", fontsize=8)
        a.tick_params(labelsize=7)
        a = axs[r, 1]
        a.imshow(np.where(thick >= 3, 1.0, 0.0), cmap="Greys", vmin=0, vmax=1.6, origin="lower")
        a.imshow(np.ma.masked_invalid(gbz.astype(float)), cmap="viridis", vmin=15, vmax=NZ, origin="lower", alpha=0.75)
        a.set_xticks([]); a.set_yticks([])
    axs[0, 0].set_title("xz slice y = 64: black = gas, red = GB", fontsize=9)
    axs[0, 1].set_title("top view: grey = covered, purple = GB", fontsize=9)
    fig.tight_layout()
    fig.savefig(R("results/fig3_xz_slices.png"), dpi=150)
print("figures written")
