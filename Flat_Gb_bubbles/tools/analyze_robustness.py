"""analyze_robustness.py - is the coalescence path A1/Atot(coverage) independent of GB motion for other
contact angles and bubble numbers?

Each case = one bubbles-only run (dfBias = 0, GB fixed) and one coupled run (GB free to migrate), started
from the same initial condition; only the bias differs inside a case.

  ref      theta = 120 deg (bGB 2.0), N = 33   runs/bubbles_only          vs runs/coupled_df0.1          (df 0.1)
  theta128 theta = 128 deg (bGB 1.5), N = 33   runs/theta128_bubbles_only vs runs/theta128_coupled_df0.1 (df 0.1)
  N22      theta = 120 deg (bGB 2.0), N = 22   runs/N22_bubbles_only      vs runs/N22_coupled_df0.1      (df 0.1)
  theta104 theta = 104 deg (bGB 3.0), N = 33   runs/theta104_bubbles_only vs runs/theta104_coupled_df0.3 (df 0.3;
           at df 0.1 the free GB is lattice-pinned, runs/theta104_gb_only_df0.1)
  df05     theta = 120, N = 33, df = 0.5 (GB ~5x faster)          runs/bubbles_only     vs runs/coupled_df0.5
  L10      theta = 120, N = 33, LGB = 10, df = 0.1 (GB 10x faster, same thermodynamics)
                                                                 runs/L10_bubbles_only vs runs/L10_coupled_df0.1

Noise floor: runs/bubbles_only vs runs/L10_bubbles_only - both GB fixed (df = 0), they differ only in the
grain-field relaxation rate, so their path difference is the size of a "kinetics-only" change.

Writes results/fig4_percolation_robustness.png, fig5_all_paths.png, fig6_gb_motion_new_cases.png,
fig7_shift_vs_gb_motion.png and results/robustness_summary.json.
"""
import csv, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(ROOT, *p)
NX = NY = 128
LEVELS = (0.25, 0.5, 0.9)   # A1/Atot thresholds at which the coverage of the two runs is compared
SCAN = np.round(np.arange(0.10, 0.951, 0.05), 2)   # finer scan for the mean horizontal shift between two paths

CASES = [
    dict(key="ref", title="reference: $\\theta$ = 120°, N = 33", df=0.1, N=33,
         bo="bubbles_only", cp="coupled_df0.1", go="runs/gb_only_df0.1"),
    dict(key="theta128", title="contact angle changed: $\\theta$ = 128°, N = 33", df=0.1, N=33,
         bo="theta128_bubbles_only", cp="theta128_coupled_df0.1", go="runs/theta128_gb_only_df0.1"),
    dict(key="N22", title="bubble number changed: $\\theta$ = 120°, N = 22", df=0.1, N=22,
         bo="N22_bubbles_only", cp="N22_coupled_df0.1", go="runs/gb_only_df0.1"),
    dict(key="theta104", title="$\\theta$ = 104°, N = 33 ($\\Delta f$ = 0.3: GB pinned at 0.1)", df=0.3, N=33,
         bo="theta104_bubbles_only", cp="theta104_coupled_df0.3", go="runs/theta104_gb_only_df0.3"),
    dict(key="df05", title="strong test A: $\\theta$ = 120°, N = 33, $\\Delta f$ = 0.5", df=0.5, N=33,
         bo="bubbles_only", cp="coupled_df0.5", go="runs/gb_only_df0.5"),
    dict(key="L10", title="strong test B: $\\theta$ = 120°, N = 33, $L_{GB}$ = 10, $\\Delta f$ = 0.1", df=0.1, N=33,
         bo="L10_bubbles_only", cp="L10_coupled_df0.1", go="runs/L10_gb_only_df0.1", extra_bo="bubbles_only"),
]


def load(fn):
    if not os.path.exists(fn):
        return None
    raw = open(fn, "rb").read().replace(b"\x00", b"").decode().splitlines()
    rows = list(csv.DictReader(raw))
    if not rows:
        return None
    return {k: np.array([float(r[k]) if r[k] not in ("", "nan", None) else np.nan for r in rows]) for k in rows[0]}


def first(t, cond):
    i = np.flatnonzero(cond)
    return (float(t[i[0]]), int(i[0])) if len(i) else (None, None)


def cov_at_levels(m, levels):
    """coverage at which A1/Atot first reaches each level (None if never)"""
    out = []
    for lv in levels:
        _, i = first(m["time_dimensionless"], m["largest_area_fraction"] >= lv)
        out.append(np.nan if i is None else float(m["coverage"][i]))
    return np.array(out)


def path_shift(a, b):
    """horizontal distance coverage_b - coverage_a between two A1/Atot(coverage) paths, scanned over SCAN"""
    d = cov_at_levels(b, SCAN) - cov_at_levels(a, SCAN)
    d = d[np.isfinite(d)]
    # single-level values can jump by a whole coalescence step when a plateau sits right at a scan level,
    # so the median / interquartile range is the robust summary; min/max are kept for completeness
    return {"mean": float(d.mean()), "median": float(np.median(d)), "p25": float(np.percentile(d, 25)),
            "p75": float(np.percentile(d, 75)), "mean_abs": float(abs(d).mean()), "max_abs": float(abs(d).max()),
            "min": float(d.min()), "max": float(d.max()), "n_levels": int(len(d))}


def path_stats(m, N):
    t, nb, a1, cov = m["time_dimensionless"], m["number_of_projected_bubbles"], m["largest_area_fraction"], m["coverage"]
    s = {"t_end": float(t[-1]), "cov_max": float(cov.max()), "n_end": int(nb[-1])}
    tf, i = first(t, nb < N)
    s["t_first_merger"], s["cov_first_merger"] = tf, (float(cov[i]) if i is not None else None)
    ts, i = first(t, nb == 1)
    s["t_single"], s["cov_single"] = ts, (float(cov[i]) if i is not None else None)
    for lv in LEVELS:
        tl, i = first(t, a1 >= lv)
        s[f"t_A1_{lv}"], s[f"cov_A1_{lv}"] = tl, (float(cov[i]) if i is not None else None)
    s["max_abs_mass_drift"] = float(np.nanmax(abs(m["mass_drift"])))
    F = m["free_energy"]
    s["free_energy_increases"] = int((np.diff(F) > 1e-9 * abs(F[1:])).sum())
    s["nonfinite"] = int(np.nansum(m["nonfinite"]))
    s["max_count_above_N"] = int(max(0, np.nanmax(nb) - N))   # >0: extra (new / split) projected gas pockets appeared
    return s


S, D = {}, {}
for c in CASES:
    bo, cp = load(R("runs", c["bo"], "bubble_metrics.csv")), load(R("runs", c["cp"], "bubble_metrics.csv"))
    cpg, go = load(R("runs", c["cp"], "gb_bubble.csv")), load(R(c["go"], "gb_position.csv"))
    D[c["key"]] = (bo, cp, cpg, go)
    e = {"df": c["df"], "N": c["N"]}
    if bo is not None: e["bubbles_only"] = path_stats(bo, c["N"])
    if cp is not None: e["coupled"] = path_stats(cp, c["N"])
    if go is not None:
        tg, zg = go["time"], go["gb_z_mean"]
        m = (tg >= 5) & (zg < 55)
        e["free_gb_velocity"] = float(np.polyfit(tg[m], zg[m], 1)[0]) if m.sum() > 10 else None
    if bo is not None and cp is not None:
        e["path_shift"] = path_shift(bo, cp)
        e["delta_cov_at_A1"] = {str(lv): (None if e["coupled"][f"cov_A1_{lv}"] is None or e["bubbles_only"][f"cov_A1_{lv}"] is None
                                          else e["coupled"][f"cov_A1_{lv}"] - e["bubbles_only"][f"cov_A1_{lv}"]) for lv in LEVELS}
        # percolation window of the coupled run: first merger -> A1/Atot >= 0.9; how far did the GB move in it?
        t0, t1 = e["coupled"]["t_first_merger"], e["coupled"]["t_A1_0.9"]
        if t0 is not None and t1 is not None and cpg is not None:
            tc, zc = cpg["time"], cpg["gb_real_mean"]
            e["window"] = [t0, t1]
            e["gb_z_window"] = [float(np.interp(t0, tc, zc)), float(np.interp(t1, tc, zc))]
            e["gb_displacement_in_window"] = e["gb_z_window"][1] - e["gb_z_window"][0]
            if go is not None:
                e["free_gb_displacement_in_window"] = float(np.interp(t1, go["time"], go["gb_z_mean"]) - np.interp(t0, go["time"], go["gb_z_mean"]))
                e["drag_ratio_in_window"] = e["gb_displacement_in_window"] / e["free_gb_displacement_in_window"]
        if cpg is not None:
            e["coupled_gb_z_end"] = float(cpg["gb_real_mean"][~np.isnan(cpg["gb_real_mean"])][-1]) if np.isfinite(cpg["gb_real_mean"]).any() else None
            e["coupled_t_gb_vanished"] = first(cpg["time"], cpg["cols_real_gb"] == 0)[0]
            if "window" in e:
                e["bubble_centroid_z_window"] = [float(np.interp(x, cpg["time"], cpg["bub_z_mean"])) for x in e["window"]]
            e["cols_gb_above_bubble_max"] = int(np.nanmax(cpg["cols_gb_above_bubble"]))
            # detachment: covered columns in which the real GB lies above the bubble top (bubble left behind)
            fr = cpg["cols_gb_above_bubble"] / np.maximum(cpg["cols_covered"], 1)
            e["detached_fraction_max"] = float(np.nanmax(fr))
            e["detached_fraction_end"] = float(fr[-1])
            if "window" in e:
                w = (cpg["time"] >= e["window"][0]) & (cpg["time"] <= e["window"][1])
                e["detached_fraction_max_in_window"] = float(np.nanmax(fr[w])) if w.any() else None
    S[c["key"]] = e
_a, _b = load(R("runs/bubbles_only/bubble_metrics.csv")), load(R("runs/L10_bubbles_only/bubble_metrics.csv"))
S["noise_floor"] = {"what": "bubbles_only (LGB 1) vs L10_bubbles_only (LGB 10), both GB fixed",
                    "path_shift": path_shift(_a, _b),
                    "delta_cov_at_A1": {str(lv): float(x) for lv, x in zip(LEVELS, cov_at_levels(_b, LEVELS) - cov_at_levels(_a, LEVELS))}}
json.dump(S, open(R("results/robustness_summary.json"), "w"), indent=1)


def fmt(x, n=3):
    return "-" if x is None else (f"{x:.{n}f}" if isinstance(x, float) else str(x))


for k, e in S.items():
    if k == "noise_floor":
        print("== noise floor (", e["what"], "):", {a: fmt(b) for a, b in e["path_shift"].items()}); continue
    print(f"== {k}  df={e['df']}  N={e['N']}  free GB v={fmt(e.get('free_gb_velocity'))}")
    for r in ("bubbles_only", "coupled"):
        if r in e:
            s = e[r]
            print(f"   {r:13s} t_end {s['t_end']:6.1f}  first merger t={fmt(s['t_first_merger'],1)} cov={fmt(s['cov_first_merger'])}"
                  f"  A1>=.25 cov={fmt(s['cov_A1_0.25'])}  A1>=.5 cov={fmt(s['cov_A1_0.5'])}  A1>=.9 cov={fmt(s['cov_A1_0.9'])}"
                  f"  single t={fmt(s['t_single'],1)} cov={fmt(s['cov_single'])}  cov_max={s['cov_max']:.3f} n_end={s['n_end']}"
                  f"  drift={s['max_abs_mass_drift']:.1e} dF+={s['free_energy_increases']} nonfin={s['nonfinite']}")
    if "delta_cov_at_A1" in e:
        print("   delta cov (coupled - bo) at A1 levels:", {a: fmt(b) for a, b in e["delta_cov_at_A1"].items()},
              " window", e.get("window"), " GB moved", fmt(e.get("gb_displacement_in_window"), 2),
              " free GB would move", fmt(e.get("free_gb_displacement_in_window"), 2))
        print("   path shift (coupled - fixed):", {a: fmt(b) for a, b in e["path_shift"].items()},
              " detached frac max", fmt(e.get("detached_fraction_max")), " in window", fmt(e.get("detached_fraction_max_in_window")))

# ------------------------------------------------------------------ Fig 4: one panel per case
G, RD = "#2ca02c", "#d62728"
fig, axs = plt.subplots(2, 3, figsize=(16.5, 9.4), sharex=True, sharey=True)
for ax, c in zip(axs.flat, CASES):
    bo, cp, cpg, go = D[c["key"]]
    e = S[c["key"]]
    if c.get("extra_bo"):
        xb = load(R("runs", c["extra_bo"], "bubble_metrics.csv"))
        ax.plot(xb["coverage"], xb["largest_area_fraction"], color="#888888", lw=1.1, ls=":",
                label="bubbles only, $L_{GB}$ = 1 (reference)")
    if bo is not None:
        ax.plot(bo["coverage"], bo["largest_area_fraction"], color=G, lw=1.9, label="bubbles only ($\\Delta f$ = 0, GB fixed)")
    if cp is not None:
        ax.plot(cp["coverage"], cp["largest_area_fraction"], color=RD, lw=1.9, ls="--",
                label=f"GB + bubbles ($\\Delta f$ = {c['df']}, GB moves)")
    for d, col, dy in ((bo, G, -11), (cp, RD, 5)):
        if d is None: continue
        t = d["time_dimensionless"]
        for tm in (5, 10, 20):
            i = np.searchsorted(t, tm)
            if i < len(t) and abs(t[i] - tm) < 0.06:
                ax.plot(d["coverage"][i], d["largest_area_fraction"][i], "o", ms=3.5, color=col)
                ax.annotate(f"t={tm}", (d["coverage"][i], d["largest_area_fraction"][i]), fontsize=7, color=col,
                            xytext=(4, dy), textcoords="offset points")
    txt = []
    if "gb_displacement_in_window" in e:
        txt.append(f"GB moved {e['gb_displacement_in_window']:.1f} cells during\npercolation (t = {e['window'][0]:.1f}-{e['window'][1]:.1f})")
    if "delta_cov_at_A1" in e and e["delta_cov_at_A1"]["0.5"] is not None:
        txt.append(f"coverage at $A_1/A_{{tot}}$ = 0.5:\n{e['bubbles_only']['cov_A1_0.5']:.3f} fixed, {e['coupled']['cov_A1_0.5']:.3f} moving")
    if txt:
        ax.text(0.98, 0.04, "\n".join(txt), transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#bbbbbb"))
    ax.set_title(c["title"], fontsize=10)
    ax.grid(alpha=0.3); ax.set_ylim(0, 1.04); ax.set_xlim(0.1, 1.01)
    ax.legend(loc="upper left", fontsize=8)
for ax in axs[1]: ax.set_xlabel("projected coverage  $A_{tot}/A_{GB}$")
for ax in axs[:, 0]: ax.set_ylabel("largest bubble area / total  $A_1/A_{tot}$")
fig.suptitle("Is the coalescence path independent of GB motion? (each panel: same IC, GB fixed vs GB migrating)", fontsize=11)
fig.tight_layout()
fig.savefig(R("results/fig4_percolation_robustness.png"), dpi=170)

# ------------------------------------------------------------------ Fig 5: all paths on one axis
CC = {"ref": "#1f77b4", "theta128": "#9467bd", "N22": "#ff7f0e", "theta104": "#8c564b", "df05": "#e377c2", "L10": "#17becf"}
LB = {"ref": "θ = 120°, N = 33", "theta128": "θ = 128°, N = 33", "N22": "θ = 120°, N = 22", "theta104": "θ = 104°, N = 33",
      "df05": "θ = 120°, N = 33", "L10": "θ = 120°, N = 33, $L_{GB}$ = 10"}
LB5 = dict(LB, theta104="θ = 104°, N = 33, Δf = 0.3", df05="θ = 120°, N = 33, Δf = 0.5")
fig, ax = plt.subplots(figsize=(8, 5.6))
for c in CASES:
    bo, cp, _, _ = D[c["key"]]
    if c["key"] == "df05": bo = None          # same GB-fixed run as the reference
    if bo is not None:
        ax.plot(bo["coverage"], bo["largest_area_fraction"], color=CC[c["key"]], lw=2, label=f"{LB5[c['key']]}: GB fixed")
    if cp is not None:
        ax.plot(cp["coverage"], cp["largest_area_fraction"], color=CC[c["key"]], lw=1.4, ls="--", label=f"{LB5[c['key']]}: GB moving")
ax.set_xlabel("projected coverage  $A_{tot}/A_{GB}$"); ax.set_ylabel("$A_1/A_{tot}$")
ax.set_ylim(0, 1.04); ax.set_xlim(0.1, 1.01); ax.grid(alpha=0.3); ax.legend(fontsize=7.8, loc="upper left", ncol=1)
ax.set_title("All coalescence paths (solid: GB fixed, dashed: GB moving)", fontsize=10.5)
fig.tight_layout()
fig.savefig(R("results/fig5_all_paths.png"), dpi=170)

# ------------------------------------------------------------------ Fig 6: did the GB really move?
fig, ax = plt.subplots(figsize=(8, 5.2))
for c in CASES[1:]:
    _, _, cpg, go = D[c["key"]]
    if go is not None:
        ax.plot(go["time"], go["gb_z_mean"], color=CC[c["key"]], lw=1.2, ls=":")
    if cpg is not None:
        ax.plot(cpg["time"], cpg["gb_real_mean"], color=CC[c["key"]], lw=2,
                label=f"{LB[c['key']]}, $\\Delta f$ = {c['df']}: GB with bubbles")
ax.plot([], [], color="k", ls=":", label="free GB (no bubbles), same θ and Δf")
ax.set_xlabel("time"); ax.set_ylabel("GB height z (cells; mean over columns where a GB exists)")
ax.set_xlim(0, 60); ax.set_ylim(15, 64); ax.grid(alpha=0.3); ax.legend(fontsize=8)
ax.set_title("GB position in the new coupled runs", fontsize=10.5)
fig.tight_layout()
fig.savefig(R("results/fig6_gb_motion_new_cases.png"), dpi=170)

# ------------------------------------------------------------------ Fig 7: does the path shift grow with GB motion?
NAME = {"ref": "reference\nθ 120°, N 33", "theta128": "θ 128°", "N22": "N 22", "theta104": "θ 104°, Δf 0.3",
        "df05": "Δf 0.5", "L10": "$L_{GB}$ 10"}
OFF = {"ref": (7, 5), "theta128": (7, -12), "N22": (7, 4), "theta104": (7, 4), "df05": (7, 4), "L10": (-9, 5)}
fig, ax = plt.subplots(figsize=(7.6, 4.8))
nf = S["noise_floor"]["path_shift"]
ax.axhspan(-nf["max_abs"], nf["max_abs"], color="#cfd4da", alpha=0.6, lw=0)
ax.text(2.97, -nf["max_abs"], "largest kinetics-only shift (GB fixed, $L_{GB}$ 1 vs 10)  ", fontsize=7.5, color="#555555",
        va="top", ha="right")
ax.axhline(0, color="#888888", lw=0.8)
for c in CASES:
    e = S[c["key"]]
    if "gb_displacement_in_window" not in e: continue
    x, ps = e["gb_displacement_in_window"], e["path_shift"]
    hollow = c["key"] == "df05"     # bias also changes the gas thermodynamics, not a pure GB-motion test
    ax.plot([x, x], [ps["p25"], ps["p75"]], color="#1f4e79", lw=2, alpha=0.6, solid_capstyle="round")
    ax.plot(x, ps["median"], "o", ms=8, mfc="white" if hollow else "#1f4e79", mec="#1f4e79", mew=1.6)
    ax.annotate(NAME[c["key"]], (x, ps["median"]), xytext=OFF[c["key"]], textcoords="offset points", fontsize=8,
                ha="right" if OFF[c["key"]][0] < 0 else "left")
ax.set_xlabel("GB displacement during coalescence (first merger → $A_1/A_{tot}$ = 0.9), cells")
ax.set_ylabel("path shift: coverage(GB moving) − coverage(GB fixed)")
ax.set_xlim(0, 3); ax.grid(alpha=0.3)
ax.set_ylim(-0.045, 0.045)
ax.set_title("Path shift vs GB motion (dot: median over $A_1/A_{tot}$ = 0.10-0.95; bar: 25-75 %)\n"
             "hollow: Δf = 0.5 also changes the gas thermodynamics", fontsize=9.5)
fig.tight_layout()
fig.savefig(R("results/fig7_shift_vs_gb_motion.png"), dpi=170)
print("done")
