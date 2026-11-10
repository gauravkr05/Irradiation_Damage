"""Analysis of the percolation study. Reads only simulation output:
runs/N*_bGB*/bubble_metrics.csv (every 0.1) and runs/*/thick/t*.npy (column-thickness maps every 1.0).

Definitions (identical to the original study; nothing tuned):
  covered column : >= 3 cells with c >= 0.5 ; components 8-connected, periodic
  first merger   : first 0.1 output with fewer components than at t = 0   ("percolation begin")
  A1/Atot > 0.5  : first output where largest component > half the covered area
  single bubble  : first output with exactly 1 component                ("percolation end")
  spanning       : first thickness map (1.0 spacing) in which one component wraps the periodic box
Event check: between consecutive thickness maps each component at t is mapped to the components at
t+1 it overlaps.  A component with no successor = dissolved (Ostwald), a successor fed by >1 = merger.
"""
import glob, os, csv, json, sys
import numpy as np
from collections import deque, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f"{ROOT}/results"; os.makedirs(OUT, exist_ok=True)
THETA = {"1.2": 148.0, "1.5": 128.4, "2.0": 120.1, "3.0": 103.9}   # measured in the original study (+-6 deg)
NB = [(dy, dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dy or dx]


def load_metrics(d):
    rows = list(csv.DictReader(open(f"{d}/bubble_metrics.csv")))
    k = ["time_dimensionless", "number_of_projected_bubbles", "coverage", "largest_area_fraction",
         "porosity", "mass_drift", "free_energy"]
    return {c: np.array([float(r[c]) for r in rows]) for c in k}


def label(mask):
    """8-connected periodic labelling; returns labels and per-component 'wraps' flag."""
    H, W = mask.shape
    lab = -np.ones(mask.shape, int); off = {}
    wraps = []
    for s in zip(*np.nonzero(mask)):
        if lab[s] >= 0: continue
        cid = len(wraps); w = False
        lab[s] = cid; off[s] = (0, 0); q = deque([s])
        while q:
            y, x = q.popleft(); oy, ox = off[(y, x)]
            for dy, dx in NB:
                yy, xx = y + dy, x + dx
                oyy, oxx = oy + (yy // H), ox + (xx // W)   # unwrapped image index
                yy %= H; xx %= W
                if not mask[yy, xx]: continue
                if lab[yy, xx] < 0:
                    lab[yy, xx] = cid; off[(yy, xx)] = (oyy, oxx); q.append((yy, xx))
                elif off[(yy, xx)] != (oyy, oxx):
                    w = True   # same cell reached in two different periodic images -> component wraps
        wraps.append(w)
    return lab, wraps


def thick_series(d):
    fs = sorted(glob.glob(f"{d}/thick/t*.npy"))
    return [(float(os.path.basename(f)[1:-4]), np.load(f)) for f in fs]


def event_check(series):
    """Classify every decrease in component count between consecutive maps."""
    merges = dissol = 0; first_dissolution = None; log = []
    prev = None
    for t, th in series:
        lab, wr = label(th >= 3)
        if prev is not None:
            pt, plab, pn = prev
            n = lab.max() + 1
            succ = defaultdict(set); pred = defaultdict(set)
            both = (plab >= 0) & (lab >= 0)
            for a, b in set(zip(plab[both].tolist(), lab[both].tolist())):
                succ[a].add(b); pred[b].add(a)
            gone = [a for a in range(pn) if not succ[a]]
            mrg = sum(len(v) - 1 for v in pred.values() if len(v) > 1)
            born = [b for b in range(n) if not pred[b]]
            merges += mrg; dissol += len(gone)
            if gone and first_dissolution is None: first_dissolution = t
            if gone or mrg or born:
                log.append(dict(t0=pt, t1=t, n0=int(pn), n1=int(n), merged=mrg, dissolved=len(gone), new=len(born),
                                dissolved_sizes=[int((plab == a).sum()) for a in gone]))
        prev = (t, lab, lab.max() + 1)
    return merges, dissol, first_dissolution, log


def first(cond, t):
    i = np.nonzero(cond)[0]
    return int(i[0]) if len(i) else None


def analyse(d):
    m = load_metrics(d); t = m["time_dimensionless"]; n = m["number_of_projected_bubbles"]; cov = m["coverage"]
    N0 = int(n[0]); res = dict(case=os.path.basename(d), N=N0, tmax=float(t[-1]),
                               n_end=int(n[-1]), cov_end=float(cov[-1]), por_end=float(m["porosity"][-1]),
                               max_drift=float(np.abs(m["mass_drift"]).max()),
                               F_incr=int((np.diff(m["free_energy"]) > 1e-9 * np.abs(m["free_energy"][1:])).sum()),
                               n_increase=int((np.diff(n) > 0).sum()))
    for key, cond in (("first_merger", n < N0), ("half", m["largest_area_fraction"] > 0.5), ("single", n == 1)):
        i = first(cond, t)
        res[f"t_{key}"] = float(t[i]) if i is not None else None
        res[f"cov_{key}"] = float(cov[i]) if i is not None else None
        res[f"covprev_{key}"] = float(cov[i - 1]) if i else None
    # does it stay single?
    i = first(n == 1, t)
    res["single_stays"] = bool((n[i:] == 1).all()) if i is not None else None
    series = thick_series(d)
    res["t_span"] = res["cov_span"] = None
    for tt, th in series:
        _, wr = label(th >= 3)
        if any(wr):
            res["t_span"] = tt; res["cov_span"] = float(cov[np.argmin(np.abs(t - tt))]); break
    mg, ds, fd, log = event_check(series)
    res.update(merge_events=mg, dissolution_events=ds, first_dissolution_t=fd)
    json.dump(log, open(f"{OUT}/events_{res['case']}.json", "w"), indent=1)
    return res, m


if __name__ == "__main__":
    allres, allm = [], {}
    for d in sorted(glob.glob(f"{ROOT}/runs/N*_bGB*")):
        r, m = analyse(d); r["theta"] = THETA[r["case"].split("bGB")[1]]; allres.append(r); allm[r["case"]] = m
    keys = list(allres[0].keys())
    with open(f"{OUT}/events_summary.csv", "w") as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(allres)
    for r in allres:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
