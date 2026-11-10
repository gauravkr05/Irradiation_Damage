"""Detection-rule robustness of 'percolation begin' (first merger) and 'percolation end' (single
component) using the saved thickness maps (spacing 1.0). Rule k: covered column = >= k cells with
c >= 0.5 (k = 3 is the original rule). Reports, per run and rule, the first map time with fewer
components than at t = 0 and the first map time with exactly one component, plus the coverage at
that map under the rule itself and under the original (>= 3) rule. Output: results/rules.csv"""
import glob, os, csv
import numpy as np
from analyze import ROOT, OUT, label

rows = []
for d in sorted(glob.glob(f"{ROOT}/runs/N*_bGB*")):
    case = os.path.basename(d)
    maps = [(float(os.path.basename(f)[1:-4]), np.load(f)) for f in sorted(glob.glob(f"{d}/thick/t*.npy"))]
    for k in (1, 3, 5):
        n0 = None; ev = {}
        for t, th in maps:
            lab, _ = label(th >= k); n = lab.max() + 1
            if n0 is None: n0 = n
            if "begin" not in ev and n < n0: ev["begin"] = (t, (th >= k).mean(), (th >= 3).mean())
            if "end" not in ev and n == 1: ev["end"] = (t, (th >= k).mean(), (th >= 3).mean()); break
        r = dict(case=case, rule=f">={k} cells", n0=int(n0), last_map_t=maps[-1][0])
        for e in ("begin", "end"):
            t, cs, c3 = ev.get(e, (None, None, None))
            r[f"t_{e}"] = t; r[f"cov_own_{e}"] = None if cs is None else round(float(cs), 4)
            r[f"cov_orig_{e}"] = None if c3 is None else round(float(c3), 4)
        rows.append(r); print(r)
with open(f"{OUT}/rules.csv", "w") as fh:
    w = csv.DictWriter(fh, list(rows[0].keys())); w.writeheader(); w.writerows(rows)
