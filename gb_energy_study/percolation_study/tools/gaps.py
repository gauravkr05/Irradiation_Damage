"""Minimum edge-to-edge gap (periodic, in cells) between distinct covered components, per thickness map.
Covered = >= 3 cells with c >= 0.5 (same as the C++ metric). Output: results/gaps.csv"""
import glob, os, csv, sys
import numpy as np
from analyze import ROOT, OUT, label

STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 1   # use every STEP-th map
rows = []
for d in sorted(glob.glob(f"{ROOT}/runs/N*_bGB*")):
    case = os.path.basename(d)
    for f in sorted(glob.glob(f"{d}/thick/t*.npy"))[::STEP]:
        t = float(os.path.basename(f)[1:-4]); m = np.load(f) >= 3
        lab, _ = label(m); n = lab.max() + 1
        if n < 2:
            rows.append(dict(case=case, t=t, n=int(n), min_gap="")); continue
        # boundary cells of each component
        nb = np.zeros_like(m)
        for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nb |= ~np.roll(m, (dy, dx), (0, 1))
        bnd = m & nb
        ys, xs = np.nonzero(bnd); ls_ = lab[ys, xs]
        best = np.inf
        for k in range(n):
            a = ls_ == k; b = ~a
            if not b.any(): continue
            dy = np.abs(ys[a][:, None] - ys[b][None, :]); dy = np.minimum(dy, 128 - dy)
            dx = np.abs(xs[a][:, None] - xs[b][None, :]); dx = np.minimum(dx, 128 - dx)
            best = min(best, np.sqrt(dy ** 2 + dx ** 2).min())
        rows.append(dict(case=case, t=t, n=int(n), min_gap=round(float(best) - 1, 3)))   # -1: centre-to-centre -> empty cells between
    print(case, rows[-1])
with open(f"{OUT}/gaps.csv", "w") as fh:
    w = csv.DictWriter(fh, ["case", "t", "n", "min_gap"]); w.writeheader(); w.writerows(rows)
