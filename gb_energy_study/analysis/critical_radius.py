"""Classical estimate of the critical radius for the supersaturated matrix.
Grain matrix with n relaxed: n^2 = 1 - 2 bS c^2  ->  g_m(c) = c^4/4 - c^2/2 - (1-2bS c^2)^2/4
Bubble (n = 0):                                       g_b(c) = c^4/4 - c^2/2
Parallel-tangent driving force per unit bubble volume:
  dg = A0 [ g_m(cinf) + g_m'(cinf)(c_b - cinf) - g_b(c_b) ],  g_b'(c_b) = g_m'(cinf)
R* = 2 gamma_PG / dg   (sphere; for a GB lens the critical radius of CURVATURE is the same)"""
import numpy as np, json
from scipy.optimize import brentq
A0, bS = 4.0, 1.0
gm  = lambda c: c**4/4 - c**2/2 - (1-2*bS*c**2)**2/4
gmp = lambda c: c**3 - c + 2*bS*c*(1-2*bS*c**2)
gb  = lambda c: c**4/4 - c**2/2
gbp = lambda c: c**3 - c
gPG = json.load(open('interface_table.json'))['gPG_lattice']
print(f"gamma_PG (lattice) = {gPG:.4f}")
print(f"matrix spinodal (g_m''=0): c = {np.sqrt(1/9):.4f}  (with bS=1)")
print(" cinf   mu_inf   c_b     dg_v    R*_sphere(cells)  R*(nm)")
out={}
for cinf in (0.05,0.10,0.15,0.20,0.25,0.30):
    s = gmp(cinf)
    cb = brentq(lambda c: gbp(c)-s, 1.0, 1.5)
    dg = A0*(gm(cinf) + s*(cb-cinf) - gb(cb))
    R = 2*gPG/dg
    out[cinf]=R
    print(f" {cinf:.2f}  {A0*s:.4f}  {cb:.4f}  {dg:.4f}   {R:6.2f}           {10*R:6.1f}")
# matrix concentration the initial bubbles leave behind if they hold all the gas -> drives nothing:
