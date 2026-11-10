"""Measure the dihedral angle (through the bubble) of a 2-D lens on the GB.
Method A: theta = 4 atan(h/a)        h = half-thickness on the axis, a = half-width at GB
Method B: circle fit to each cap's c=0.5 contour, excluding |z-ZGB| < zcut
          theta = 2 acos(-d/rho), d = fitted centre height relative to GB plane.
Both assume circular caps (true at equilibrium in 2-D)."""
import numpy as np, sys, glob, csv, os
NX, NZ = 160, 96; ZGB = NZ//2 - 0.5
def crossings_1d(v, level=0.5):
    out=[]
    for i in range(len(v)-1):
        if (v[i]-level)*(v[i+1]-level) < 0: out.append(i + (v[i]-level)/(v[i]-v[i+1]))
    return out
def contour_points(c, level=0.5):
    pts=[]
    for z in range(NZ):                      # horizontal crossings
        for x in crossings_1d(c[z,:]): pts.append((x, z))
    for x in range(NX):                      # vertical crossings
        for z in crossings_1d(c[:,x]): pts.append((x, z))
    return np.array(pts)
def fit_circle(p):
    A = np.c_[2*p[:,0], 2*p[:,1], np.ones(len(p))]; b = (p**2).sum(1)
    (xc, zc, k), *_ = np.linalg.lstsq(A, b, rcond=None)
    return xc, zc, np.sqrt(k + xc**2 + zc**2)
def measure(c, zcut=4.0):
    xc = NX/2
    col = c[:, int(round(xc))]; zz = crossings_1d(col)
    h = 0.5*((max(zz) - ZGB) + (ZGB - min(zz)))
    r1 = crossings_1d(c[int(ZGB-0.5),:]); r2 = crossings_1d(c[int(ZGB+0.5),:])
    a = 0.25*((max(r1)-min(r1)) + (max(r2)-min(r2)))
    thA = np.degrees(4*np.arctan(h/a))
    p = contour_points(c); thB=[]
    for sgn in (+1,-1):
        q = p[sgn*(p[:,1]-ZGB) > zcut]
        x0, z0, rho = fit_circle(q)
        d = sgn*(z0 - ZGB)
        thB.append(np.degrees(2*np.arccos(np.clip(-d/rho,-1,1))))
    return thA, np.mean(thB), h, a, thB
if __name__ == "__main__":
    for d in sys.argv[1:]:
        rows = list(csv.DictReader(open(os.path.join(d,'snapshots.csv'))))
        for r in rows:
            c = np.fromfile(os.path.join(d,r['c_file']),dtype=np.float32).reshape(NZ,NX)
            thA, thB, h, a, tb = measure(c)
            print(f"{d}  t={float(r['time_dimensionless']):7.1f}  h={h:6.2f} a={a:6.2f}  theta_A={thA:6.1f}  theta_B={thB:6.1f} (top {tb[0]:.1f}, bottom {tb[1]:.1f})")
