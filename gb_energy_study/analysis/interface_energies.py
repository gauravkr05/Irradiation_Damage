"""
interface_energies.py

Equilibrium 1-D interface profiles and excess free energies for the
free-energy functional used in phasefield_gb_bubbles3d.cpp:

  F = int [ A0 f(c,n1,n2) + j/2 |grad c|^2 + jGB/2 (|grad n1|^2+|grad n2|^2) ]

  f = c^4/4 - c^2/2 + n1^4/4 - n1^2/2 + n2^4/4 - n2^2/2
    + bGB/2 n1^2 n2^2 + bS c^2 (n1^2+n2^2)

Bulk states (all have f = -1/4):
  grain 1 : (c,n1,n2) = (0,1,0)    grain 2 : (0,0,1)    bubble : (1,0,0)
so the flat-interface coexistence is at mu = 0 and each excess energy is
  sigma = int [ A0 (f + 1/4) + gradient terms ] dz.

Two discretisations are evaluated:
  * lattice  : dx = 1 cell, nearest-neighbour differences, i.e. exactly
               the energy whose variational derivative the C++ code uses
  * continuum: dx = 0.05 cell (grid-converged reference)

Profiles are obtained by gradient flow with fixed end states, run until the
energy change per step is at round-off. Units: energy per cell^2 in the
dimensionless units of the code (A0, j, jGB as given).
"""
import numpy as np


def f_bulk(c, a, b, bS, bGB):
    return (c**4 / 4 - c**2 / 2 + a**4 / 4 - a**2 / 2 + b**4 / 4 - b**2 / 2
            + bGB / 2 * a**2 * b**2 + bS * c**2 * (a**2 + b**2))


def dfdc(c, a, b, bS):
    return c**3 - c + 2 * bS * c * (a**2 + b**2)


def dfda(c, a, b, bS, bGB):
    return a**3 - a + bGB * a * b**2 + 2 * bS * c**2 * a


def excess(c, a, b, P, dx):
    A0, bS, bGB, j, jGB = P["A0"], P["bS"], P["bGB"], P["j"], P["jGB"]
    bulk = np.sum(A0 * (f_bulk(c, a, b, bS, bGB) + 0.25)) * dx
    grad = (j / 2 * np.sum(np.diff(c)**2) + jGB / 2 * (np.sum(np.diff(a)**2) + np.sum(np.diff(b)**2))) / dx
    return bulk + grad


def lap(u, dx):
    L = np.zeros_like(u)
    L[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / dx**2
    return L


def relax(c, a, b, P, dx, which, tol=1e-13, maxit=4_000_000):
    """gradient flow; `which` = set of fields that evolve; ends are fixed."""
    A0, bS, bGB, j, jGB = P["A0"], P["bS"], P["bGB"], P["j"], P["jGB"]
    kmax = 4 / dx**2
    dt = 0.4 / (A0 * 6 + max(j, jGB) * kmax)
    E_old = excess(c, a, b, P, dx)
    for it in range(maxit):
        if "c" in which:
            gc = A0 * dfdc(c, a, b, bS) - j * lap(c, dx)
        if "a" in which:
            ga = A0 * dfda(c, a, b, bS, bGB) - jGB * lap(a, dx)
        if "b" in which:
            gb = A0 * dfda(c, b, a, bS, bGB) - jGB * lap(b, dx)
        if "c" in which:
            c[1:-1] -= dt * gc[1:-1]
        if "a" in which:
            a[1:-1] -= dt * ga[1:-1]
        if "b" in which:
            b[1:-1] -= dt * gb[1:-1]
        if it % 2000 == 0:
            E = excess(c, a, b, P, dx)
            if abs(E - E_old) < tol:
                return E, it
            E_old = E
    return excess(c, a, b, P, dx), maxit


def gb_energy(P, dx, L=40.0, z0=0.0):
    z = np.arange(-L / 2, L / 2 + dx / 2, dx) - z0
    a = 0.5 * (1 - np.tanh(z / 1.5))
    b = 0.5 * (1 + np.tanh(z / 1.5))
    c = np.zeros_like(z)
    E, it = relax(c, a, b, P, dx, {"a", "b"})
    # width: distance over which n1 goes 0.9 -> 0.1
    zz = z + z0
    w = np.interp(0.1, a[::-1], zz[::-1]) - np.interp(0.9, a[::-1], zz[::-1])
    return E, w, (zz, c, a, b)


def pg_energy(P, dx, L=40.0, z0=0.0):
    """bubble (z<0) / grain 1 (z>0); n2 = 0 everywhere."""
    z = np.arange(-L / 2, L / 2 + dx / 2, dx) - z0
    c = 0.5 * (1 - np.tanh(z / 1.5))
    a = 0.5 * (1 + np.tanh(z / 1.5))
    b = np.zeros_like(z)
    E, it = relax(c, a, b, P, dx, {"c", "a"})
    zz = z + z0
    w = np.interp(0.1, c[::-1], zz[::-1]) - np.interp(0.9, c[::-1], zz[::-1])
    return E, w, (zz, c, a, b)


BASE = dict(A0=4.0, bS=1.0, j=2.0, jGB=2.0)


def theta_young(gGB, gPG):
    """dihedral angle measured through the bubble: gGB = 2 gPG cos(theta/2)."""
    r = gGB / (2 * gPG)
    if r >= 1:
        return 0.0
    return np.degrees(2 * np.arccos(r))


if __name__ == "__main__":
    import sys
    print("=== bubble/grain interface, baseline bS=1, j=2, jGB=2 (independent of bGB) ===")
    for dx in (1.0, 0.05):
        for z0 in (0.0, 0.25, 0.5):
            P = dict(BASE, bGB=1.0)
            E, w, _ = pg_energy(P, dx, z0=z0 if dx == 1.0 else 0.0)
            print(f"dx={dx:5.2f} offset={z0:4.2f}  gamma_PG = {E:.5f}   10-90 width = {w:.3f} cells")
            if dx != 1.0:
                break
