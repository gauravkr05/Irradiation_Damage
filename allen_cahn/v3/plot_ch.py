import os

import numpy as np
import matplotlib
matplotlib.use("Agg")   # no GUI window -> avoids the GTK/Gdk segfault
import matplotlib.pyplot as plt

tag = "ch"
field = "c"
save_times = [0, 1, 5, 200, 400]
nx = 512
dx = 1.0
block_lo, block_hi = (nx - 50) // 2, (nx - 50) // 2 + 50


def load(tm):
    return np.loadtxt("%s_time_%d/%s_time_%d.dat" % (tag, tm, field, tm), skiprows=1)


# ---------- mass check : the whole point of Cahn-Hilliard ----------
print("integral of %s   (dx = %g)" % (field, dx))
for tm in save_times:
    d = load(tm)
    total = d[:, 1].sum() * dx
    print("  time = %-5d  sum = %8.3f   mean = %6.4f" % (tm, total, total / (nx * dx)))
print("")

# ---------- all times, full domain ----------
plt.figure(figsize=(9, 5))
for tm in save_times:
    d = load(tm)
    plt.plot(d[:, 0], d[:, 1], label="time = %d" % tm)

plt.axhline(0.0, color="0.7", ls="--", lw=1)
plt.axhline(1.0, color="0.7", ls="--", lw=1)
plt.xlabel("x")
plt.ylabel(r"$c$")
plt.title(r"1D Cahn-Hilliard,  $\kappa = 1$, $M = 1$   (I.C.: background 0.2, block of 50)")
plt.xlim(0, nx)
plt.ylim(-0.15, 1.15)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("%s_profiles.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

# ---------- zoom on the block, to see the interfaces + depletion zone ----------
plt.figure(figsize=(7, 5))
for tm in save_times:
    d = load(tm)
    plt.plot(d[:, 0], d[:, 1], label="time = %d" % tm)

plt.axhline(0.2, color="0.7", ls=":", lw=1)
plt.xlabel("x")
plt.ylabel(r"$c$")
plt.title("Zoom on the block (interfaces and depletion zone)")
plt.xlim(block_lo - 60, block_hi + 60)
plt.ylim(-0.15, 1.15)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("%s_zoom.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

# ---------- one figure per saved time ----------
for tm in save_times:
    folder = "%s_time_%d" % (tag, tm)
    d = load(tm)
    total = d[:, 1].sum() * dx

    plt.figure(figsize=(9, 5))
    plt.plot(d[:, 0], d[:, 1])
    plt.xlabel("x")
    plt.ylabel(r"$c$")
    plt.title(r"1D Cahn-Hilliard,  time = %d   ($\int c\,dx$ = %.3f)" % (tm, total))
    plt.xlim(0, nx)
    plt.ylim(-0.15, 1.15)
    plt.grid(True, alpha=0.3)
    plt.savefig("%s/%s_time_%d.png" % (folder, field, tm), dpi=150, bbox_inches="tight")
    plt.close()

# ---------- Allen-Cahn vs Cahn-Hilliard, if the ac run is also there ----------
ac_file = "ac_time_400/eta_time_400.dat"
if os.path.exists(ac_file):
    ac = np.loadtxt(ac_file, skiprows=1)
    ch = load(400)

    plt.figure(figsize=(9, 5))
    plt.plot(ac[:, 0], ac[:, 1], color="tab:red",
             label=r"Allen-Cahn,  $\int\eta\,dx$ = %.1f" % (ac[:, 1].sum() * dx))
    plt.plot(ch[:, 0], ch[:, 1], color="tab:blue",
             label=r"Cahn-Hilliard,  $\int c\,dx$ = %.1f" % (ch[:, 1].sum() * dx))
    plt.xlabel("x")
    plt.ylabel("order parameter")
    plt.title("time = 400 : same initial condition, two dynamics")
    plt.xlim(0, nx)
    plt.ylim(-0.15, 1.15)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig("ac_vs_ch_t400.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("wrote ac_vs_ch_t400.png")

print("plots written for %s" % tag)
