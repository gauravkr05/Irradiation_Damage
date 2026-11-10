import numpy as np
import matplotlib
matplotlib.use("Agg")   # no GUI window -> avoids the GTK/Gdk segfault
import matplotlib.pyplot as plt

tag = "ac"
save_times = [0, 1, 5, 200, 400]
nx = 512
block_lo, block_hi = (nx - 50) // 2, (nx - 50) // 2 + 50


def load(tm):
    return np.loadtxt("%s_time_%d/eta_time_%d.dat" % (tag, tm, tm), skiprows=1)


# ---------- all times, full domain ----------
plt.figure(figsize=(9, 5))
for tm in save_times:
    d = load(tm)
    plt.plot(d[:, 0], d[:, 1], label="time = %d" % tm)

plt.xlabel("x")
plt.ylabel(r"$\eta$")
plt.title(r"1D Allen-Cahn,  $\kappa = 1$, $L = 1$   (I.C.: background 0.2, block of 50)")
plt.xlim(0, nx)
plt.ylim(-0.05, 1.05)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("%s_profiles.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

# ---------- zoom on the block, to see the interfaces ----------
plt.figure(figsize=(7, 5))
for tm in save_times:
    d = load(tm)
    plt.plot(d[:, 0], d[:, 1], label="time = %d" % tm)

plt.xlabel("x")
plt.ylabel(r"$\eta$")
plt.title("Zoom on the block (interfaces)")
plt.xlim(block_lo - 25, block_hi + 25)
plt.ylim(-0.05, 1.05)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("%s_zoom.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

# ---------- one figure per saved time ----------
for tm in save_times:
    folder = "%s_time_%d" % (tag, tm)
    d = load(tm)

    plt.figure(figsize=(9, 5))
    plt.plot(d[:, 0], d[:, 1])
    plt.xlabel("x")
    plt.ylabel(r"$\eta$")
    plt.title(r"1D Allen-Cahn,  time = %d" % tm)
    plt.xlim(0, nx)
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    plt.savefig("%s/eta_time_%d.png" % (folder, tm), dpi=150, bbox_inches="tight")
    plt.close()

print("plots written for %s" % tag)
