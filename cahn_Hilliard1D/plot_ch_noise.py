import numpy as np
import matplotlib
matplotlib.use("Agg")   # no GUI window -> avoids the GTK/Gdk segfault
import matplotlib.pyplot as plt

tag = "ch_noise"
save_times = [0, 10, 50, 100, 200, 400, 800]
label = "1D Cahn-Hilliard,  nx = 512,  I.C. = 0.5 + noise"


def load(tm):
    return np.loadtxt("%s_time_%d/c_time_%d.dat" % (tag, tm, tm), skiprows=1)


# ---------- stacked panels: one row per time (clearest at nx = 512) ----------
fig, ax = plt.subplots(len(save_times), 1, figsize=(11, 2.0*len(save_times)),
                       sharex=True)
for a, tm in zip(ax, save_times):
    d = load(tm)
    a.plot(d[:, 0], d[:, 1], lw=1.4)
    a.set_ylim(-0.05, 1.05)
    a.set_ylabel("c")
    a.grid(True, alpha=0.3)
    a.text(0.995, 0.86, "t = %d" % tm, transform=a.transAxes,
           ha="right", va="top", fontsize=10,
           bbox=dict(fc="white", ec="0.7", alpha=0.9))
ax[-1].set_xlabel("x")
plt.suptitle(label, y=0.995)
plt.tight_layout()
plt.savefig("%s_panels.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

# ---------- all times overlaid ----------
plt.figure(figsize=(11, 5.5))
for tm in save_times:
    d = load(tm)
    plt.plot(d[:, 0], d[:, 1], lw=1.2, label="t = %d" % tm)
plt.xlabel("x")
plt.ylabel("c")
plt.title(label)
plt.ylim(-0.05, 1.05)
plt.legend(ncol=7, fontsize=8, loc="upper center")
plt.grid(True, alpha=0.3)
plt.savefig("%s_profiles.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

# ---------- one figure per saved time ----------
for tm in save_times:
    folder = "%s_time_%d" % (tag, tm)
    d = load(tm)
    plt.figure(figsize=(11, 4))
    plt.plot(d[:, 0], d[:, 1], lw=1.4)
    plt.xlabel("x")
    plt.ylabel("c")
    plt.title("%s,  time = %d" % (label, tm))
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    plt.savefig("%s/c_time_%d.png" % (folder, tm), dpi=150, bbox_inches="tight")
    plt.close()

# ---------- t = 0 on its own scale (noise is only +/- 0.005) ----------
d = load(0)
plt.figure(figsize=(11, 4))
plt.plot(d[:, 0], d[:, 1], lw=1.0)
plt.xlabel("x")
plt.ylabel("c")
plt.title("Initial condition:  c = 0.5 + noise   (zoomed)")
plt.grid(True, alpha=0.3)
plt.savefig("%s_IC_zoom.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

# ---------- number of domains vs time (coarsening) ----------
print(" time   domains (c crosses 0.5)   mass")
for tm in save_times:
    d = load(tm)[:, 1]
    crossings = np.sum(np.diff(np.sign(d - 0.5)) != 0)
    print("%6d   %-22d %.4f" % (tm, crossings // 2, d.sum()))

print("plots written for %s" % tag)
