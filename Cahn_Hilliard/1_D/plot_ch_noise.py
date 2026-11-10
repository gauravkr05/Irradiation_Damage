import numpy as np
import matplotlib
matplotlib.use("Agg")   # no GUI window -> avoids the GTK/Gdk segfault
import matplotlib.pyplot as plt

# ---- Cahn-Hilliard, I.C. = 0.5 + noise ----
tag = "ch_noise"
label = "Cahn-Hilliard,  I.C. = 0.5 + noise"
save_times = [0, 50, 200, 800]

# ---------- all times on one figure ----------
plt.figure(figsize=(7, 5))
for tm in save_times:
    data = np.loadtxt("%s_time_%d/c_time_%d.dat" % (tag, tm, tm), skiprows=1)
    plt.plot(data[:, 0], data[:, 1], label="time = %d" % tm)

plt.xlabel("x")
plt.ylabel("c")
plt.title(label)
plt.ylim(-0.05, 1.05)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("%s_profiles.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

# ---------- one figure per saved time ----------
for tm in save_times:
    folder = "%s_time_%d" % (tag, tm)
    data = np.loadtxt("%s/c_time_%d.dat" % (folder, tm), skiprows=1)

    plt.figure(figsize=(7, 5))
    plt.plot(data[:, 0], data[:, 1])
    plt.xlabel("x")
    plt.ylabel("c")
    plt.title("%s,  time = %d" % (label, tm))
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    plt.savefig("%s/c_time_%d.png" % (folder, tm), dpi=150, bbox_inches="tight")
    plt.close()

# ---------- time 0 on its own axis (noise is only +/- 0.005) ----------
data = np.loadtxt("%s_time_0/c_time_0.dat" % tag, skiprows=1)
plt.figure(figsize=(7, 5))
plt.plot(data[:, 0], data[:, 1])
plt.xlabel("x")
plt.ylabel("c")
plt.title("Initial condition:  c = 0.5 + noise   (zoomed)")
plt.grid(True, alpha=0.3)
plt.savefig("%s_IC_zoom.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

print("plots written for %s" % tag)
