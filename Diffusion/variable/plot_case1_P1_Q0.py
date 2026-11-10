import numpy as np
import matplotlib
matplotlib.use("Agg")   # no GUI window -> avoids the GTK/Gdk segfault
import matplotlib.pyplot as plt

# ---- CASE 1 : D(c) = 1 + 0*c^2 ----
tag = "P1_Q0"
P, Q = 1.0, 0.0
save_times = [0, 200, 500, 1000]

# ---------- all times on one figure ----------
plt.figure(figsize=(7, 5))
for tm in save_times:
    data = np.loadtxt("%s_time_%d/c_time_%d.dat" % (tag, tm, tm), skiprows=1)
    plt.plot(data[:, 0], data[:, 1], label="time = %d" % tm)

plt.xlabel("x")
plt.ylabel("c")
plt.title("1D Diffusion,  D(c) = %g + %g$c^2$" % (P, Q))
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
    plt.title("D(c) = %g + %g$c^2$,  time = %d" % (P, Q, tm))
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    plt.savefig("%s/c_time_%d.png" % (folder, tm), dpi=150, bbox_inches="tight")
    plt.close()

print("plots written for %s" % tag)
