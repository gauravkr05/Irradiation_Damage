import numpy as np
import matplotlib
matplotlib.use("Agg")   # no GUI window -> avoids the GTK/Gdk segfault
import matplotlib.pyplot as plt

# ---- compare all three D(c) cases at ONE time ----
compare_time = 200      # change to 0, 500 or 1000 to compare at another time

cases = [
    ("P1_Q0",  1.0,  0.0),
    ("P1_Q1",  1.0,  1.0),
    ("P1_Q10", 1.0, 10.0),
]

plt.figure(figsize=(7, 5))
for tag, P, Q in cases:
    fname = "%s_time_%d/c_time_%d.dat" % (tag, compare_time, compare_time)
    data = np.loadtxt(fname, skiprows=1)
    plt.plot(data[:, 0], data[:, 1],
             label="D(c) = %g + %g$c^2$" % (P, Q))

plt.xlabel("x")
plt.ylabel("c")
plt.title("Comparison of D(c) at time = %d" % compare_time)
plt.ylim(-0.05, 1.05)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("compare_time_%d.png" % compare_time, dpi=150, bbox_inches="tight")
plt.close()

print("wrote compare_time_%d.png" % compare_time)
