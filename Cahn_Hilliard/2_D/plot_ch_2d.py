import numpy as np
import matplotlib
matplotlib.use("Agg")   # no GUI window -> avoids the GTK/Gdk segfault
import matplotlib.pyplot as plt

tag = "ch2d"
save_times = [0, 50, 200, 800]
nx, ny = 100, 100


def load(tm):
    """read the x y c file and reshape into an (nx, ny) array"""
    d = np.loadtxt("%s_time_%d/c_time_%d.dat" % (tag, tm, tm), skiprows=1)
    return d[:, 2].reshape(nx, ny)


# ---------- one panel per save time, all in one figure ----------
fig, ax = plt.subplots(1, len(save_times), figsize=(4*len(save_times), 4.2))
for a, tm in zip(ax, save_times):
    field = load(tm)
    im = a.imshow(field.T, origin="lower", cmap="coolwarm",
                  vmin=0.0, vmax=1.0, extent=[0, nx, 0, ny])
    a.set_title("time = %d" % tm)
    a.set_xlabel("x")
a.figure.colorbar(im, ax=ax, fraction=0.02, pad=0.02, label="c")
ax[0].set_ylabel("y")
plt.suptitle("2D Cahn-Hilliard,  I.C. = 0.5 + noise")
plt.savefig("%s_panels.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

# ---------- one figure per save time ----------
for tm in save_times:
    field = load(tm)
    plt.figure(figsize=(5.5, 4.6))
    plt.imshow(field.T, origin="lower", cmap="coolwarm",
               vmin=0.0, vmax=1.0, extent=[0, nx, 0, ny])
    plt.colorbar(label="c")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("2D Cahn-Hilliard,  time = %d" % tm)
    plt.savefig("%s_time_%d/c_time_%d.png" % (tag, tm, tm),
                dpi=150, bbox_inches="tight")
    plt.close()

# ---------- t = 0 on its own scale (noise is only +/- 0.005) ----------
field = load(0)
plt.figure(figsize=(5.5, 4.6))
plt.imshow(field.T, origin="lower", cmap="coolwarm", extent=[0, nx, 0, ny])
plt.colorbar(label="c")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Initial condition:  c = 0.5 + noise   (zoomed)")
plt.savefig("%s_IC_zoom.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

# ---------- middle-row slice, to compare against the 1D run ----------
plt.figure(figsize=(7, 5))
for tm in save_times:
    field = load(tm)
    plt.plot(np.arange(nx), field[:, ny//2], label="time = %d" % tm)
plt.xlabel("x")
plt.ylabel("c")
plt.title("2D Cahn-Hilliard,  middle-row slice (y = %d)" % (ny//2))
plt.ylim(-0.05, 1.05)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("%s_slice.png" % tag, dpi=150, bbox_inches="tight")
plt.close()

print("plots written for %s" % tag)
