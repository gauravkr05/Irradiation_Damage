import numpy as np, glob, re, os
import matplotlib
matplotlib.use("Agg")   # no GUI window -> avoids the GTK/Gdk segfault
import matplotlib.pyplot as plt

folder = "output_ac_fd_2d"

# read grid size from the inputs file the code writes
p = {}
for line in open(os.path.join(folder, "inputs_taken.dat")):
    k, v = line.split("=")
    p[k.strip()] = float(v)
NX, NY = int(p["NX"]), int(p["NY"])

files = sorted(glob.glob(os.path.join(folder, "time*.dat")),
               key=lambda f: float(re.findall(r"time([\d.]+)\.dat", f)[0]))
times = [float(re.findall(r"time([\d.]+)\.dat", f)[0]) for f in files]

def load(f):
    return np.loadtxt(f)[:, 2].reshape(NX, NY)

# ---------- panel of snapshots ----------
pick = files[:6]
fig, ax = plt.subplots(1, len(pick), figsize=(3.4*len(pick), 3.8))
for a, f in zip(ax, pick):
    im = a.imshow(load(f).T, origin="lower", cmap="coolwarm", vmin=0, vmax=1,
                  extent=[0, NX, 0, NY])
    a.set_title("t = %g" % float(re.findall(r"time([\d.]+)\.dat", f)[0]))
    a.set_xlabel("x")
ax[0].set_ylabel("y")
fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02, label=r"$\eta$")
plt.suptitle("2D Allen-Cahn (explicit FD): circular grain shrinking")
plt.savefig("ac2d_panels.png", dpi=150, bbox_inches="tight")
plt.close()

# ---------- slice through the centre : this is the prof's 1D picture ----------
plt.figure(figsize=(9, 5))
for f in pick:
    t = float(re.findall(r"time([\d.]+)\.dat", f)[0])
    plt.plot(np.arange(NX), load(f)[:, NY//2], label="t = %g" % t)
plt.axhline(0.2, color="grey", ls=":", lw=1)
plt.xlabel("x"); plt.ylabel(r"$\eta$")
plt.title(r"Slice through the centre (y = NY/2) -- diameter 50 at t=0")
plt.ylim(-0.05, 1.05)
plt.legend(); plt.grid(True, alpha=0.3)
plt.savefig("ac2d_slice.png", dpi=150, bbox_inches="tight")
plt.close()

# ---------- radius vs time, against the curvature law ----------
R = []
for f in files:
    a = (load(f) > 0.5).sum()
    R.append(np.sqrt(a/np.pi) if a > 0 else 0.0)
tt = np.array(times); R = np.array(R)
m = tt <= 160
plt.figure(figsize=(7, 5))
plt.plot(tt[m], R[m]**2, "o-", label=r"measured $R^2$")
plt.plot(tt[m], np.maximum(625 - 4*tt[m], 0), "--",
         label=r"theory $R^2 = R_0^2 - 4L\kappa t$")
plt.xlabel("time"); plt.ylabel(r"$R^2$")
plt.title("Curvature-driven shrinkage")
plt.legend(); plt.grid(True, alpha=0.3)
plt.savefig("ac2d_radius.png", dpi=150, bbox_inches="tight")
plt.close()

print("wrote ac2d_panels.png, ac2d_slice.png, ac2d_radius.png")
