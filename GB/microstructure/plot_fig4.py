# ============================================================
#  plot_fig4.py  —  run with:  python plot_fig4.py
#  Needs: pip install numpy matplotlib
#  Input : fig4_out/phi_map.dat   (written by fig4_microstructure.cpp)
#          columns: x  y  phi  dominant_orientation
#  Output: fig4_out/fig4_microstructure.png
#          fig4_out/fig4_orientation.png
# ============================================================
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("fig4_out/phi_map.dat")          # blank lines are ignored by loadtxt
n = int(round(np.sqrt(len(data))))                 # square grid: 200 x 200

# the file is written x-major (all y for one x, then a blank line),
# so reshape to (nx, ny) and transpose to get rows = y, columns = x
phi    = data[:, 2].reshape(n, n).T
orient = data[:, 3].reshape(n, n).T

extent = [data[:, 0].min(), data[:, 0].max(),
          data[:, 1].min(), data[:, 1].max()]

# ---- grey-level microstructure: bright grains, dark boundary lines ----
plt.figure(figsize=(6.5, 5.8))
im = plt.imshow(phi, cmap='gray', origin='lower', extent=extent, vmin=0.6, vmax=1.0)
plt.colorbar(im, label=r'$\phi$')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Microstructure: two straight grain boundaries')
plt.tight_layout()
plt.savefig("fig4_out/fig4_microstructure.png", dpi=150)

# ---- same cell coloured by which orientation dominates ----
plt.figure(figsize=(6.5, 5.8))
im = plt.imshow(orient, cmap='coolwarm', origin='lower', extent=extent)
plt.colorbar(im, label='orientation index', ticks=[0, 1])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Dominant orientation (grain 1 vs grain 2)')
plt.tight_layout()
plt.savefig("fig4_out/fig4_orientation.png", dpi=150)
plt.show()