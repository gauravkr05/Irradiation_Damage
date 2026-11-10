# ============================================================
#  plot_fig3.py  —  run with:  python plot_fig3.py
#  Needs: pip install numpy matplotlib
#  Input : fig3_out/phi_profile.dat   (written by fig3_phi_profile.cpp)
#  Output: fig3_out/fig3_phi.png
# ============================================================
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("fig3_out/phi_profile.dat")      # columns: x  phi
x, phi = data[:, 0], data[:, 1]

plt.figure(figsize=(7, 5))
plt.plot(x, phi, 'o-', ms=4, lw=1.5, color='darkgreen', label=r'$\phi(r)$')

plt.xlim(160, 240)
plt.ylim(0.65, 1.05)
plt.xlabel('x')
plt.ylabel(r'$\phi = \sum_i \eta_i^2$')
plt.title(r'Visualisation function $\phi(r)$ across the boundary')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("fig3_out/fig3_phi.png", dpi=150)
plt.show()