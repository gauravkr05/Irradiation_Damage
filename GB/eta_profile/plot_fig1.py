# ============================================================
#  plot_fig1.py  —  run with:  python plot_fig1.py
#  Needs: pip install numpy matplotlib
#  Input : fig1_out/eta_profiles.dat   (written by fig1_eta_profiles.cpp)
#  Output: fig1_out/fig1_eta_profiles.png  (and shows the plot on screen)
# ============================================================
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("fig1_out/eta_profiles.dat")     # columns: x  eta1  eta2
x, eta1, eta2 = data[:, 0], data[:, 1], data[:, 2]

plt.figure(figsize=(7, 5))
plt.plot(x, eta1, 'o-', ms=4, lw=1.5, label=r'$\eta_1$')
plt.plot(x, eta2, 's-', ms=4, lw=1.5, label=r'$\eta_2$')

plt.xlim(160, 240)          # zoom on the middle boundary; full cell is 0 to 400
plt.ylim(-0.1, 1.1)
plt.xlabel('x')
plt.ylabel(r'$\eta_i$')
plt.title('Equilibrium profiles of the field variables across a straight GB')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("fig1_out/fig1_eta_profiles.png", dpi=150)
plt.show()