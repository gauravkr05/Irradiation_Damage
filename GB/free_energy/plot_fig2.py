# ============================================================
#  plot_fig2.py  —  run with:  python plot_fig2.py
#  Needs: pip install numpy matplotlib
#  Input : fig2_out/free_energy.dat   (written by fig2_free_energy.cpp)
#  Output: fig2_out/fig2_free_energy.png
# ============================================================
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("fig2_out/free_energy.dat")      # columns: x  f  f-fmin
x, f = data[:, 0], data[:, 1]

fmin = -0.25                # -alpha/2 + beta/4 for alpha = beta = 1

plt.figure(figsize=(7, 5))
plt.plot(x, f, 'o-', ms=4, lw=1.5, color='darkred', label='f')
plt.axhline(fmin, ls='--', lw=0.8, color='black', label=r'$f_{min} = -0.25$')

plt.xlim(160, 240)
plt.xlabel('x')
plt.ylabel(r'$f(\eta_1,\eta_2)$')
plt.title('Local free energy density across the grain boundary')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("fig2_out/fig2_free_energy.png", dpi=150)
plt.show()