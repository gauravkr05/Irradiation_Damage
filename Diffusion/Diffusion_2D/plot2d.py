import numpy as np
import matplotlib.pyplot as plt

save_times = [0, 200, 500, 1000]

for tm in save_times:
    folder = "time_%d" % tm
    data = np.loadtxt("%s/c_time_%d.dat" % (folder, tm), skiprows=1)

    nx = len(np.unique(data[:, 0]))
    ny = len(np.unique(data[:, 1]))
    c = data[:, 2].reshape(nx, ny)

    plt.figure(figsize=(6, 5))
    plt.imshow(c.T, origin="lower", cmap="viridis", vmin=0.0, vmax=1.0,
               extent=[data[:, 0].min(), data[:, 0].max(),
                       data[:, 1].min(), data[:, 1].max()])
    plt.colorbar(label="c")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("2D Diffusion at time = %d" % tm)
    plt.savefig("%s/c_time_%d.png" % (folder, tm), dpi=150, bbox_inches="tight")
    plt.close()
