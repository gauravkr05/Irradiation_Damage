import numpy as np
import matplotlib.pyplot as plt

save_times = [0, 200, 500, 1000]

plt.figure(figsize=(7, 5))
for tm in save_times:
    data = np.loadtxt("time_%d/c_time_%d.dat" % (tm, tm), skiprows=1)
    plt.plot(data[:, 0], data[:, 1], label="time = %d" % tm)

plt.xlabel("x")
plt.ylabel("c")
plt.title("1D Diffusion")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("diffusion_profiles.png", dpi=150, bbox_inches="tight")
plt.show()
import numpy as np
import matplotlib.pyplot as plt

save_times = [0, 200, 500, 1000]

for tm in save_times:
    folder = "time_%d" % tm
    data = np.loadtxt("%s/c_time_%d.dat" % (folder, tm), skiprows=1)

    plt.figure(figsize=(7, 5))
    plt.plot(data[:, 0], data[:, 1])
    plt.xlabel("x")
    plt.ylabel("c")
    plt.title("1D Diffusion at time = %d" % tm)
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    plt.savefig("%s/c_time_%d.png" % (folder, tm), dpi=150, bbox_inches="tight")
    plt.close()