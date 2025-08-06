import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("nsx_output.txt")
E_obs = data[:, 0]
Insx = data[:, 1]
Ibb = data[:, 2]
plt.figure(figsize=(8, 8))

plt.plot(E_obs, Ibb / E_obs, "k--", lw=2.0)
plt.plot(E_obs, Insx / E_obs, "b-", lw=2.0)

ax = plt.gca()
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylim([9.0e25,4.0e29])
ax.set_ylabel("Photon specific intensity")
_ = ax.set_xlabel("Energy [keV]")
# plt.savefig('nsx_intensity_plot.png')
plt.show()
