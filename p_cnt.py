import matplotlib.pyplot as plt
import numpy as np


data = np.loadtxt("total_counts.txt", usecols=2)
data = data.reshape((-1, 32))

fig, ax = plt.subplots()
im = ax.pcolor(data)

# fig.colorbar(im, ax=ax, label="Interactive colorbar")

plt.show()
