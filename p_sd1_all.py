import numpy as np
import matplotlib.pyplot as plt
import sys

plt.figure(figsize=(12, 7))

for sid, s in enumerate(["a", "b", "c", "d", "e", "f"]):
    dat1 = np.loadtxt(f"sd1{s}.txt", skiprows=2)
    dat2 = np.loadtxt(f"apjlab5968/SD1{s}_test_IM.txt")

    plt.subplot(2, 3, sid + 1)

    plt.plot(dat2[:, 0], dat1, ".", ms=1)
    plt.plot(dat2[:, 0], dat2[:, 1], ".", ms=1)

    if True:
        yhat = dat1
        y = dat2[:, 1]
        rel_err = (yhat - y) / np.median(y) * 1e3
        plt.twinx()
        plt.plot(dat2[:, 0], rel_err, "r--")
        plt.ylim(-5, 5)
        plt.hlines([-1, 0, 1], 0, 1, color="r", lw=0.5)
        plt.yticks(color="r")
    plt.title(f"SD1{s}")

plt.tight_layout()
if sys.argv[-1] == "s":
    plt.savefig("p_sd1_all.png", dpi=400)
    print("Saved ./p_sd1_all.png")
else:
    plt.show()
