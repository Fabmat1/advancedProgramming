import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

galaxypos = [
    (0., 0., 0.),
    (28.9307, 147.77176, 91.3073),
    (159.32645, -27.599352, -80.43161),
    (-29.184427, 266.3879, -95.67256),
    (-10.797834, 148.63385, -17.152555),
    (-44.74542, 58.259293, -40.10336),
    (-82.00299, 89.14701, -168.42316),
    (39.252796, 20.820791, 19.300589),
    (162.94783, -146.49446, 32.405037),
    (170.69516, -58.414986, -39.244144)
]
galaxypos = np.array(galaxypos)

# Task 1 & 2
arr = np.loadtxt(
    r"GalaxyFromIllustrisTNG50Dark_DM_Subhalo852966.txt"
)

pos = arr[:, 0:3]
m = arr[:, 3]

fig, ax = plt.subplots(1, 3, figsize=(15, 5))

h0 = ax[0].hist2d(pos[:, 0], pos[:, 1], weights=m, bins=600, cmin=0, range=np.array([(-300, 300), (-300, 300)]),
                  cmap="viridis")
h1 = ax[1].hist2d(pos[:, 0], pos[:, 2], weights=m, bins=600, cmin=0, range=np.array([(-300, 300), (-300, 300)]),
                  cmap="viridis")
h2 = ax[2].hist2d(pos[:, 1], pos[:, 2], weights=m, bins=600, cmin=0, range=np.array([(-300, 300), (-300, 300)]),
                  cmap="viridis")

ax[0].scatter(galaxypos[:, 0], galaxypos[:, 1], marker="x", c="red")
ax[1].scatter(galaxypos[:, 0], galaxypos[:, 2], marker="x", c="red")
ax[2].scatter(galaxypos[:, 1], galaxypos[:, 2], marker="x", c="red")

fig.suptitle("Task 1 & 2: 2d slices of the dark matter density")
ax[0].set_xlabel("x")
ax[1].set_xlabel("x")
ax[2].set_xlabel("y")
ax[0].set_ylabel("y")
ax[1].set_ylabel("z")
ax[2].set_ylabel("z")


cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])

# Sqrt normalization
norm = mpl.colors.PowerNorm(1/2, vmin=np.min(np.concatenate([h0[0].flatten(), h1[0].flatten(), h2[0].flatten()])),
                          vmax=np.max(np.concatenate([h0[0].flatten(), h1[0].flatten(), h2[0].flatten()])))

h0[3].norm = norm
h1[3].norm = norm
h2[3].norm = norm

for a in ax:
    a.set_aspect("equal", adjustable='box')
    a.set_xlim(-300, 300)
    a.set_ylim(-300, 300)


fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="viridis"),
             cax=cax, label='$M_{\odot}$kpc$^{-2}$')

plt.show()
