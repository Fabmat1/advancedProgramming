import copy
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from astropy.cosmology import Planck15
import astropy.units as u


def plot_galaxy_2d(pos, m, saveas, figtitle, do_lookback_slice=False, lkbk_times=None, lkbk_range=None):
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))

    # Set "bad pixels", meaning pixels outside of the range of the normalization to be the 0 color of the colormap
    my_cmap = copy.copy(mpl.cm.get_cmap('viridis'))
    my_cmap.set_bad(my_cmap.colors[0])

    if do_lookback_slice:
        mask = np.logical_and(lkbk_times > lkbk_range[0], lkbk_times < lkbk_range[1])
        m = m[mask]
        h0 = ax[0].hist2d(pos[:, 0][mask], pos[:, 1][mask], weights=m, bins=600, cmin=0,
                          range=np.array([(-25, 25), (-25, 25)]),
                          cmap=my_cmap)
        h1 = ax[1].hist2d(pos[:, 0][mask], pos[:, 2][mask], weights=m, bins=600, cmin=0,
                          range=np.array([(-25, 25), (-25, 25)]),
                          cmap=my_cmap)
        h2 = ax[2].hist2d(pos[:, 1][mask], pos[:, 2][mask], weights=m, bins=600, cmin=0,
                          range=np.array([(-25, 25), (-25, 25)]),
                          cmap=my_cmap)

    else:
        h0 = ax[0].hist2d(pos[:, 0], pos[:, 1], weights=m, bins=600, cmin=0, range=np.array([(-25, 25), (-25, 25)]),
                          cmap=my_cmap)
        h1 = ax[1].hist2d(pos[:, 0], pos[:, 2], weights=m, bins=600, cmin=0, range=np.array([(-25, 25), (-25, 25)]),
                          cmap=my_cmap)
        h2 = ax[2].hist2d(pos[:, 1], pos[:, 2], weights=m, bins=600, cmin=0, range=np.array([(-25, 25), (-25, 25)]),
                          cmap=my_cmap)

    fig.suptitle(figtitle)
    ax[0].set_xlabel("x")
    ax[1].set_xlabel("x")
    ax[2].set_xlabel("y")
    ax[0].set_ylabel("y")
    ax[1].set_ylabel("z")
    ax[2].set_ylabel("z")
    ax[0].set_aspect("equal", adjustable='box')
    ax[1].set_aspect("equal", adjustable='box')
    ax[2].set_aspect("equal", adjustable='box')

    cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])

    # Log normalization
    normmatrix = np.concatenate([h0[0].flatten(), h1[0].flatten(), h2[0].flatten()])
    norm = mpl.colors.LogNorm(vmin=np.min(normmatrix[normmatrix != 0]),
                              vmax=np.max(normmatrix))

    # A loop. So space-efficient.
    for h in [h0, h1, h2]:
        h[3].norm = norm

    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="viridis"), cax=cax, label=r'$M_{\odot}$kpc$^{-2}$')

    plt.savefig(saveas, dpi=500)
    plt.close()  # plt.show()


def pol_cyl_hist(pos, m, saveas, title, do_lookback_slice=False, lkbk_times=None, lkbk_range=None):
    if do_lookback_slice:
        mask = np.logical_and(lkbk_times > lkbk_range[0], lkbk_times < lkbk_range[1])

        x = pos[:, 0][mask]
        y = pos[:, 1][mask]
        z = pos[:, 2][mask]

        radial_dist = np.sqrt(y ** 2 + z ** 2)

        hist, bin_edges = np.histogram(radial_dist, bins=20, weights=m[mask])
    else:
        x = pos[:, 0]
        y = pos[:, 1]
        z = pos[:, 2]

        radial_dist = np.sqrt(y ** 2 + z ** 2)

        hist, bin_edges = np.histogram(radial_dist, bins=20, weights=m)

    hist /= 2 * np.pi * (bin_edges[1:] + bin_edges[:-1]) / 2 * (bin_edges[1] - bin_edges[0])

    plt.figure(figsize=(4.8*16/9, 4.8))
    plt.bar(bin_edges[:-1], hist, width=np.diff(bin_edges), align='edge', color="salmon", edgecolor="darkred")
    plt.title(title)
    plt.xlabel("Radius $r$ [kpc]")
    plt.ylabel("Surface density [$M_{\odot}$kpc$^{-2}$]")
    plt.semilogy()
    plt.tight_layout()
    plt.savefig(saveas, dpi=300)
    plt.close()  # plt.show()


def plot_SFR(lkbk_time, init_mass, saveas, n_bins):
    dt = np.ptp(lkbk_time.value) / n_bins

    plt.figure(figsize=(4.8*16/9, 4.8))
    plt.hist(lkbk_time.value, weights=init_mass / (dt * u.Gyr).to(u.yr), bins=n_bins, color="lightblue",
             edgecolor="navy")
    plt.title("Star formation history for the simulated galaxy")
    plt.ylabel("SFR [$M_{\odot}$yr$^{-1}$]")
    plt.xlabel("Lookback time [Gyr]")
    plt.xlim(14, 0)
    plt.savefig(saveas, dpi=300)
    plt.close()  # plt.show()


if __name__ == "__main__":
    arr = np.loadtxt(
        r"GalaxyFromIllustrisTNG50_Stars_Subhalo521803.txt"
    )

    pos = arr[:, 0:3]
    m = arr[:, 3]

    plot_galaxy_2d(pos, m, "task_1_plot.png", "Task 1: 2d slices of the stellar density, with non aligned axes.")

    arr = np.loadtxt(
        r"../output/rotated_stars.csv"  # Calculate this before with the C++ script
    )

    pos = arr[:, 0:3]
    m = arr[:, 3]

    plot_galaxy_2d(pos, m, "task_2_plot.png", "Task 2: 2d slices of the stellar density, with aligned axes.")

    pol_cyl_hist(pos, m, "task_5_plot.png", "Polar/cylindrical shell histogram of the stellar mass distribution")

    redshift = arr[:, 5]
    init_mass = arr[:, 4]

    lkbk_time = Planck15.lookback_time(redshift).to(u.Gyr)

    plot_SFR(lkbk_time, init_mass, "Task_7_plot.png", 50)

    plot_galaxy_2d(pos, m, "task_8_plot.png",
                   "Task 8: 2d slices of the stellar density, for lookback times 0Gyr < t < 2Gyr.",
                   True, lkbk_time.value, (0, 2))

    plot_galaxy_2d(pos, m, "task_9.1_plot.png",
                   "Task 9.1: 2d slices of the stellar density, for lookback times 2Gyr < t < 4Gyr.",
                   True, lkbk_time.value, (2, 4))

    plot_galaxy_2d(pos, m, "task_9.2_plot.png",
                   "Task 9.2: 2d slices of the stellar density, for lookback times 4Gyr < t < 6Gyr.",
                   True, lkbk_time.value, (4, 6))

    plot_galaxy_2d(pos, m, "task_9.3_plot.png",
                   "Task 9.3: 2d slices of the stellar density, for lookback times 6Gyr < t < 10Gyr.",
                   True, lkbk_time.value, (6, 10))

    pol_cyl_hist(pos, m, "task_10.1_plot.png", "Polar/cylindrical shell histogram of the stellar mass distribution,\n"
                                               " for lookback times 0Gyr < t < 2Gyr.", True, lkbk_time.value, (0, 2))

    pol_cyl_hist(pos, m, "task_10.2_plot.png", "Polar/cylindrical shell histogram of the stellar mass distribution,\n"
                                               " for lookback times 2Gyr < t < 4Gyr.", True, lkbk_time.value, (2, 4))

    pol_cyl_hist(pos, m, "task_10.3_plot.png", "Polar/cylindrical shell histogram of the stellar mass distribution,\n"
                                               " for lookback times 4Gyr < t < 6Gyr.", True, lkbk_time.value, (4, 6))

    pol_cyl_hist(pos, m, "task_10.4_plot.png", "Polar/cylindrical shell histogram of the stellar mass distribution,\n"
                                               " for lookback times 6Gyr < t < 10Gyr.", True, lkbk_time.value, (6, 10))
