import numpy as np
from matplotlib import pyplot as plt
from astropy.cosmology import Planck15
import astropy.units as u
from scipy.interpolate import interp1d

data = np.genfromtxt("../output/DM_output.txt", delimiter="\t")

# 200 times total critical density at z=0
rho200 = 200 * Planck15.critical_density0

r = data[:, 0][1:]
m = data[:, 1][1:]
v = data[:, 2][1:]
rho = data[:, 3][1:]
rho_mean = data[:, 4][1:]

# Rho is in units of M_sun/kpc^3 and needs to be converted to g/cm^3 to match Planck15.critical_density0
conversion_factor = (u.M_sun / u.kpc ** 3).to(u.g / u.cm ** 3)
rho = rho * conversion_factor
rho_mean = rho_mean * conversion_factor

plt.plot(r, rho, color="crimson", zorder=5)
plt.axhline(rho200.value, linestyle="--", color="grey", zorder=1)
plt.yscale("log")
plt.xscale("log")
plt.title("Plot representing task 4 and 5")
plt.xlabel("Radius $r$ [kpc]")
plt.ylabel(r"Density $\rho$ $\left[\frac{g}{cm^3}\right]$")
plt.legend([r"Density $\rho$", r"$200\rho_{crit}$"])
plt.tight_layout()
plt.savefig("task45.pdf")
plt.show()

cumulative_mass = np.array([np.sum(m[:i + 1]) for i in range(len(m))])
cumulative_mass = np.concatenate([np.array([0]), cumulative_mass])
r = np.concatenate([np.array([0]), r])

plt.plot(r, cumulative_mass/(4/3*np.pi*np.power(r, 3)), color="navy", zorder=5)
plt.title("Plot representing task 6")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Radius $r$ [kpc]")
plt.ylabel(r"Cumulative mass inside radius $M_{tot}$ [$M_{\odot}$]")
plt.tight_layout()
plt.savefig("task6.pdf")
plt.show()

# Linearly interpolate for R200
interp_func = interp1d(rho_mean, r[1:])
R200 = interp_func(rho200)
interp_func = interp1d(r, cumulative_mass)
M200 = interp_func(R200)

print("R200 =", np.round(R200, 2), "kpc")
print('R200 (as calculated by Martin Sparre)', r[np.where(rho_mean > rho200.value)][-1])
print("M200 =", "{:.2E}".format(M200), "M_sol")
