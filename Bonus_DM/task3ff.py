import numpy as np
from matplotlib import pyplot as plt
from astropy.cosmology import Planck15
import astropy.units as u

data = np.genfromtxt("../output/DM_output.txt", delimiter="\t")

# 200 times total critical density at z=0
R200 = 200*Planck15.critical_density0

r = data[:, 0][1:]
m = data[:, 1][1:]
v = data[:, 2][1:]
rho = data[:, 3][1:]

# Rho is in units of M_sun/kpc^3 and needs to be converted to g/cm^3 to match Planck15.critical_density0
conversion_factor = (u.M_sun / u.kpc**3).to(u.g / u.cm**3)
rho = rho * conversion_factor


plt.plot(r, rho, color="navy")
plt.yscale("log")
plt.xscale("log")
plt.title("Plot for Task 4")
plt.xlabel("Radius $r$ [kpc]")
plt.ylabel(r"Density $\rho$ $\left[\frac{M_{\odot}}{kpc^3}\right]$")
plt.tight_layout()
plt.show()

plt.plot(r, rho, color="crimson", zorder=5)
plt.axhline(R200.value, linestyle="--", color="grey", zorder=1)
plt.yscale("log")
plt.xscale("log")
plt.title("Plot representing task 5 and 6")
plt.xlabel("Radius $r$ [kpc]")
plt.ylabel(r"Density $\rho$ $\left[\frac{g}{cm^3}\right]$")
plt.legend([r"Density $\rho$", r"$200\rho_{crit}$"])
plt.tight_layout()
plt.show()

