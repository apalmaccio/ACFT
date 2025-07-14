import numpy as np
from boson_mass import m
from astropy import constants as const

P_DE = 1e-9  # Pa
rho = 8.7e-27  # kg/m^3
a = P_DE / rho**2  # m^5 kg^-1 s^-2
a_s = a * m**3 / (2 * np.pi * const.hbar.value**2)

print(f"Scattering length |a_s|: {a_s:.2e} m")
print(f"a: {a:.2e} m^5 kg^-1 s^-2")
