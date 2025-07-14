import numpy as np
from astropy import constants as const
from scattering_length import a_s, m
from vev import v

# Non-relativistic BEC T_c
n = 8.7e-27 / m  # m^-3
T_c = (const.hbar.value**2 * n**(2/3)) / (const.k_B.value * m)

print(f"Critical temperature T_c: {T_c:.2e} K")
