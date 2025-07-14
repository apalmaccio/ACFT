import numpy as np
from astropy import constants as const

# De Broglie wavelength for galactic core
lambda_dB = 1e3 * 3.08568e16  # 1 kpc in m
v_gal = 2e5  # 200 km/s
hbar = const.hbar.value
m = hbar / (lambda_dB * v_gal)  # kg
m_eV = m * const.c.value**2 / 1.602e-19  # Convert to eV

print(f"Boson mass: {m:.2e} kg ({m_eV:.2e} eV)")
