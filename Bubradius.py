import numpy as np
from astropy import constants as const

rho = 8.7e-27
rho_early = 1e17
delta_V = 0.1 * rho_early * const.c.value**2
R_b = np.sqrt(3 * delta_V / (8 * np.pi * const.G.value * rho))

print(f"Bubble radius R_b: {R_b:.2e} m")
