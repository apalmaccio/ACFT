import numpy as np
from boson_mass import m

rho_crit = 8.7e-27  # kg/m^3
v = np.sqrt(2 * rho_crit / m)

print(f"VEV v: {v:.2e} m^-3/2")
