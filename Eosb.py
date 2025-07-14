import numpy as np
from scattering_length import a_s, m

b = 16 * np.pi * a_s**3 / (3 * m)

print(f"b: {b:.2e} m^3 kg^-1")
print(f"Stability check b*rho: {b * 8.7e-27:.2e}")
