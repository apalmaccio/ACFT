import numpy as np
import matplotlib.pyplot as plt
from parameters import R_b, sigma_scatter_si  # Bubble radius and scattering

# Intensity I(r) ~ exp(-σ r) / r^2 for finite sources
def intensity(r):
    if r > R_b:
        return 0  # Beyond bubble
    return np.exp(-sigma_scatter_si * r) / r**2  # Dimming

# Integrate over radius for total sky brightness
r_vals = np.linspace(1e25, R_b, 100)
I_vals = [intensity(r) for r in r_vals]

# Plot
plt.loglog(r_vals, I_vals)
plt.axvline(R_b, color='r', label='Bubble Edge')
plt.xlabel('Distance (m)')
plt.ylabel('Relative Intensity')
plt.title("Olbers' Paradox Resolution in AFCT")
plt.legend()
plt.show()

# Total brightness (finite integral)
total_I = np.trapz(I_vals, r_vals)
print(f"Total sky brightness (arbitrary units): {total_I:.2e} (finite, resolving paradox)")
