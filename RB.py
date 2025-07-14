import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const

# Parameters
rho = 8.7e-27  # Critical density (kg/m^3)
rho_early = 100  # Adjusted early density (kg/m^3)
delta_V = 0.1 * rho_early * const.c.value**2  # Vacuum energy release (J/m^3)

# Bubble radius
R_b = np.sqrt(3 * delta_V / (8 * np.pi * const.G.value * rho))

# Error propagation (assuming 10% uncertainty in rho, rho_early)
sigma_rho = 0.1 * rho
sigma_rho_early = 0.1 * rho_early
sigma_delta_V = 0.1 * delta_V  # From rho_early
sigma_R_b = R_b * 0.5 * np.sqrt((sigma_delta_V / delta_V)**2 + (sigma_rho / rho)**2)

print(f"Bubble radius R_b: {R_b:.2e} m ± {sigma_R_b:.2e} m")
print(f"Matches target 4.3e26 m: {abs(R_b - 4.3e26) / 4.3e26 * 100:.2f}% error")

# Plot R_b vs rho_early to show sensitivity
rho_early_vals = np.logspace(0, 17, 100)
R_b_vals = [np.sqrt(3 * 0.1 * rho_e * const.c.value**2 / (8 * np.pi * const.G.value * rho)) for rho_e in rho_early_vals]

plt.loglog(rho_early_vals, R_b_vals, label='R_b')
plt.axhline(4.3e26, color='r', linestyle='--', label='Target R_b = 4.3e26 m')
plt.axvline(100, color='g', linestyle='--', label='rho_early = 100 kg/m^3')
plt.xlabel('Early Density ρ_early (kg/m^3)')
plt.ylabel('Bubble Radius R_b (m)')
plt.title('Bubble Radius vs Early Density in AFCT')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.savefig('bubble_radius_plot.png')
plt.show()
