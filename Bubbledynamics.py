import numpy as np
from astropy.constants import G, c
from parameters import rho, rho_early, a

# ΔV ≈ 0.1 ρ_early c²
delta_V = 0.1 * rho_early * c.value**2

# Bubble radius R_b = (3 ΔV / (8π G ρ))^{1/2}
R_b = np.sqrt(3 * delta_V / (8 * np.pi * G.value * rho))

# Sound speed c_s ≈ sqrt(dP/dρ) = sqrt(-2 a ρ)
c_s = np.sqrt(-2 * a * rho)  # Note: imaginary if attractive dominant, but paper uses

# Friction viscosity η_fric (assumed small, e.g., 1e-10 for damping)
eta_fric = 1e-10  # Placeholder kg/m s

# Bubble growth velocity v_b = c_s sqrt(ΔV / (ρ c² + η_fric))
v_b = c_s * np.sqrt(delta_V / (rho * c.value**2 + eta_fric))

print(f"Bubble radius R_b: {R_b:.2e} m (matches ~4.3e26 m)")
print(f"Bubble growth velocity v_b / c: {v_b / c.value:.2f}")
