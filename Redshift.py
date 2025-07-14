import numpy as np
from scipy.optimize import fsolve
from astropy.constants import c
from parameters import alpha, omega, sigma_scatter_si

# Full 1 + z ≈ n(z) * γ_d * (1 + σ z)
# n(z) = 1 + α exp(-z/(1+z))
# γ_d = exp(ω χ / (2 c))  # Unit-adjusted integral

def z_from_chi(chi_val):
    def equation(z):
        n_z = 1 + alpha * np.exp(-z / (1 + z))
        gamma_d = np.exp(omega * chi_val / (2 * c.value))
        return 1 + z - n_z * gamma_d * (1 + sigma_scatter_si * chi_val)  # σ z ≈ σ χ (approx)
    z_guess = 1  # Initial guess
    return fsolve(equation, z_guess)[0]

def chi_from_z(z_val):
    def equation(chi):
        return z_from_chi(chi) - z_val
    chi_guess = 1e26  # m
    return fsolve(equation, chi_guess)[0]

# Example: χ at z=1089 (matches table 4.29e26 m)
z_cmb = 1089
chi_cmb = chi_from_z(z_cmb)
print(f"χ(z=1089): {chi_cmb:.2e} m")

# Plot z vs χ
chi_vals = np.linspace(1e25, 5e26, 100)
z_vals = [z_from_chi(chi) for chi in chi_vals]
import matplotlib.pyplot as plt
plt.plot(chi_vals, z_vals)
plt.xlabel('χ (m)')
plt.ylabel('z')
plt.title('Redshift in AFCT')
plt.show()
