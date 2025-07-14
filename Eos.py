import numpy as np
import matplotlib.pyplot as plt
from parameters import rho, a, gamma_prime  # Import from parameters.py

# EOS: P = -a ρ² + (γ'/2) ρ³
def eos_pressure(rho_val):
    return -a * rho_val**2 + (gamma_prime / 2) * rho_val**3

# Density range around cosmic ρ
rho_vals = np.logspace(-30, -20, 100)  # kg/m³

P_vals = eos_pressure(rho_vals)

# Plot P vs ρ
plt.figure()
plt.loglog(rho_vals, np.abs(P_vals))  # Absolute for visibility, as P negative
plt.xlabel('Density ρ (kg/m³)')
plt.ylabel('|Pressure P| (Pa)')
plt.title('Van der Waals EOS in AFCT')
plt.axvline(rho, color='r', linestyle='--', label='Cosmic ρ')
plt.legend()
plt.savefig('eos_plot.png')
plt.show()

# Example at cosmic ρ (P ≈ -10^{-9} Pa)
print(f"P at ρ = {rho:.2e}: {eos_pressure(rho):.2e} Pa")
