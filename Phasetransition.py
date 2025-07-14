import numpy as np
from scipy.integrate import quad
from parameters import T_c, lambda_val, gamma, v  # Add lambda_val = -1/f_a^2 ≈ -1e-24, gamma>0≈1e-50

# Potential V(Φ) from Lagrangian
def potential(phi):
    return m**2 * phi**2 + (lambda_val / 2) * phi**4 + (gamma / 6) * phi**6

# Critical density at transition (minima)
def critical_density(T):
    if T > T_c:
        return 0  # No condensate
    else:
        # Approximate from mean-field: ρ_c ≈ -3 m^4 / (4 lambda) for Mexican hat
        return m * v**2 / 2 * (1 - (T / T_c)**2)  # Simple Landau approximation

# Bubble nucleation rate (supercooling)
def nucleation_rate(T):
    # Action S ≈ (ΔV)^{1/2} / T^2 or similar; placeholder for GW burst
    delta_V = potential(0) - min(potential(np.linspace(0, v, 100)))  # False to true vacuum
    return np.exp(-delta_V / T**4)  # Euclidean action approx

# Integrate over temperature for transition
T_vals = np.linspace(1e9, T_c + 1e9, 100)
rho_vals = [critical_density(T) for T in T_vals]

# Plot
import matplotlib.pyplot as plt
plt.plot(T_vals, rho_vals)
plt.axvline(T_c, color='r', label='T_c')
plt.xlabel('Temperature (K)')
plt.ylabel('Condensate Density (kg/m³)')
plt.title('BEC Phase Transition in AFCT')
plt.legend()
plt.show()

# GW burst frequency from transition ~10^{-3} Hz
f_gw = 1e-3 / (T_c / 1e10)  # Scaled
print(f"Phase transition GW frequency: {f_gw:.2e} Hz")
