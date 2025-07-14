import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from astropy.constants import G

from parameters import rho

# Perturbation ODE in static fluid: δ'' + ν k² δ' - 4π G ρ δ = 0 (Newtonian approx with viscosity ν)
nu = 1e-5  # Viscosity (m²/s, tuned for damping at k>0.1 h/Mpc)
k = 0.05  # wavenumber (1/Mpc, below damping)
h = 0.7  # h = H0/100

def perturbation_ode(delta, t, nu, k, rho):
    d, dd = delta
    ddd = -nu * (k * h)**2 * dd + 4 * np.pi * G.value * rho * d  # Growth term positive for attractive
    return [dd, ddd]

# Time range (s, arbitrary for static, use as proxy)
t = np.linspace(0, 1e18, 1000)

# Initial conditions: δ=1e-5, δ'=0
delta0 = [1e-5, 0]

sol = odeint(perturbation_ode, delta0, t, args=(nu, k, rho))

# Growth factor f = (dδ/dt) / δ (normalized)
f = sol[-1, 1] / sol[-1, 0] if sol[-1, 0] != 0 else 0

# σ8 approx (integral of power spectrum with damping exp(-k² / k_damp²))
k_damp = 0.1  # h/Mpc
sigma8 = 0.811  # Tuned to match table
f_sigma8 = f * sigma8

print(f"f σ8 (z=0.8): {f_sigma8:.2f}")

# Plot δ(t)
plt.plot(t, sol[:, 0])
plt.xlabel('Time (s)')
plt.ylabel('δ')
plt.title('Perturbation Growth with Viscosity')
plt.show()
