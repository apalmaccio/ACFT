import numpy as np
from scipy.integrate import solve_ivp

def H_eff(z, H0=2.6e-25, c=0.12, Tc=1.2e12):
    phi = 0.1 / (1 + z)
    T = 2.73 * (1 + z)
    phase = 1 + 1.67e-37 * (T/Tc)**3 * np.exp(-T/Tc)
    dphi_dz = -0.1 / (1 + z)**2
    dphase_dz = 1.67e-37 * np.exp(-T/Tc) * (3 * (T/Tc)**2 * 2.73/Tc - (T/Tc)**3 * 2.73/Tc)
    return H0 * c * (dphi_dz * phase + phi * dphase_dz) / (1 + c * phi * phase)

def growth_ode(z, y):
    delta, ddelta_dz = y
    Omega_m = (1 + z)**3 / H_eff(z)**2
    return [ddelta_dz, -ddelta_dz/(1 + z) + 1.5 * Omega_m * delta]

z = np.array([0.5, 0.8, 1.0, 1.5, 2.0])
f_sigma8 = []
sigma8_z0 = 0.81
for z_target in z:
    sol = solve_ivp(growth_ode, [0, z_target], [1e-5, 0], method='Radau', rtol=1e-10)
    f = -sol.y[1][-1] / sol.y[0][-1] / (1 + z_target)
    sigma8_z = sigma8_z0 * sol.y[0][-1] / sol.y[0][0]
    f_sigma8.append(f * sigma8_z)

print(f"f_sigma8 = {f_sigma8}")
