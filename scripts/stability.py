import numpy as np
from scipy.integrate import solve_ivp

def perturb_ode(k, t, y, g_phi=4.1e-19, G=6.674e-11):
    delta_phi, ddelta_phi = y
    delta_rho_m = 1e-5  # Perturbation amplitude
    return [ddelta_phi, g_phi * delta_rho_m - k**2 * delta_phi]

k_values = np.logspace(-2, 0, 10)  # k = 0.01 to 1 h/Mpc
max_pert = []
for k in k_values:
    sol = solve_ivp(perturb_ode, [0, 1e10], [0, 0], args=(k,), method='Radau', rtol=1e-10)
    max_pert.append(max(abs(sol.y[0])))
print(f"Max perturbation = {max(max_pert)}")
