import numpy as np

kappa_phi0 = 0.12
beta = 1.67e-37
T0 = 2.73
Tc = 1.2e12
z = 1e4
phase_term = 1 + beta * (T0 * (1 + z) / Tc)**3 * np.exp(-T0 * (1 + z) / Tc)
n_z = 1 + kappa_phi0 / (1 + z) * phase_term
v_phase = 3e8 / n_z
print(f"v_phase = {v_phase} m/s")
