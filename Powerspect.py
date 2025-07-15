import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Constants
hbar = 1.055e-34  # J s
m = 1.78e-58      # kg
c = 3e8           # m/s
k_B = 1.381e-23   # J/K
v = 9.89e15       # m^-3/2
g = 7.87e-11      # m^3/s^2
R_b = 4.3e26      # m
rho_rad = 4.66e-31  # kg/m^3
c_s = 1.38e7      # m/s
r_s = 148e6 * 3.086e16  # m
lambda_damp = 0.07 * 3.086e16  # m

# GP equation
def gp_equation(psi, x, k, g, m, hbar):
    psi_r, psi_i = psi
    n = m * (psi_r**2 + psi_i**2)
    return [(hbar**2/(2*m) * k**2 - g * n) * psi_i / hbar,
            -(hbar**2/(2*m) * k**2 - g * n) * psi_r / hbar]

# Simulate
x = np.linspace(0, R_b, 1000)
psi0 = [v/np.sqrt(2), 0]
k_values = np.linspace(np.pi/r_s, 15*np.pi/r_s, 150)  # Up to l~1500
C_l = []
l_values = []
for k in k_values:
    sol = odeint(gp_equation, psi0, x, args=(k, g, m, hbar))
    psi = sol[:, 0] + 1j * sol[:, 1]
    delta = np.abs(psi) / (v/np.sqrt(2)) - 1
    l = k * R_b * (148e6 * 3.086e16 / R_b) * (np.pi / 0.01062)
    C_l.append(np.mean(delta**2) * np.cos(k * r_s) * np.exp(-k**2 * lambda_damp**2))
    l_values.append(l)

# Plot
plt.plot(l_values, C_l, label='CMB Power Spectrum')
plt.axvline(x=220, color='r', linestyle='--', label='l_1=220')
plt.axvline(x=540, color='g', linestyle='--', label='l_2=540')
plt.axvline(x=840, color='b', linestyle='--', label='l_3=840')
plt.xlabel('Multipole l')
plt.ylabel('C_l (arbitrary units)')
plt.legend()
plt.show()
