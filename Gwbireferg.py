import numpy as np
import matplotlib.pyplot as plt
from parameters import omega  # Vorticity ~10^{-18} s^{-1}

# GW strain with birefringence
def gw_strain(f, t, asymmetry=0.01):
    h_plus = np.sin(2 * np.pi * f * t)
    h_cross = np.cos(2 * np.pi * f * t) * (1 + asymmetry)  # Vorticity effect
    return h_plus, h_cross

# Frequencies for stochastic background ~10^{-8} Hz
f_bg = 1e-8
t = np.linspace(0, 10, 1000)

h_plus, h_cross = gw_strain(f_bg, t)

# Plot polarization
plt.plot(t, h_plus, label='h+')
plt.plot(t, h_cross, label='h× (with 1% asymmetry)')
plt.xlabel('Time (s)')
plt.ylabel('Strain')
plt.title('GW Birefringence in AFCT')
plt.legend()
plt.show()

# Burst Omega_GW h^2 ~10^{-10}
Omega_gw = 1e-10 * (omega / 1e-18)**2  # Scaled by vorticity
print(f"GW burst Ω_GW h²: {Omega_gw:.2e}")
