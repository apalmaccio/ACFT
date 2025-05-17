import numpy as np
from scipy.integrate import odeint

def dR_dT(R, T, phi0=1.22e19, Tc=1.2e12, lambda_phi=1.67e-37):
    phi = phi0 * (T / Tc)**-1
    sigma = lambda_phi * phi0**4 * np.exp(-T / Tc)
    dsigma_dT = -sigma / Tc
    num = R**2 * (dsigma_dT + lambda_phi * phi**4)
    denom = R * (2 * sigma + R * 0.5 * lambda_phi * phi**4)
    return -num / denom

T = np.linspace(1.2e12, 1e9, 1000)
R = odeint(dR_dT, 1e4, T)
R_f = np.trapz(R[:, 0], T) * 3e8  # Convert to meters
print(f"R_f = {R_f} m")
