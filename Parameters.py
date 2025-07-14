import numpy as np
from astropy import constants as const

# Boson mass (m ≈ 10^{-22} eV → kg)
m_eV = 1e-22
m = m_eV * 1.60217662e-19 / const.c.value**2  # kg ≈ 1.78e-58

# Condensate density (ρ ≈ 10^{-26} kg/m³)
rho = 1e-26  # kg/m³

# vev v from ρ = m v² / 2
v = np.sqrt(2 * rho / m)  # ≈ 1.06e16 m^{-3/2}

# Scattering length |a_s| ≈ 8.1 × 10^{-64} m
a_s = 8.1e-64  # m

# EOS parameter a = 2π ħ² |a_s| / m³ ≈ 10^{43} m⁵ kg⁻¹ s⁻²
a = 2 * np.pi * const.hbar.value**2 * np.abs(a_s) / m**3

# EOS parameter b ≈ 16π |a_s|³ / (3 m) ≈ 5 × 10^{-182} m³ kg⁻¹
b = 16 * np.pi * a_s**3 / (3 * m)

# γ' ≈ γ / m³ ≈ 10^{129} m⁹ kg⁻² s⁻²
gamma_prime = 1e129

# Critical temperature T_c ≈ 10^{10} K
T_c = 1e10  # K

# Early density (adjusted to match R_b; paper's 10^{17} is inconsistent)
rho_early = 114  # kg/m³

# Aether couplings α ≈ c1 + c3 ≈ 0.1
alpha = 0.1

# Vorticity magnitude ω ≈ 10^{-18} s^{-1} (adjusted to ~10^{-17} for matching high-z)
omega = 1e-17  # s^{-1}

# Scattering opacity σ_scatter ≈ 10^{-18} Mpc^{-1}
sigma_scatter = 1e-18  # Mpc^{-1}
sigma_scatter_si = sigma_scatter / 3.08568e22  # Convert to m^{-1}

# Print key parameters
print(f"Boson mass m: {m:.2e} kg")
print(f"vev v: {v:.2e} m^{-3/2}")
print(f"a: {a:.2e} m⁵ kg⁻¹ s⁻²")
print(f"b: {b:.2e} m³ kg⁻¹")
print(f"gamma_prime: {gamma_prime:.2e} m⁹ kg⁻² s⁻²")
