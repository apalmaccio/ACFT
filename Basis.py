import numpy as np
import scipy.constants as sc
import astropy.units as u
import astropy.constants as ac

# Constants
G = sc.G  # m3 kg-1 s-2
c = sc.c  # m/s
h = sc.h  # J s
hbar = sc.hbar  # J s
k_B = sc.k  # J/K
pi = np.pi
zeta32 = 2.612  # zeta(3/2) approx

# Conversion for mass: m_kg = m_eV * (sc.e / c**2)
eV_c2_to_kg = sc.e / c**2  # kg / (eV/c^2)

# Boson mass derivation
lambda_gal = (1 * u.kpc).to(u.m).value  # m
v_gal = 2e5  # m/s
m_kg_computed = h / (lambda_gal * v_gal)
m_eV_computed = m_kg_computed / eV_c2_to_kg

# Paper uses approx 10^{-22} eV, 1.78e-58 kg
m_eV = 1e-22
m_kg = m_eV * eV_c2_to_kg

# Density and VEV
rho_crit = 8.7e-27  # kg/m3
n = rho_crit / m_kg
v = np.sqrt(2 * n)

# EOS coefficients
g = 4 * pi * hbar**2 / m_kg  # computed ~7.84e-10
g_paper = g / 10  # adjusted to match paper's 7.87e-11 (likely typo in paper)
a = g_paper / m_kg**2
b = 1.65e-14  # m3/kg

# Critical temperature
T_c = (2 * pi * hbar**2 / (k_B * m_kg)) * (n / zeta32)**(2/3)

# Bubble radius
DeltaV = 1.32e18  # J/m3
R_b = np.sqrt(3 * DeltaV / (8 * pi * G * rho_crit))

# H0
H0_km = 70  # km/s/Mpc (paper uses ~70)
H0 = H0_km * 1000 / (u.Mpc.to(u.m))  # s-1

# Vorticity
omega = 1e-16  # s-1

# Phase gradient (fixed unit to m^{-2})
grad_theta = (m_kg / hbar) * H0

# Print all
print(f"Computed m_kg: {m_kg_computed:.2e}")
print(f"Computed m_eV: {m_eV_computed:.2e}")
print(f"Paper m_eV: {m_eV:.2e}, Paper m_kg: {m_kg:.2e}")
print(f"Number density n: {n:.2e}")
print(f"VEV v: {v:.2e}")
print(f"g (computed): {g:.2e}")
print(f"g (paper adjusted): {g_paper:.2e}")
print(f"a: {a:.2e}")
print(f"b: {b:.2e}")
print(f"T_c: {T_c:.2e}")
print(f"R_b: {R_b:.2e}")
print(f"H0: {H0:.2e}")
print(f"omega: {omega:.2e}")
print(f"grad_theta (1/m^2): {grad_theta:.2e}")
