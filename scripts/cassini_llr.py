import numpy as np

G = 6.674e-11  # m^3 kg^-1 s^-2
g_phi = 4.1e-19  # m^2 kg^-1 s^2
G_eff = G + g_phi**2 / (4 * np.pi)
M_sun = 1.989e30  # kg
c = 3e8  # m/s
R_sun = 6.96e8  # m
AU = 1.496e11  # m
delta_t_GR = (2 * G * M_sun / c**3) * np.log(4 * AU**2 / R_sun**2)
delta_t_ACFT = delta_t_GR * (G_eff / G)
shapiro_deviation = (delta_t_ACFT - delta_t_GR) * 1e6
print(f"Shapiro delay deviation = {shapiro_deviation} μs")

M_earth = 5.972e24  # kg
r_moon = 3.844e8  # m
a_extra = (g_phi**2 * M_earth) / (4 * np.pi * r_moon**2)
print(f"LLR extra acceleration = {a_extra} m/s^2")
