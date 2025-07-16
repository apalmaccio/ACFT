import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, quad
from scipy.interpolate import interp1d
import pandas as pd

# Physical constants
c = 299792458  # m/s
G = 6.67430e-11  # m^3/kg/s^2
h = 6.62607015e-34  # J*s
hbar = h / (2 * np.pi)
k_B = 1.380649e-23  # J/K
eV = 1.602176634e-19  # J
Mpc = 3.0857e22  # m
yr = 365.25 * 24 * 3600  # s

# AFCT Parameters from the paper
m_boson = 1e-22 * eV  # boson mass in J
m_boson_kg = 1.78e-58  # boson mass in kg
rho_crit = 8.7e-27  # kg/m^3
v = 9.89e15  # m^(-3/2) - vacuum expectation value
T_c = 2.00e34  # K - critical temperature
R_b = 5.21e26  # m - bubble radius
omega = 1e-16  # s^(-1) - vorticity
H0 = 70  # km/s/Mpc
H0_SI = H0 * 1000 / Mpc  # s^(-1)
grad_theta = (m_boson_kg / hbar) * H0_SI  # m^(-2)
alpha_correction = 0.05  # quantum correction term

# CMB parameters
T_CMB_0 = 2.725  # K
z_CMB = 1089

# Simulation functions

def redshift_AFCT(distance, grad_theta=grad_theta):
    """Calculate redshift in AFCT from phase gradient"""
    z = np.exp((hbar / (m_boson_kg * c)) * grad_theta * distance) - 1
    return z

def luminosity_distance_AFCT(z, alpha=alpha_correction):
    """Calculate luminosity distance in AFCT"""
    # Including quantum correction term
    integral_part = (m_boson_kg * c / (hbar * grad_theta)) * np.log(1 + z)
    correction = alpha * z / (1 + z)
    d_L = integral_part * (1 + correction)
    return d_L

def luminosity_distance_LCDM(z, Omega_m=0.3, Omega_L=0.7):
    """Standard ΛCDM luminosity distance for comparison"""
    def integrand(z_prime):
        return 1 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)
    
    integral, _ = quad(integrand, 0, z)
    d_L = (c / H0_SI) * (1 + z) * integral
    return d_L

def angular_diameter_distance_AFCT(z):
    """Angular diameter distance in AFCT"""
    d_L = luminosity_distance_AFCT(z)
    d_A = d_L / (1 + z)**2
    return d_A

def CMB_temperature_AFCT(z, gamma_int_ratio=0.01):
    """CMB temperature evolution in AFCT"""
    delta_T = (gamma_int_ratio) * 1e-4 * np.log(1 + z)
    T = T_CMB_0 * (1 + z) * (1 + delta_T)
    return T

def matter_power_spectrum_enhancement(k, k_screen=0.1):
    """Fifth-force enhancement of matter power spectrum"""
    # Enhancement below screening scale
    enhancement = 1 + 0.15 * np.exp(-(k/k_screen)**2)
    return enhancement

def surface_brightness_dimming_AFCT(z):
    """Surface brightness dimming in AFCT"""
    # Standard (1+z)^4 plus additional BEC effects
    standard_dimming = (1 + z)**4
    bec_excess = 1 + 0.05 * z / (1 + z/10)  # 5-10% excess at high z
    return standard_dimming * bec_excess

def acoustic_peak_positions():
    """CMB acoustic peak positions in AFCT"""
    # From the paper
    peaks = {
        'l1': 220,
        'l2': 540,
        'l3': 840,
        'l4': 1150,
        'l5': 1450
    }
    return peaks

# Run simulations and create plots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('AFCT Cosmology Simulations', fontsize=16)

# 1. Luminosity Distance Comparison
z_range = np.logspace(-2, 0.5, 100)
d_L_AFCT = np.array([luminosity_distance_AFCT(z) for z in z_range])
d_L_LCDM = np.array([luminosity_distance_LCDM(z) for z in z_range])

ax1 = axes[0, 0]
ax1.plot(z_range, d_L_AFCT/1e9/Mpc, 'b-', label='AFCT', linewidth=2)
ax1.plot(z_range, d_L_LCDM/1e9/Mpc, 'r--', label='ΛCDM', linewidth=2)
ax1.set_xlabel('Redshift z')
ax1.set_ylabel('Luminosity Distance (Gpc)')
ax1.set_title('Luminosity Distance Comparison')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xscale('log')

# 2. Relative Deviation in Luminosity Distance
ax2 = axes[0, 1]
deviation = (d_L_AFCT - d_L_LCDM) / d_L_LCDM * 100
ax2.plot(z_range, deviation, 'g-', linewidth=2)
ax2.set_xlabel('Redshift z')
ax2.set_ylabel('Deviation (%)')
ax2.set_title('AFCT vs ΛCDM Distance Deviation')
ax2.grid(True, alpha=0.3)
ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)

# Add specific points mentioned in the paper
z_points = [0.5, 1.0, 2.0]
for z in z_points:
    dev = (luminosity_distance_AFCT(z) - luminosity_distance_LCDM(z)) / luminosity_distance_LCDM(z) * 100
    ax2.plot(z, dev, 'ro', markersize=8)
    ax2.annotate(f'z={z}: {dev:.1f}%', (z, dev), xytext=(5, 5), 
                 textcoords='offset points', fontsize=8)

# 3. CMB Temperature Evolution
z_cmb_range = np.logspace(0, 3.1, 100)
T_cmb = np.array([CMB_temperature_AFCT(z) for z in z_cmb_range])

ax3 = axes[0, 2]
ax3.plot(z_cmb_range, T_cmb, 'r-', linewidth=2)
ax3.axhline(y=T_CMB_0 * (1 + z_CMB), color='k', linestyle='--', 
            label=f'Standard at z={z_CMB}')
ax3.set_xlabel('Redshift z')
ax3.set_ylabel('CMB Temperature (K)')
ax3.set_title('CMB Temperature Evolution')
ax3.set_xscale('log')
ax3.grid(True, alpha=0.3)
ax3.legend()

# 4. Matter Power Spectrum Enhancement
k_range = np.logspace(-2, 1, 100)  # h/Mpc
P_enhancement = matter_power_spectrum_enhancement(k_range)

ax4 = axes[1, 0]
ax4.plot(k_range, (P_enhancement - 1) * 100, 'purple', linewidth=2)
ax4.set_xlabel('k (h/Mpc)')
ax4.set_ylabel('Enhancement (%)')
ax4.set_title('Matter Power Spectrum Enhancement')
ax4.set_xscale('log')
ax4.grid(True, alpha=0.3)
ax4.axvline(x=0.1, color='k', linestyle='--', alpha=0.5, 
            label='Screening scale')
ax4.legend()

# 5. Surface Brightness Dimming
z_sb_range = np.logspace(-1, 1.2, 100)
sb_standard = (1 + z_sb_range)**4
sb_afct = surface_brightness_dimming_AFCT(z_sb_range)

ax5 = axes[1, 1]
ax5.plot(z_sb_range, sb_afct/sb_standard - 1, 'orange', linewidth=2)
ax5.set_xlabel('Redshift z')
ax5.set_ylabel('Excess Dimming Factor')
ax5.set_title('Surface Brightness Excess Dimming')
ax5.set_xscale('log')
ax5.grid(True, alpha=0.3)
ax5.axhline(y=0, color='k', linestyle='--', alpha=0.5)

# 6. CMB Acoustic Peaks
peaks = acoustic_peak_positions()
l_values = list(peaks.values())
peak_names = list(peaks.keys())

ax6 = axes[1, 2]
ax6.bar(range(len(l_values)), l_values, color='cyan', alpha=0.7)
ax6.set_xticks(range(len(l_values)))
ax6.set_xticklabels(peak_names)
ax6.set_ylabel('Multipole l')
ax6.set_title('CMB Acoustic Peak Positions')
ax6.grid(True, alpha=0.3, axis='y')

# Add theoretical values as horizontal lines
theoretical_peaks = [220, 540, 840, 1150, 1450]
for i, l in enumerate(theoretical_peaks):
    ax6.axhline(y=l, color='red', linestyle='--', alpha=0.5, linewidth=1)

plt.tight_layout()
plt.show()

# Generate summary statistics
print("AFCT Cosmology Simulation Results")
print("=" * 50)
print(f"\nFundamental Parameters:")
print(f"Boson mass: {m_boson/eV:.2e} eV")
print(f"Critical density: {rho_crit:.2e} kg/m³")
print(f"Phase gradient: {grad_theta:.2e} m⁻²")
print(f"Hubble-like parameter H₀: {H0} km/s/Mpc")

print(f"\nKey Predictions at Selected Redshifts:")
for z in [0.1, 0.5, 1.0, 2.0, 5.0]:
    d_afct = luminosity_distance_AFCT(z) / 1e9 / Mpc
    d_lcdm = luminosity_distance_LCDM(z) / 1e9 / Mpc
    deviation = (d_afct - d_lcdm) / d_lcdm * 100
    print(f"\nz = {z}:")
    print(f"  AFCT d_L = {d_afct:.3f} Gpc")
    print(f"  ΛCDM d_L = {d_lcdm:.3f} Gpc")
    print(f"  Deviation = {deviation:.1f}%")

print(f"\nCMB Predictions:")
print(f"Temperature at z={z_CMB}: {CMB_temperature_AFCT(z_CMB):.1f} K")
print(f"Expected: {T_CMB_0 * (1 + z_CMB):.1f} K")
print(f"Deviation: {(CMB_temperature_AFCT(z_CMB)/(T_CMB_0*(1+z_CMB)) - 1)*100:.3f}%")

print(f"\nMatter Power Spectrum:")
print(f"Enhancement at k=0.01 h/Mpc: {(matter_power_spectrum_enhancement(0.01)-1)*100:.1f}%")
print(f"Enhancement at k=0.1 h/Mpc: {(matter_power_spectrum_enhancement(0.1)-1)*100:.1f}%")

print(f"\nTestable Predictions for Upcoming Surveys:")
print(f"1. DESI (2025-2030): {deviation:.1f}% deviation in d_L at z=0.5")
print(f"2. Euclid: Enhanced cosmic shear by ~{(matter_power_spectrum_enhancement(0.05)-1)*100:.0f}%")
print(f"3. JWST: Surface brightness excess of {(surface_brightness_dimming_AFCT(10)/(1+10)**4 - 1)*100:.0f}% at z=10")
print(f"4. LISA: GW birefringence ~1% at f~10⁻⁴ Hz")
print(f"5. Rubin Observatory: Non-achromatic time dilation scatter 5-10%")
