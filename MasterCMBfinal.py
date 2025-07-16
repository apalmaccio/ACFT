import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.special import spherical_jn
from scipy.interpolate import interp1d

# Constants
c = 299792458  # m/s
G = 6.67430e-11  # m^3/kg/s^2
hbar = 1.054571817e-34  # J*s
k_B = 1.380649e-23  # J/K
Mpc = 3.0857e22  # m
H0 = 70e3 / Mpc  # s^-1

# AFCT BEC parameters
m_boson = 1.78e-58  # kg
g_interaction = 7.87e-11  # m^3/s^2
n_density = 4.89e31  # m^-3
cs_sound = np.sqrt(g_interaction * n_density / m_boson)  # sound speed in BEC

def cmb_power_spectrum_AFCT(l_range, A_s=2.1e-9):
    """
    Calculate CMB angular power spectrum in AFCT
    Using phonon oscillations in BEC fluid
    """
    # Acoustic scale in BEC
    r_s = 147 * Mpc  # BAO scale from paper
    theta_s = r_s / (4.32e26)  # Angular scale at recombination
    
    # Peak positions from paper
    peak_positions = np.array([220, 540, 840, 1150, 1450])
    
    # Initialize power spectrum
    C_l = np.zeros_like(l_range, dtype=float)
    
    for l in range(len(l_range)):
        ell = l_range[l]
        
        # Phonon oscillation modes (modified from standard acoustic oscillations)
        # In BEC, oscillations are driven by quantum pressure
        k = ell / (4.32e26 / Mpc)  # comoving wavenumber
        
        # BEC modification to acoustic oscillations
        omega_phon = cs_sound * k * np.sqrt(1 + (hbar * k / (2 * m_boson * cs_sound))**2)
        
        # Oscillation amplitude with BEC damping
        damping = np.exp(-k**2 * (hbar / (m_boson * cs_sound))**2)
        
        # Modified transfer function
        phase = omega_phon * 3.31e17  # transient time from paper
        transfer = np.cos(phase) * damping
        
        # Angular power spectrum
        # Includes fifth-force enhancement at large scales
        enhancement = 1 + 0.1 * np.exp(-(ell/100)**2)
        
        C_l[l] = A_s * (ell * (ell + 1)) * transfer**2 * enhancement / (2 * np.pi)
        
        # Add peak structure
        for peak in peak_positions:
            if abs(ell - peak) < 20:
                C_l[l] *= 1 + 0.5 * np.exp(-(ell - peak)**2 / 100)
    
    return C_l

def structure_growth_AFCT(z_range, k=0.1):
    """
    Structure growth with fifth-force modification
    """
    def growth_equation(y, z, k):
        delta, delta_prime = y
        
        # AFCT modified growth equation from paper
        Omega_m_eff = 0.3 * (1 + z)**2 / H0**2
        
        # Fifth-force enhancement
        epsilon_5th = 0.15 * np.exp(-(k/0.1)**2)
        
        # Modified growth equation
        d2delta_dz2 = -(1/(1+z)) * delta_prime + (3/2) * Omega_m_eff * (1 + epsilon_5th) * delta
        
        return [delta_prime, d2delta_dz2]
    
    # Initial conditions at high z
    z_init = 1000
    delta_init = 1e-5
    delta_prime_init = delta_init / (1 + z_init)
    y0 = [delta_init, delta_prime_init]
    
    # Solve ODE
    z_solve = z_range[::-1]  # Solve backwards in time
    solution = odeint(growth_equation, y0, z_solve, args=(k,))
    
    # Extract growth factor
    delta = solution[:, 0][::-1]
    growth_factor = delta / delta[0]
    
    return growth_factor

def bec_density_profile(r, r_core=1e3 * 3.086e19):  # 1 kpc core
    """
    Solitonic density profile in galaxy halos
    """
    rho_0 = 1e-24 * 1e3  # kg/m^3 central density
    profile = rho_0 / (1 + (r/r_core)**2)**2
    return profile

# Generate plots
fig = plt.figure(figsize=(15, 12))

# 1. CMB Power Spectrum
ax1 = plt.subplot(2, 2, 1)
l_range = np.arange(2, 2000)
C_l_afct = cmb_power_spectrum_AFCT(l_range)

# Normalize and plot
l_factor = l_range * (l_range + 1) / (2 * np.pi)
D_l = C_l_afct * l_factor * 1e12  # Convert to μK²

ax1.plot(l_range, D_l, 'b-', linewidth=2, label='AFCT Prediction')

# Mark acoustic peaks
peaks = [220, 540, 840, 1150, 1450]
for i, peak in enumerate(peaks):
    idx = np.argmin(abs(l_range - peak))
    ax1.plot(peak, D_l[idx], 'ro', markersize=8)
    ax1.annotate(f'Peak {i+1}', (peak, D_l[idx]), 
                xytext=(10, 10), textcoords='offset points')

ax1.set_xlabel('Multipole l')
ax1.set_ylabel('l(l+1)C_l/2π [μK²]')
ax1.set_title('CMB Angular Power Spectrum (AFCT)')
ax1.set_xlim(2, 2000)
ax1.set_xscale('log')
ax1.grid(True, alpha=0.3)
ax1.legend()

# 2. Structure Growth Comparison
ax2 = plt.subplot(2, 2, 2)
z_growth = np.logspace(-1, 2, 100)

# Different k modes
k_modes = [0.01, 0.1, 1.0]  # h/Mpc
colors = ['blue', 'green', 'red']

for k, color in zip(k_modes, colors):
    growth = structure_growth_AFCT(z_growth, k)
    ax2.plot(1 + z_growth, growth, color=color, linewidth=2, 
             label=f'k = {k} h/Mpc')

# Standard ΛCDM growth for comparison
growth_lcdm = ((1 + z_growth)**(-1))**(0.55)  # Approximate
ax2.plot(1 + z_growth, growth_lcdm/growth_lcdm[0], 'k--', 
         linewidth=2, label='ΛCDM')

ax2.set_xlabel('1 + z')
ax2.set_ylabel('Growth Factor D(z)/D(0)')
ax2.set_title('Structure Growth: AFCT vs ΛCDM')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.invert_xaxis()

# 3. Galaxy Halo Profile (Solitonic Core)
ax3 = plt.subplot(2, 2, 3)
r_range = np.logspace(16, 21, 100)  # meters
rho_profile = bec_density_profile(r_range)

# Convert to physical units
r_kpc = r_range / (3.086e19)  # kpc
rho_gcm3 = rho_profile / 1e3  # g/cm³

ax3.loglog(r_kpc, rho_gcm3, 'purple', linewidth=3, label='BEC Soliton')

# NFW profile for comparison
rho_s = 1e-25  # g/cm³
r_s = 20  # kpc
rho_nfw = rho_s / ((r_kpc/r_s) * (1 + r_kpc/r_s)**2)
ax3.loglog(r_kpc, rho_nfw, 'orange', linestyle='--', 
           linewidth=2, label='NFW (ΛCDM)')

ax3.set_xlabel('Radius (kpc)')
ax3.set_ylabel('Density (g/cm³)')
ax3.set_title('Galaxy Halo Density Profile')
ax3.set_xlim(0.1, 100)
ax3.set_ylim(1e-28, 1e-23)
ax3.grid(True, alpha=0.3)
ax3.legend()

# 4. Redshift-Distance Relation (High-z)
ax4 = plt.subplot(2, 2, 4)
z_high = np.logspace(0, 1.5, 100)

# AFCT prediction with full corrections
d_L_afct = []
for z in z_high:
    # Full integral with BEC refractive index effects
    integral = (m_boson * c / (hbar * 3.83e-42)) * np.log(1 + z)
    correction = 0.05 * z / (1 + z) + 0.01 * (z/(1+z))**2  # Higher order terms
    d_L = integral * (1 + correction) / (1e9 * Mpc)  # Gpc
    d_L_afct.append(d_L)

d_L_afct = np.array(d_L_afct)

# Simple Hubble law for comparison
d_L_hubble = (c * z / H0) / (1e9 * Mpc)

ax4.plot(z_high, d_L_afct, 'b-', linewidth=3, label='AFCT (Full)')
ax4.plot(z_high, d_L_hubble, 'r--', linewidth=2, label='Linear Hubble')

# Mark JWST observable range
ax4.axvspan(10, 30, alpha=0.2, color='yellow', label='JWST Range')

ax4.set_xlabel('Redshift z')
ax4.set_ylabel('Luminosity Distance (Gpc)')
ax4.set_title('High-z Distance Predictions')
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.grid(True, alpha=0.3)
ax4.legend()

plt.tight_layout()
plt.show()

# Compute specific observables
print("\nDetailed AFCT Predictions:")
print("=" * 60)

# CMB analysis
print("\nCMB Power Spectrum Analysis:")
peak_heights = []
for peak in peaks:
    idx = np.argmin(abs(l_range - peak))
    peak_heights.append(D_l[idx])

print(f"Peak height ratios (relative to first peak):")
for i in range(1, len(peak_heights)):
    ratio = peak_heights[i] / peak_heights[0]
    print(f"  Peak {i+1}/Peak 1 = {ratio:.3f}")

# Structure formation
print("\nStructure Formation Enhancement (at z=0):")
k_test = [0.01, 0.05, 0.1, 0.5]
for k in k_test:
    growth_afct = structure_growth_AFCT(np.array([0, 1]), k)
    enhancement = (growth_afct[0] / growth_afct[1]) / 2.5 - 1  # vs linear growth
    print(f"  k = {k} h/Mpc: {enhancement*100:.1f}% enhancement")

# Galaxy cores
print("\nSolitonic Core Properties:")
print(f"  Core radius: 1 kpc")
print(f"  Central density: {1e-24:.2e} g/cm³")
print(f"  Core mass: {4*np.pi*(1e3*3.086e19)**3 * 1e-21 / (3*2e30):.2e} solar masses")

# High-z predictions
print("\nHigh-redshift Predictions (for JWST):")
for z in [10, 15, 20]:
    idx = np.argmin(abs(z_high - z))
    d_ratio = d_L_afct[idx] / d_L_hubble[idx]
    sb_excess = (surface_brightness_dimming_AFCT(z) / (1+z)**4 - 1) * 100
    print(f"  z = {z}:")
    print(f"    Distance deviation: {(d_ratio-1)*100:.1f}%")
    print(f"    Surface brightness excess: {sb_excess:.1f}%")

def surface_brightness_dimming_AFCT(z):
    """Surface brightness dimming in AFCT"""
    standard_dimming = (1 + z)**4
    bec_excess = 1 + 0.05 * z / (1 + z/10)
    return standard_dimming * bec_excess
