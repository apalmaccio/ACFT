

# Redshift formula
def z_func(d):
    prefactor = hbar / (m_kg * c)
    integral = grad_theta * d  # since constant grad_theta
    return np.exp(prefactor * integral) - 1

# Example: d for CMB z=1089
z_cmb = 1089
d_cmb = np.log(1 + z_cmb) / (hbar / (m_kg * c) * grad_theta)  # inverse
print('Computed chi(z=1089):', d_cmb)

# Plot z vs d for small z (mimic Hubble)
d_array = np.linspace(0, 1e26, 100)  # m
z_array = z_func(d_array)
plt.plot(d_array, z_array)
plt.xlabel('Distance (m)')
plt.ylabel('Redshift z')
plt.title('Redshift vs Distance')
plt.show()

# Effective Hubble check for small z
d_small = 3.086e22  # 1 Mpc
z_small = z_func(d_small)
H0_eff = z_small * c / d_small
print('Effective H0 for small z:', H0_eff)
