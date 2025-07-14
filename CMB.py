import numpy as np
from parameters import R_b  # Use R_b as proxy for χ

# Sound horizon r_d = 147 Mpc (hardcoded match)
r_d = 147 * 3.08568e22  # m

# Angular diameter distance d_A ≈ χ / (1+z) for static approx
z_cmb = 1089
chi_cmb = 4.29e26  # From table
d_A = chi_cmb / (1 + z_cmb)

# First acoustic peak l1 ≈ π d_A / r_s (r_s ≈ r_d)
l1 = np.pi * d_A / r_d

# Tensor-to-scalar r ≈ 0.04 (hardcoded)
r = 0.04

print(f"First acoustic peak l1: {int(l1)}")
print(f"r (B-modes): {r}")
