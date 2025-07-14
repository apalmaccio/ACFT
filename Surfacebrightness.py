import numpy as np
import matplotlib.pyplot as plt

# SB ∝ (1+z)^{-4} from energy loss, bandwidth, dilation, dilution
def sb_dimming(z):
    return (1 + z)**(-4)

# z range
z_vals = np.linspace(0, 10, 100)

sb_vals = sb_dimming(z_vals)

# Plot
plt.plot(z_vals, sb_vals)
plt.xlabel('z')
plt.ylabel('Relative Surface Brightness')
plt.title('SB Dimming in AFCT')
plt.yscale('log')
plt.show()

# Example at z=10
print(f"SB at z=10: {sb_dimming(10):.2e}")
