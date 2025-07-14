import numpy as np
import matplotlib.pyplot as plt
from parameters import alpha  # Aether coupling ~0.1

# Time dilation factor with chromaticity
def time_dilation(z, wavelength):
    base = 1 + z
    scatter = 0.075 * (wavelength / 500e-9)  # 5-10% UV vs IR (wavelength in m)
    return base * (1 + scatter * np.random.uniform(-0.05, 0.05))  # Random scatter

# z range
z_vals = np.linspace(0, 10, 100)

# UV (300nm) vs IR (1000nm)
td_uv = [time_dilation(z, 300e-9) for z in z_vals]
td_ir = [time_dilation(z, 1000e-9) for z in z_vals]

# Plot
plt.plot(z_vals, td_uv, label='UV')
plt.plot(z_vals, td_ir, label='IR')
plt.xlabel('z')
plt.ylabel('Time Dilation Factor')
plt.title('Chromatic Time Dilation in AFCT')
plt.legend()
plt.show()

# Scatter at z=10
scatter_pct = abs(td_uv[-1] - td_ir[-1]) / np.mean([td_uv[-1], td_ir[-1]]) * 100
print(f"Chromatic scatter at z=10: {scatter_pct:.1f}%")
