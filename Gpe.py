import numpy as np
from scipy.fft import fft, ifft, fftfreq
import matplotlib.pyplot as plt
from parameters import m, rho, a_s  # Reuse from previous parameters.py (adjust if needed)

# Grid setup (1D for simplicity; extend to 2D/3D as needed)
N = 512
L = 100  # Box size (arbitrary units, scaled to cosmic scales)
x = np.linspace(-L/2, L/2, N)
dx = x[1] - x[0]
k = 2 * np.pi * fftfreq(N, dx)

# Potential: Harmonic trap + self-interaction (from paper's V(Φ))
def potential(x, psi):
    return 0.5 * x**2  # Harmonic for trap; add cosmic terms if needed
    # Interaction strength g ≈ 4π ħ² |a_s| / m (from scattering length)
    g = 4 * np.pi * (6.626e-34 / (2 * np.pi))**2 * abs(a_s) / m  # ħ = hbar
    return potential(x, psi) + g * np.abs(psi)**2

# Split-step Fourier method for GPE evolution
def evolve_gpe(psi, dt, steps):
    densities = []
    for _ in range(steps):
        # Kinetic step (Fourier space)
        psi_k = fft(psi)
        psi_k *= np.exp(-1j * (k**2 / (2 * m)) * (dt / 2))  # Half-step
        psi = ifft(psi_k)
        
        # Potential step (real space)
        V = 0.5 * x**2  # Example harmonic
        g = 4 * np.pi * (1.0545718e-34)**2 * abs(a_s) / m  # hbar in SI
        psi *= np.exp(-1j * (V + g * np.abs(psi)**2) * dt)
        
        # Kinetic half-step again
        psi_k = fft(psi)
        psi_k *= np.exp(-1j * (k**2 / (2 * m)) * (dt / 2))
        psi = ifft(psi_k)
        
        # Normalize to particle number ~ rho * volume
        norm = np.sqrt(np.trapz(np.abs(psi)**2, x))
        psi /= norm
        
        densities.append(np.abs(psi)**2)
    
    return x, densities

# Initial Gaussian wavefunction (ground state approximation)
psi_init = np.exp(-x**2 / 2) * np.pi**(-0.25)  # Normalized

# Evolve
dt = 0.01
steps = 100
x, densities = evolve_gpe(psi_init, dt, steps)

# Plot final density (should peak at ~rho)
plt.plot(x, densities[-1])
plt.xlabel('Position')
plt.ylabel('Density |ψ|²')
plt.title('BEC Density Evolution in AFCT')
plt.savefig('bec_density.png')
plt.show()

# Check average density ~10^{-26} kg/m³ (scaled)
avg_rho = np.mean(densities[-1]) / (L / np.sqrt(m))  # Rough scaling
print(f"Average density: {avg_rho:.2e} kg/m³")
