import numpy as np
from scipy.fft import fftn, ifftn
from scipy.interpolate import RegularGridInterpolator

# Parameters
L = 50.0  # Mpc/h
N = 64  # 64^3 particles
M_tot = 1e14  # M_sun/h
m_p = M_tot / N**3
G = 4.302e-3  # (pc/M_sun) (km/s)^2
g_phi = 4.1e-19 * 3.086e22  # Mpc kg^-1 s^2
dx = L / N
k2 = np.fft.fftfreq(N, d=dx) * 2 * np.pi
kx, ky, kz = np.meshgrid(k2, k2, k2)
k_sq = kx**2 + ky**2 + kz**2
k_sq[0, 0, 0] = 1e-10

# Initial conditions
np.random.seed(42)
pos = np.random.uniform(0, L, (N**3, 3))
vel = np.zeros_like(pos)
rho = np.ones((N, N, N)) * M_tot / L**3

def solve_poisson(rho):
    rho_k = fftn(rho)
    phi_k = -4 * np.pi * G * rho_k / k_sq * (1 + g_phi**2 / (4 * np.pi * G))
    phi = np.real(ifftn(phi_k))
    return phi

def compute_force(phi):
    grad_x = np.gradient(phi, dx, axis=0)
    grad_y = np.gradient(phi, dx, axis=1)
    grad_z = np.gradient(phi, dx, axis=2)
    interp = RegularGridInterpolator((np.arange(N)*dx, np.arange(N)*dx, np.arange(N)*dx),
                                    np.stack([grad_x, grad_y, grad_z], axis=-1))
    forces = interp(pos)
    return -forces

dt = 0.01  # Gyr
n_steps = 100
for step in range(n_steps):
    rho = np.zeros((N, N, N))
    for p in pos:
        i, j, k = (p / dx).astype(int) % N
        rho[i, j, k] += m_p / dx**3
    rho -= np.mean(rho)
    phi = solve_poisson(rho)
    acc = compute_force(phi)
    vel += acc * dt
    pos = (pos + vel * dt) % L

from scipy.ndimage import label
rho_smooth = np.abs(rho) > 10 * np.std(rho)
labels, n_halos = label(rho_smooth)
masses = []
for i in range(1, n_halos + 1):
    mass = np.sum(rho[labels == i]) * dx**3
    if mass > 1e12:
        masses.append(mass)

bins = np.logspace(12, 14, 10)
hmf, _ = np.histogram(masses, bins=bins, density=True)
hmf /= np.diff(bins)
np.savetxt('../data/synthetic_sdss_hmf.csv', np.column_stack((bins[:-1], hmf)), delimiter=',')
print(f"HMF: {hmf}")
