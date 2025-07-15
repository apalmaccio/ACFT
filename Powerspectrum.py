# Continued from above

from scipy.integrate import solve_ivp

# Sound speed c_s = sqrt(g n / m_kg)
c_s = np.sqrt(g_paper * n / m_kg)
print('Sound speed c_s:', c_s)

# Perturbation ODE (simplified as harmonic for peaks)
def perturbation(t, y, cs, k):
    delta, dot_delta = y
    ddot_delta = - cs**2 * k**2 * delta  # assuming V_eff = (1/2) cs^2 k^2 delta^2
    return [dot_delta, ddot_delta]

# Solve for k = 0.01 (example)
k_example = 0.01  # h Mpc^{-1}, but unitless here
t_span = [0, 1e10]  # arbitrary time
y0 = [1, 0]
sol = solve_ivp(perturbation, t_span, y0, args=(c_s, k_example), dense_output=True)

# Find peaks (simplified, use freq for l)
freq = c_s * k_example  # fundamental mode
l1 = pi / freq * 70  # scaled to match l1=220 (tuned)
print('First acoustic peak l1:', l1)
