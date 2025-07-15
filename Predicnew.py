# Continued from above

# Growth ODE
def growth_ode(z, y):
    delta, ddelta_dz = y
    omz = (1 + z)**2  # simplified Ω_m(z)
    d2delta_dz2 = -1/(1 + z) * ddelta_dz + 1.5 * omz * delta
    return [ddelta_dz, d2delta_dz2]

z_span = [0, 10]
y0 = [1, 0]
sol_growth = solve_ivp(growth_ode, z_span, y0, dense_output=True)
print('Growth at z=0:', sol_growth.y[0][-1])
