import numpy as np

# Prediction 1: SB at high z with 5-10% excess due to vorticity lensing
def sb_excess(z):
    base = (1 + z)**(-4)
    excess = 1 + 0.075 * (z > 10)  # 7.5% average excess for z>10
    return base * excess

print(f"SB at z=12 with excess: {sb_excess(12):.2e}")

# Prediction 6: GW burst from phase transition
f_gw = 1e-3  # Hz
Omega_gw_h2 = 1e-10
print(f"GW burst: f={f_gw} Hz, Ω_GW h²={Omega_gw_h2:.2e}")
