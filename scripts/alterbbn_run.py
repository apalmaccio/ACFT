import numpy as np

# Mock AlterBBN output (replace with actual AlterBBN call)
def run_bbn(params):
    T_burst = params['T_burst']
    mu = params['mu']
    T = params['T']
    # Placeholder: Actual AlterBBN call would go here
    Yp = 0.2458
    Li_H = 1.39e-10
    return Yp, Li_H

params = {'T_burst': 1.2e12, 'mu': 1.1e-8, 'T': 1e6}
Yp, Li_H = run_bbn(params)
print(f"Y_p = {Yp}, Li/H = {Li_H}")
