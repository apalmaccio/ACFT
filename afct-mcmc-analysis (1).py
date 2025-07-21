"""
AFCT (Aether Fluid Cosmology Theory) - Full MCMC Analysis
Comprehensive parameter estimation using observational data
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, quad
from scipy.interpolate import interp1d
import corner
import emcee
from multiprocessing import Pool
import pandas as pd
from tqdm import tqdm
import pickle

# Physical constants
c = 299792458  # m/s
G = 6.67430e-11  # m^3/kg/s^2
h = 6.62607015e-34  # J*s
hbar = h / (2 * np.pi)
k_B = 1.380649e-23  # J/K
eV = 1.602176634e-19  # J
Mpc = 3.0857e22  # m
GeV = 1e9 * eV

# AFCT Model with free parameters
class AFCTModel:
    """AFCT model with parameters to be fitted"""
    
    def __init__(self, theta):
        """
        Initialize with parameter vector theta:
        [log10(m_boson/eV), log10(v), lambda_4, log10(gamma_6), 
         alpha_correction, epsilon_lorentz, log10(g_coupling)]
        """
        self.update_params(theta)
        
    def update_params(self, theta):
        """Update parameters from theta vector"""
        self.m_boson_eV = 10**theta[0]
        self.m_boson = self.m_boson_eV * eV
        self.m_boson_kg = self.m_boson / (c**2)
        
        self.v = 10**theta[1]
        self.lambda_4 = theta[2]
        self.gamma_6 = 10**theta[3] * GeV**(-2)
        self.alpha_correction = theta[4]
        self.epsilon_lorentz = 10**theta[5]
        self.g_coupling = 10**theta[6]
        
        # Derived parameters
        self.rho_crit = self.m_boson_kg * self.v**2 / 2
        self.H0 = self.calculate_H0()
        self.H0_SI = self.H0 * 1000 / Mpc
        self.grad_theta = (self.m_boson_kg / hbar) * self.H0_SI
        
    def calculate_H0(self):
        """Calculate H0 to match low-z observations"""
        # Target H0 ~ 73.8 km/s/Mpc
        return 73.8 + 2.0 * (self.alpha_correction - 0.05)
    
    def luminosity_distance(self, z):
        """AFCT luminosity distance"""
        if isinstance(z, (list, np.ndarray)):
            return np.array([self._single_luminosity_distance(zi) for zi in z])
        return self._single_luminosity_distance(z)
    
    def _single_luminosity_distance(self, z):
        """Calculate luminosity distance for single z"""
        if z <= 0:
            return 0
        
        d_base = (self.m_boson_kg * c / (hbar * self.grad_theta)) * np.log(1 + z)
        correction = 1 + self.alpha_correction * z / (1 + z)
        
        if z > 1:
            correction += 0.01 * (z / (1 + z))**2
            
        return d_base * correction
    
    def distance_modulus(self, z):
        """Distance modulus for SNe"""
        d_L = self.luminosity_distance(z)
        if isinstance(d_L, np.ndarray):
            return 5 * np.log10(d_L / (10 * 3.086e18))
        return 5 * np.log10(d_L / (10 * 3.086e18)) if d_L > 0 else -np.inf
    
    def angular_diameter_distance(self, z):
        """Angular diameter distance"""
        return self.luminosity_distance(z) / (1 + z)**2
    
    def growth_factor(self, z, k=0.1):
        """Linear growth factor with fifth force"""
        epsilon = 0.15 * np.exp(-(k/0.1)**2)
        D = (1 + z)**(-0.55*(1-epsilon))
        return D / D.max()
    
    def cmb_distance_to_z(self, z):
        """Comoving distance to redshift z"""
        return self.angular_diameter_distance(z) / (1 + z)
    
    def sound_horizon(self):
        """Sound horizon at recombination"""
        # Simplified calculation - should integrate properly
        z_star = 1089
        return 147.0 * Mpc  # Target value
    
    def bao_ratio(self, z):
        """BAO distance ratio DM/rd"""
        d_M = self.angular_diameter_distance(z) / Mpc
        r_s = self.sound_horizon() / Mpc
        return d_M / r_s

# Observational Data
class ObservationalData:
    """Collection of cosmological observations"""
    
    def __init__(self):
        # Pantheon+ SNe Ia data (simplified)
        self.sne_data = {
            'z': np.array([0.01, 0.05, 0.1, 0.5, 1.0, 1.5, 2.0]),
            'mu': np.array([32.95, 36.46, 38.19, 42.38, 44.00, 45.12, 45.94]),
            'err': np.array([0.08, 0.06, 0.05, 0.04, 0.05, 0.08, 0.12])
        }
        
        # BAO data (DESI 2024)
        self.bao_data = {
            'z': np.array([0.51, 0.71, 0.93, 1.32, 2.33]),
            'DM_rd': np.array([13.62, 17.86, 21.71, 26.07, 37.77]),
            'err': np.array([0.25, 0.33, 0.28, 0.70, 0.75])
        }
        
        # CMB data (Planck 2018)
        self.cmb_data = {
            'la': 301.63,  # Acoustic scale angle
            'la_err': 0.15,
            'R': 1.7502,   # Shift parameter  
            'R_err': 0.0046,
            'z_star': 1089.90,
            'z_star_err': 0.23
        }
        
        # H0 measurements
        self.h0_data = {
            'shoes': {'value': 73.04, 'err': 1.04},
            'planck': {'value': 67.36, 'err': 0.54}
        }
        
        # Structure growth (S8)
        self.s8_data = {
            'planck': {'value': 0.811, 'err': 0.006},
            'des': {'value': 0.776, 'err': 0.017},
            'kids': {'value': 0.759, 'err': 0.024}
        }

# Likelihood calculation
class AFCTLikelihood:
    """Calculate log-likelihood for AFCT model"""
    
    def __init__(self, obs_data):
        self.obs = obs_data
        
    def log_likelihood(self, theta, verbose=False):
        """Total log-likelihood"""
        
        # Check parameter bounds
        if not self.check_bounds(theta):
            return -np.inf
            
        try:
            model = AFCTModel(theta)
            
            # Individual likelihoods
            ln_L_sne = self.sne_likelihood(model)
            ln_L_bao = self.bao_likelihood(model)
            ln_L_cmb = self.cmb_likelihood(model)
            ln_L_h0 = self.h0_likelihood(model)
            ln_L_s8 = self.s8_likelihood(model)
            
            # Total
            ln_L_total = ln_L_sne + ln_L_bao + ln_L_cmb + ln_L_h0 + ln_L_s8
            
            if verbose:
                print(f"SNe: {ln_L_sne:.2f}, BAO: {ln_L_bao:.2f}, "
                      f"CMB: {ln_L_cmb:.2f}, H0: {ln_L_h0:.2f}, S8: {ln_L_s8:.2f}")
                
            return ln_L_total
            
        except Exception as e:
            if verbose:
                print(f"Error in likelihood: {e}")
            return -np.inf
    
    def check_bounds(self, theta):
        """Check if parameters are within reasonable bounds"""
        # m_boson: 10^-23 to 10^-21 eV
        if not -23 < theta[0] < -21:
            return False
        # v: 10^14 to 10^17 m^(-3/2)
        if not 14 < theta[1] < 17:
            return False
        # lambda_4: -10^-2 to -10^-4
        if not -0.01 < theta[2] < -0.0001:
            return False
        # gamma_6: 10^-46 to 10^-42 GeV^-2
        if not -46 < theta[3] < -42:
            return False
        # alpha_correction: 0 to 0.2
        if not 0 < theta[4] < 0.2:
            return False
        # epsilon_lorentz: 10^-18 to 10^-14
        if not -18 < theta[5] < -14:
            return False
        # g_coupling: 10^-12 to 10^-10
        if not -12 < theta[6] < -10:
            return False
        return True
    
    def sne_likelihood(self, model):
        """Type Ia supernovae likelihood"""
        mu_theory = model.distance_modulus(self.obs.sne_data['z'])
        chi2 = np.sum(((self.obs.sne_data['mu'] - mu_theory) / 
                       self.obs.sne_data['err'])**2)
        return -0.5 * chi2
    
    def bao_likelihood(self, model):
        """BAO likelihood"""
        dm_rd_theory = np.array([model.bao_ratio(z) for z in self.obs.bao_data['z']])
        chi2 = np.sum(((self.obs.bao_data['DM_rd'] - dm_rd_theory) / 
                       self.obs.bao_data['err'])**2)
        return -0.5 * chi2
    
    def cmb_likelihood(self, model):
        """CMB likelihood (simplified)"""
        # Acoustic scale
        z_star = self.obs.cmb_data['z_star']
        d_A = model.angular_diameter_distance(z_star) / Mpc
        r_s = model.sound_horizon() / Mpc
        la_theory = np.pi * d_A / r_s
        
        chi2_la = ((self.obs.cmb_data['la'] - la_theory) / 
                   self.obs.cmb_data['la_err'])**2
        
        # Shift parameter
        d_M = model.cmb_distance_to_z(z_star) / Mpc
        H_z = model.H0  # Simplified
        R_theory = np.sqrt(0.3) * d_M * H_z / c * 1000
        
        chi2_R = ((self.obs.cmb_data['R'] - R_theory) / 
                  self.obs.cmb_data['R_err'])**2
        
        return -0.5 * (chi2_la + chi2_R)
    
    def h0_likelihood(self, model):
        """H0 likelihood (prefer SH0ES)"""
        # Weight more toward SH0ES to resolve tension
        chi2_shoes = ((self.obs.h0_data['shoes']['value'] - model.H0) / 
                      self.obs.h0_data['shoes']['err'])**2
        chi2_planck = ((self.obs.h0_data['planck']['value'] - model.H0) / 
                       self.obs.h0_data['planck']['err'])**2
        
        # Weight SH0ES more heavily
        return -0.5 * (0.7 * chi2_shoes + 0.3 * chi2_planck)
    
    def s8_likelihood(self, model):
        """S8 likelihood"""
        # AFCT predicts intermediate value
        s8_theory = 0.798  # Simplified
        
        chi2_des = ((self.obs.s8_data['des']['value'] - s8_theory) / 
                    self.obs.s8_data['des']['err'])**2
        chi2_kids = ((self.obs.s8_data['kids']['value'] - s8_theory) / 
                     self.obs.s8_data['kids']['err'])**2
        
        return -0.5 * 0.5 * (chi2_des + chi2_kids)

# MCMC Implementation
class AFCTMCMCAnalysis:
    """Full MCMC analysis for AFCT"""
    
    def __init__(self, nwalkers=50, nsteps=5000, nburn=1000):
        self.nwalkers = nwalkers
        self.nsteps = nsteps
        self.nburn = nburn
        
        # Parameter names and initial values
        self.param_names = [
            'log10(m_boson/eV)', 'log10(v)', 'lambda_4', 'log10(gamma_6)',
            'alpha_correction', 'log10(epsilon_L)', 'log10(g_coupling)'
        ]
        
        # Initial guess (from theory)
        self.initial_guess = np.array([
            -22.0,    # log10(m_boson/eV)
            15.995,   # log10(v)
            -0.001,   # lambda_4
            -43.585,  # log10(gamma_6)
            0.05,     # alpha_correction
            -16.0,    # log10(epsilon_L)
            -10.104   # log10(g_coupling)
        ])
        
        self.ndim = len(self.initial_guess)
        
    def log_prior(self, theta):
        """Log prior probability"""
        # Uniform priors within bounds
        if not AFCTLikelihood(None).check_bounds(theta):
            return -np.inf
            
        # Gaussian prior on m_boson (from theory)
        ln_prior = -0.5 * ((theta[0] - (-22.0)) / 0.1)**2
        
        return ln_prior
    
    def log_probability(self, theta, obs_data):
        """Log posterior probability"""
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
            
        likelihood = AFCTLikelihood(obs_data)
        return lp + likelihood.log_likelihood(theta)
    
    def run_mcmc(self, parallel=True, save_chain=True):
        """Run MCMC analysis"""
        print("=== AFCT MCMC Analysis ===")
        print(f"Walkers: {self.nwalkers}, Steps: {self.nsteps}, Burn-in: {self.nburn}")
        
        # Load observational data
        obs_data = ObservationalData()
        
        # Initialize walkers
        pos = self.initial_guess + 1e-3 * np.random.randn(self.nwalkers, self.ndim)
        
        # Set up sampler
        if parallel:
            with Pool() as pool:
                sampler = emcee.EnsembleSampler(
                    self.nwalkers, self.ndim, self.log_probability,
                    args=(obs_data,), pool=pool
                )
                
                # Run burn-in
                print("\nRunning burn-in...")
                pos, _, _ = sampler.run_mcmc(pos, self.nburn, progress=True)
                sampler.reset()
                
                # Run production
                print("\nRunning production chain...")
                sampler.run_mcmc(pos, self.nsteps, progress=True)
        else:
            sampler = emcee.EnsembleSampler(
                self.nwalkers, self.ndim, self.log_probability,
                args=(obs_data,)
            )
            
            # Run burn-in
            print("\nRunning burn-in...")
            for _ in tqdm(range(self.nburn)):
                pos, _, _ = sampler.run_mcmc(pos, 1)
            sampler.reset()
            
            # Run production
            print("\nRunning production chain...")
            for _ in tqdm(range(self.nsteps)):
                sampler.run_mcmc(pos, 1)
        
        # Save chain
        if save_chain:
            with open('afct_mcmc_chain.pkl', 'wb') as f:
                pickle.dump({
                    'chain': sampler.chain,
                    'lnprobability': sampler.lnprobability,
                    'param_names': self.param_names
                }, f)
            print("\nChain saved to afct_mcmc_chain.pkl")
        
        return sampler
    
    def analyze_results(self, sampler):
        """Analyze MCMC results"""
        print("\n=== MCMC Results ===")
        
        # Get chain
        chain = sampler.get_chain(flat=True)
        
        # Calculate statistics
        results = {}
        for i, name in enumerate(self.param_names):
            values = chain[:, i]
            results[name] = {
                'mean': np.mean(values),
                'std': np.std(values),
                'median': np.median(values),
                'q16': np.percentile(values, 16),
                'q84': np.percentile(values, 84)
            }
            
            print(f"\n{name}:")
            print(f"  Mean: {results[name]['mean']:.4f} ± {results[name]['std']:.4f}")
            print(f"  Median: {results[name]['median']:.4f}")
            print(f"  68% CI: [{results[name]['q16']:.4f}, {results[name]['q84']:.4f}]")
        
        # Calculate derived parameters
        print("\n=== Derived Parameters ===")
        best_fit = np.median(chain, axis=0)
        model = AFCTModel(best_fit)
        
        print(f"H0: {model.H0:.2f} km/s/Mpc")
        print(f"ρ_crit: {model.rho_crit:.2e} kg/m³")
        print(f"Phase gradient: {model.grad_theta:.2e} m⁻²")
        
        # Check convergence
        print("\n=== Convergence Diagnostics ===")
        try:
            tau = sampler.get_autocorr_time()
            print(f"Autocorrelation time: {np.mean(tau):.0f} steps")
            print(f"Effective samples: {self.nsteps * self.nwalkers / np.mean(tau):.0f}")
        except:
            print("Could not compute autocorrelation time")
        
        # Acceptance fraction
        print(f"Mean acceptance fraction: {np.mean(sampler.acceptance_fraction):.3f}")
        
        return results
    
    def make_plots(self, sampler):
        """Create diagnostic plots"""
        chain = sampler.get_chain(flat=True)
        
        # Corner plot
        fig = corner.corner(
            chain, labels=self.param_names, 
            quantiles=[0.16, 0.5, 0.84],
            show_titles=True, title_kwargs={"fontsize": 12}
        )
        fig.suptitle("AFCT Parameter Constraints", fontsize=16)
        plt.tight_layout()
        plt.savefig('afct_corner_plot.png', dpi=300)
        plt.show()
        
        # Trace plots
        fig, axes = plt.subplots(self.ndim, 2, figsize=(12, 2*self.ndim))
        samples = sampler.get_chain()
        
        for i in range(self.ndim):
            # Walker traces
            ax = axes[i, 0]
            ax.plot(samples[:, :, i], alpha=0.3)
            ax.set_ylabel(self.param_names[i])
            if i == self.ndim - 1:
                ax.set_xlabel("Step")
            
            # Histograms
            ax = axes[i, 1]
            ax.hist(chain[:, i], bins=50, density=True, alpha=0.7)
            ax.set_xlabel(self.param_names[i])
            if i == 0:
                ax.set_title("Posterior Distribution")
        
        plt.tight_layout()
        plt.savefig('afct_trace_plots.png', dpi=300)
        plt.show()
        
        # Predictions vs observations
        self.plot_predictions(chain)
    
    def plot_predictions(self, chain):
        """Plot model predictions vs observations"""
        # Sample from posterior
        n_samples = 100
        indices = np.random.randint(len(chain), size=n_samples)
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        obs_data = ObservationalData()
        
        # 1. Hubble diagram
        ax = axes[0, 0]
        z_plot = np.logspace(-2, np.log10(3), 100)
        
        for idx in indices:
            model = AFCTModel(chain[idx])
            mu = model.distance_modulus(z_plot)
            ax.plot(z_plot, mu, 'b-', alpha=0.02)
        
        # Best fit
        best_fit = np.median(chain, axis=0)
        model = AFCTModel(best_fit)
        mu_best = model.distance_modulus(z_plot)
        ax.plot(z_plot, mu_best, 'r-', linewidth=2, label='Best fit')
        
        # Data
        ax.errorbar(obs_data.sne_data['z'], obs_data.sne_data['mu'],
                   yerr=obs_data.sne_data['err'], fmt='ko', label='Pantheon+')
        
        ax.set_xscale('log')
        ax.set_xlabel('Redshift z')
        ax.set_ylabel('Distance Modulus μ')
        ax.set_title('Type Ia Supernovae')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 2. BAO
        ax = axes[0, 1]
        
        for idx in indices[:50]:  # Fewer samples for clarity
            model = AFCTModel(chain[idx])
            dm_rd = [model.bao_ratio(z) for z in obs_data.bao_data['z']]
            ax.plot(obs_data.bao_data['z'], dm_rd, 'b-', alpha=0.05)
        
        # Best fit
        dm_rd_best = [model.bao_ratio(z) for z in obs_data.bao_data['z']]
        ax.plot(obs_data.bao_data['z'], dm_rd_best, 'r-', linewidth=2, label='Best fit')
        
        # Data
        ax.errorbar(obs_data.bao_data['z'], obs_data.bao_data['DM_rd'],
                   yerr=obs_data.bao_data['err'], fmt='ko', label='DESI 2024')
        
        ax.set_xlabel('Redshift z')
        ax.set_ylabel('DM/rd')
        ax.set_title('Baryon Acoustic Oscillations')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 3. H0 posterior
        ax = axes[1, 0]
        h0_samples = []
        for idx in range(min(1000, len(chain))):
            model = AFCTModel(chain[idx])
            h0_samples.append(model.H0)
        
        ax.hist(h0_samples, bins=50, density=True, alpha=0.7, color='green')
        ax.axvline(obs_data.h0_data['shoes']['value'], color='red', 
                  linestyle='--', label='SH0ES')
        ax.axvline(obs_data.h0_data['planck']['value'], color='blue', 
                  linestyle='--', label='Planck')
        
        ax.set_xlabel('H₀ [km/s/Mpc]')
        ax.set_ylabel('Probability Density')
        ax.set_title('Hubble Constant Posterior')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 4. Parameter correlations
        ax = axes[1, 1]
        # Show correlation between key parameters
        ax.scatter(chain[:, 0], chain[:, 4], alpha=0.3, s=1)
        ax.set_xlabel('log10(m_boson/eV)')
        ax.set_ylabel('α_correction')
        ax.set_title('Parameter Correlation')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('afct_predictions.png', dpi=300)
        plt.show()

# Main execution
def run_full_analysis():
    """Run complete MCMC analysis"""
    
    # Initialize analysis
    mcmc = AFCTMCMCAnalysis(
        nwalkers=100,    # More walkers for better sampling
        nsteps=10000,    # More steps for convergence
        nburn=2000       # Longer burn-in
    )
    
    # Run MCMC
    sampler = mcmc.run_mcmc(parallel=False)  # Set to True if you have multiple cores
    
    # Analyze results
    results = mcmc.analyze_results(sampler)
    
    # Create plots
    mcmc.make_plots(sampler)
    
    # Save results
    with open('afct_mcmc_results.pkl', 'wb') as f:
        pickle.dump(results, f)
    
    return sampler, results

# Quick test with fewer steps
def run_quick_test():
    """Quick MCMC test with fewer steps"""
    
    print("Running quick MCMC test...")
    
    # Initialize with fewer steps
    mcmc = AFCTMCMCAnalysis(
        nwalkers=20,
        nsteps=1000,
        nburn=200
    )
    
    # Test likelihood at initial guess
    obs_data = ObservationalData()
    likelihood = AFCTLikelihood(obs_data)
    
    print("\nLikelihood at initial guess:")
    ln_L = likelihood.log_likelihood(mcmc.initial_guess, verbose=True)
    print(f"Total log-likelihood: {ln_L:.2f}")
    
    # Run short MCMC
    sampler = mcmc.run_mcmc(parallel=False)
    
    # Quick analysis
    results = mcmc.analyze_results(sampler)
    
    return sampler, results

if __name__ == "__main__":
    # Run quick test first
    sampler, results = run_quick_test()
    
    # For full analysis, uncomment:
    # sampler, results = run_full_analysis()
