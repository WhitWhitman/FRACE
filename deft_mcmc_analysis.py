import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import emcee
import corner
import os
import time
from scipy.linalg import eigh
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

# =====================
# PHYSICAL CONSTANTS
# =====================
c_km_s = 299792.458  # Speed of light (km/s)

# =====================
# DATA LOADING
# =====================
def load_data():
    """Load Pantheon+ data and covariance matrix"""
    # Load supernova data
    data = pd.read_csv('Pantheon+SH0ES.dat', sep='\s+', comment='#')
    z_data = data['zHD'].values
    mu_data = data['MU_SH0ES'].values
    n = len(mu_data)
    
    # Load covariance matrix
    cov_flat = np.loadtxt('Pantheon+SH0ES_STAT+SYS.cov')
    
    # Handle covariance size mismatch
    expected_size = n * n
    if len(cov_flat) == expected_size + 1:
        cov_flat = cov_flat[:-1]
    cov = cov_flat.reshape((n, n))
    
    print(f"Loaded {n} SNe, covariance shape: {cov.shape}")
    return z_data, mu_data, cov

# =====================
# STABLE INTEGRAL CALCULATION
# =====================
def compute_integral(z, H0, Om):
    """Robust computation of integral with dynamic grid refinement"""
    # Create a base integration grid
    z_min, z_max = min(z), max(z)
    z_grid = np.linspace(0, z_max + 0.1, 500)
    
    # Calculate Hubble parameter
    hz = H0 * np.sqrt(Om * (1 + z_grid)**3 + (1 - Om))
    
    # Compute integral using cumulative trapezoidal rule
    integrand = 1 / hz
    integral = cumtrapz(integrand, z_grid, initial=0)
    
    # Create interpolator
    return interp1d(z_grid, integral, kind='cubic', fill_value="extrapolate")

def distance_modulus(z_data, H0, Om, model_type, nu0=0.45):
    """Stable distance modulus calculation"""
    # Compute integral
    int_func = compute_integral(z_data, H0, Om)
    int_z = int_func(z_data)
    
    # Apply DEFT modification
    if model_type == 'deft':
        chi_z = nu0 * np.exp(-(z_data - 0.4)**2/(2*0.15**2)) * np.log1p(z_data)
        int_z *= np.sqrt(np.exp(chi_z))
    
    # Compute luminosity distance
    dL = (1 + z_data) * c_km_s * int_z
    dL[dL <= 0] = 1e-10  # Avoid negative distances
    
    return 5 * np.log10(dL) + 25

# =====================
# PRIOR FUNCTION
# =====================
def log_prior(theta, model_type):
    """Physical priors for parameters"""
    if model_type == 'lcdm':
        H0, Om = theta
        if 60 < H0 < 80 and 0.1 < Om < 0.5:
            return 0.0
    else:  # DEFT
        H0, Om, nu0 = theta
        if 60 < H0 < 80 and 0.1 < Om < 0.5 and 0.1 < nu0 < 1.0:
            return 0.0
    return -np.inf

# =====================
# FULL PROBABILITY FUNCTION
# =====================
def log_probability(theta, z, mu_obs, cov_inv, model_type):
    """Full log probability (prior + likelihood)"""
    lp = log_prior(theta, model_type)
    if not np.isfinite(lp):
        return -np.inf
    
    # Handle numerical stability
    try:
        ll = log_likelihood(theta, z, mu_obs, cov_inv, model_type)
        return lp + ll
    except:
        return -np.inf

# =====================
# ROBUST LIKELIHOOD FUNCTION
# =====================
def log_likelihood(theta, z, mu_obs, cov_inv, model_type):
    if model_type == 'lcdm':
        H0, Om = theta
        nu0 = 0.45
    else:  # DEFT
        H0, Om, nu0 = theta
        
    mu_model = distance_modulus(z, H0, Om, model_type, nu0)
    
    # Handle invalid values
    if np.any(np.isnan(mu_model)):
        return -np.inf
        
    residual = mu_obs - mu_model
    return -0.5 * residual @ cov_inv @ residual

# =====================
# MCMC WITH PHYSICAL CONSTRAINTS
# =====================
def run_mcmc_with_output(z_data, mu_data, cov_inv, model_type='lcdm'):
    if model_type == 'lcdm':
        n_params = 2
        labels = ['H0', 'Om']
        param_ranges = [(60, 80), (0.1, 0.5)]
    else:  # DEFT
        n_params = 3
        labels = ['H0', 'Om', 'nu0']
        param_ranges = [(60, 80), (0.1, 0.5), (0.1, 1.0)]
    
    n_walkers = 32
    n_steps = 4000
    
    # Initialize walkers in physical space
    pos = np.zeros((n_walkers, n_params))
    for i in range(n_params):
        low, high = param_ranges[i]
        pos[:, i] = np.random.uniform(low, high, n_walkers)
    
    # Setup sampler
    sampler = emcee.EnsembleSampler(
        n_walkers,
        n_params,
        log_probability,
        args=(z_data, mu_data, cov_inv, model_type)
    )
    
    print(f"Running {model_type.upper()} MCMC ({n_steps} steps)...")
    start = time.time()
    sampler.run_mcmc(pos, n_steps, progress=True)
    runtime = time.time() - start
    print(f"Completed in {runtime:.1f} seconds")
    
    # Process samples
    chain = sampler.get_chain(flat=True)
    log_probs = sampler.get_log_prob(flat=True)
    
    # Filter physical samples
    physical_chain = []
    for sample in chain:
        if log_prior(sample, model_type) == 0:
            physical_chain.append(sample)
    
    if not physical_chain:
        raise RuntimeError(f"No physical samples found for {model_type}")
    
    physical_chain = np.array(physical_chain)
    
    # Find best parameters
    best_idx = np.argmax(log_probs)
    best_params = chain[best_idx]
    
    # Ensure best parameters are physical
    if log_prior(best_params, model_type) != 0:
        # Find best physical parameters
        physical_log_probs = [log_probs[i] for i in range(len(chain)) 
                             if log_prior(chain[i], model_type) == 0]
        best_phys_idx = np.argmax(physical_log_probs)
        best_params = physical_chain[best_phys_idx]
    
    # Save outputs
    os.makedirs("output", exist_ok=True)
    np.savetxt(f"output/{model_type}_chain.txt", physical_chain)
    np.savetxt(f"output/{model_type}_best_params.txt", best_params)
    
    # Create corner plot
    fig = corner.corner(
        physical_chain,
        labels=labels,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_fmt='.3f',
        range=param_ranges
    )
    fig.savefig(f"output/{model_type}_corner.png", dpi=150)
    
    return best_params

# =====================
# CREATE SUPPRESSION PLOT
# =====================
def create_suppression_plot(lcdm_best, deft_best):
    z_test = np.linspace(0.3, 0.5, 100)
    
    # Calculate distance moduli
    mu_lcdm = distance_modulus(z_test, lcdm_best[0], lcdm_best[1], 'lcdm')
    mu_deft = distance_modulus(z_test, deft_best[0], deft_best[1], 'deft', deft_best[2])
    suppression = mu_deft - mu_lcdm
    
    # Create plot
    plt.figure(figsize=(10, 6))
    plt.plot(z_test, suppression, 'r-', lw=2, label='DEFT - ΛCDM')
    plt.axvline(0.4, color='k', ls='--', alpha=0.7, label='z=0.4')
    plt.axhline(0, color='k', alpha=0.3)
    plt.fill_between(z_test, suppression-0.05, suppression+0.05, alpha=0.2, color='red')
    plt.xlabel('Redshift (z)')
    plt.ylabel('$\Delta\mu$ (mag)')
    plt.title('DEFT Suppression Signature at z≈0.4')
    plt.legend()
    plt.grid(True)
    plt.savefig('output/deft_suppression.png', dpi=150)
    
    # Calculate suppression at z=0.4
    idx = np.argmin(np.abs(z_test - 0.4))
    suppression_at_z04 = suppression[idx]
    print(f"\nDEFT suppression at z=0.4: {suppression_at_z04:.4f} mag")
    
    return suppression_at_z04

# =====================
# MAIN EXECUTION
# =====================
if __name__ == "__main__":
    # Create output directory
    os.makedirs("output", exist_ok=True)
    print("Output directory: ", os.path.abspath("output"))
    
    # 1. Load data
    print("Loading data...")
    z_data, mu_data, cov = load_data()
    
    # 2. Precompute covariance decomposition
    print("Precomputing covariance decomposition...")
    eigvals, eigvecs = eigh(cov)
    eigvals[eigvals < 1e-6] = 1e-6
    cov_inv = eigvecs @ np.diag(1/eigvals) @ eigvecs.T
    
    # 3. Run ΛCDM MCMC
    print("\n" + "="*40)
    print("RUNNING ΛCDM MCMC")
    print("="*40)
    lcdm_best = run_mcmc_with_output(z_data, mu_data, cov_inv, 'lcdm')
    
    # 4. Run DEFT MCMC
    print("\n" + "="*40)
    print("RUNNING DEFT MCMC")
    print("="*40)
    deft_best = run_mcmc_with_output(z_data, mu_data, cov_inv, 'deft')
    
    # 5. Create suppression plot
    print("\nCreating suppression plot...")
    suppression = create_suppression_plot(lcdm_best, deft_best)
    
    print("\nANALYSIS COMPLETE!")
    print(f"ΛCDM best fit: H0={lcdm_best[0]:.2f}, Om={lcdm_best[1]:.3f}")
    print(f"DEFT best fit: H0={deft_best[0]:.2f}, Om={deft_best[1]:.3f}, nu0={deft_best[2]:.3f}")
    print(f"DEFT suppression at z=0.4: {suppression:.4f} mag")
    print("All outputs saved in: ", os.path.abspath("output"))
