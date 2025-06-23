import numpy as np
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
import pandas as pd
from scipy.optimize import minimize

def deft_chi2(params, z_data, mu_obs, sigma_mu):
    """
    Compute DEFT/GRACE model χ² for JLA supernova data with numerical safeguards
    and proper error weighting.
    
    Parameters:
    params : list or array
        Model parameters: [gamma, z0, sigma, H0, Omega_m]
    z_data : array
        JLA supernova redshifts
    mu_obs : array
        Observed distance moduli
    sigma_mu : array
        Uncertainties in distance moduli
        
    Returns:
    chi2 : float
        Properly weighted χ² value
    """
    gamma, z0, sigma, H0, Omega_m = params
    
    c_light = 299792.458  # km/s
    h = H0 / 100.0
    
    # High-res redshift grid for integration and stability
    z_max = np.max(z_data)*1.2
    z_grid = np.linspace(1e-5, z_max, 2000)  # avoid zero to prevent division by zero
    
    # Gaussian entropic correction
    chi_z = gamma * np.exp(-0.5 * ((z_grid - z0)/sigma)**2)
    
    # Compute the integral in exponent: I_exp = ∫₀^z (1 + χ(z'))/(1+z') dz'
    integrand_exp = (1 + chi_z) / (1 + z_grid)
    I_exp = cumtrapz(integrand_exp, z_grid, initial=0)
    
    # Clip exponential argument for numerical stability
    exp_arg = 3 * I_exp
    max_exp_val = 100
    exp_arg_clipped = np.minimum(exp_arg, max_exp_val)
    
    # Compute normalized Hubble parameter squared
    E_z_sq = Omega_m * (1 + z_grid)**3 + (1 - Omega_m) * np.exp(exp_arg_clipped)
    # Ensure positivity to avoid sqrt of negative
    E_z_sq = np.maximum(E_z_sq, 1e-10)
    E_z = np.sqrt(E_z_sq)
    
    # Luminosity distance integral: I_dL = ∫₀^z dz'/E(z')
    integrand_dL = 1.0 / E_z
    I_dL = cumtrapz(integrand_dL, z_grid, initial=0)
    
    # Calculate luminosity distance in Mpc
    dL_grid = (1 + z_grid) * (c_light / H0) * I_dL
    
    # Clip to avoid log10(0) issues
    dL_grid_safe = np.clip(dL_grid, 1e-8, None)
    
    # Distance modulus μ = 5 log10(dL) + 25
    mu_grid = 5 * np.log10(dL_grid_safe) + 25
    
    # Interpolate model μ to data redshifts
    mu_interp = interp1d(z_grid, mu_grid, kind='cubic', bounds_error=False, fill_value='extrapolate')
    mu_model = mu_interp(z_data)
    
    residuals = mu_model - mu_obs
    chi2 = np.sum((residuals / sigma_mu)**2)
    
    return chi2

def load_jla_data(filename="jla_mub.txt"):
    # Assumes jla_mub.txt format: z mu sigma_mu
    df = pd.read_csv(filename, delim_whitespace=True, comment='#', names=['z', 'mu', 'sigma_mu'])
    return df['z'].values, df['mu'].values, df['sigma_mu'].values

if __name__ == "__main__":
    # Load cleaned JLA data
    z_data, mu_obs, sigma_mu = load_jla_data("jla_mub.txt")
    print(f"Loaded {len(z_data)} JLA supernovae.")
    
    # Initial guess for parameters: gamma, z0, sigma, H0, Omega_m
    p0 = [0.15, 0.4, 0.15, 70.0, 0.3]
    bounds = [(0, 1), (0.2, 0.7), (0.05, 0.3), (65, 80), (0.1, 0.5)]
    
    # Minimize chi2
    result = minimize(deft_chi2, p0, args=(z_data, mu_obs, sigma_mu), bounds=bounds, method='L-BFGS-B')
    
    if result.success:
        gamma_fit, z0_fit, sigma_fit, H0_fit, Omega_m_fit = result.x
        chi2_fit = result.fun
        dof = len(z_data) - len(p0)
        print(f"Fit results:")
        print(f"  gamma0 = {gamma_fit:.4f}")
        print(f"  z0     = {z0_fit:.4f}")
        print(f"  sigma  = {sigma_fit:.4f}")
        print(f"  H0     = {H0_fit:.2f} km/s/Mpc")
        print(f"  Omega_m= {Omega_m_fit:.4f}")
        print(f"  Chi^2  = {chi2_fit:.2f}")
        print(f"  dof    = {dof}")
        print(f"  Reduced Chi^2 = {chi2_fit/dof:.2f}")
    else:
        print("Fit failed:", result.message)
