import numpy as np
import pandas as pd

# 1) Pick a few redshifts (e.g. 50 points between 0.01 and 1.3)
z_mock = np.linspace(0.01, 1.3, 50)

# 2) Compute a toy distance modulus: for low z, 
#    dL ≈ (c/H0)*z * (1 + 0.5 z) gives a rough shape
H0 = 67.4
c = 299792.458  # km/s
dL_approx = (c/H0) * z_mock * (1 + 0.5 * z_mock)

# 3) Distance modulus μ = 5*log10(dL) + 25
mu_mock = 5 * np.log10(dL_approx) + 25

# 4) Assign a constant uncertainty, e.g. σ_μ = 0.1 mag
sigma_mu_mock = np.full_like(z_mock, 0.1)

# 5) Save as a tab-separated file with no header
mock_df = pd.DataFrame({
    'z': z_mock,
    'mu': mu_mock,
    'sigma_mu': sigma_mu_mock
})
mock_df.to_csv('mock_jla.txt', sep='\t', index=False, header=False)

# Print the first few lines so you can inspect
print(mock_df.head().to_string(index=False))
