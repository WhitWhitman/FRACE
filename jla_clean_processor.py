    def test_grace_variants(self):
        """Linear GRACE fit with sensible, physical bounds and fixed H0"""

        H0_fixed = 67.40      # Planck value in km/s/Mpc
        c        = 299792.458 # speed of light

        def mu_linear_grace(z, Om, gamma0):
            """
            Distance modulus for linear GRACE
            χ(z) = γ0 * z / (1 + z)
            """
            # integrand for comoving distance with entropic correction
            def E_inv_corr(zp):
                chi = gamma0 * zp / (1.0 + zp)
                corr = np.exp(chi)
                return corr / np.sqrt(Om * (1. + zp)**3 + (1. - Om))

            mu = []
            for zi in np.atleast_1d(z):
                if zi <= 0.0:
                    mu.append(-np.inf)          # safeguard
                    continue
                dc, _ = quad(E_inv_corr, 0.0, zi)
                dL = (1. + zi) * (c / H0_fixed) * dc
                mu.append(5. * np.log10(dL) + 25.)
            return np.array(mu)

        def chi2_linear(params):
            Om, gamma0 = params
            if not (0.20 <= Om <= 0.40 and 0.0 <= gamma0 <= 5.0):
                return 1e30                    # out-of-bounds penalty
            mu_model = mu_linear_grace(self.df['z'].values, Om, gamma0)
            resid    = (self.df['mu'].values - mu_model) / self.df['sigma_mu'].values
            return np.sum(resid**2)

        # --- Run fit -------------------------------------------------------
        print("\n🔧  Fitting *Linear GRACE* with fixed H₀ = 67.40")
        init   = [0.30, 0.5]                   # Ωm, γ0
        bounds = [(0.20, 0.40), (0.0, 5.0)]

        result = minimize(chi2_linear, init, bounds=bounds, method='L-BFGS-B')

        if result.success:
            Om_best, gamma_best = result.x
            chi2 = result.fun
            dof  = len(self.df) - 2            # 2 free params (Ωm, γ0)
            print(f"   Ω_m   = {Om_best:6.3f}")
            print(f"   γ₀    = {gamma_best:6.3f}")
            print(f"   χ²    = {chi2:10.2f}")
            print(f"   dof   = {dof}")
            print(f"   χ²/dof= {chi2/dof:8.2f}")
        else:
            print("   ❌ optimiser failed:", result.message)
