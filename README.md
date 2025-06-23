# FRACE: Field-Regulated Acceleration from Cosmic Entanglement

Welcome to the official repository for the FRACE+ cosmology framework.

FRACE+ is a scalar-field-driven model designed to explain late-time cosmic acceleration through a decoherence phase transition at redshift \( z \approx 0.4 \). This repository contains the code, figures, and analysis used to identify a statistically significant suppression feature in the Pantheon+ Type Ia supernovae dataset and to reconstruct a minimal scalar field \( \chi(z) = \gamma \log(1+z) \) capable of resolving the Hubble and \( S_8 \) tensions without altering early-universe physics.

---

## 🔍 Key Papers

- **Paper I:** *A Suppression Scar at \( z \approx 0.4 \): Evidence for a Decoherence Boundary in the FRACE+ Cosmology*  
  Focuses on the anomaly in Hubble residuals and presents the initial observational evidence.

- **Paper II:** *GRACE+: A Minimal Scalar Field Reconstruction for Resolving Late-Time Cosmological Tensions*  
  Introduces the theoretical mechanism and validates the scalar field dynamics through observational consistency checks.

---

## 📄 Download Papers

You can view or download the formal PDF drafts of each paper below:

- [Paper I – FRACE+ Suppression Scar Detection](papers/Paper_I_FRACE_SuppressionScar_z0.4.pdf)  
- [Paper II – GRACE+ Scalar Field Reconstruction](papers/Paper_II_GRACE_ScalarField_Model.pdf)

---

## 📂 Repository Contents

- `fit_grace_stable_to_jla.py` – Scalar field fitting to Pantheon+ data  
- `deft_mcmc_analysis.py` – DEFT suppression feature analysis  
- `make_mock_jla.py` – Generator for mock SNe samples under ΛCDM  
- `jla_clean_processor.py` – Preprocessor for real SNe datasets  
- `deft_joint_fit_gamma4.py` – Joint fit with gamma = 4 fixed  
- `deft_joint_fit_desi_gamma4.py` – Joint DESI + Pantheon+ fit  
- `plots/` – Visualizations, corner plots, suppression diagnostics  
- `papers/` – Official PDF drafts of FRACE+ and GRACE+ papers  

---

## 🖼 Sample Figures (from `/plots`)

- ![Suppression Residual Feature](plots/deft_residual_suppression_feature_z0.4.png)  
- ![Corner Plot (H₀, ν₀, Ωₘ)](plots/corner_deft_H0_nu0_fixedOm.png)
- ![χ² Confidence Sweep](plots/z0_chi2_confidence_sweep.png)  
- ![Joint Corner: Pantheon+ + DESI](plots/pantheon_desi_corner_chi0_gamma.png)

---

## 🚀 Getting Started

Clone the repo and install dependencies:

```bash
git clone https://github.com/WhitWhitman/FRACE.git
cd FRACE
pip install -r requirements.txt
