# FRACE: Field-Regulated Acceleration from Cosmic Entanglement

Welcome to the official repository for the FRACE+ cosmology framework.

FRACE+ is a scalar-field-driven model designed to explain late-time cosmic acceleration through a decoherence phase transition at redshift \( z \approx 0.4 \). This repository contains the code and analysis used to identify a statistically significant suppression feature in the Pantheon+ Type Ia supernovae dataset, and to reconstruct a minimal scalar field \( \chi(z) = \gamma \log(1+z) \) capable of resolving the Hubble and \( S_8 \) tensions without altering early-universe physics.

---

## 🔍 Key Papers

- **Paper I:** *A Suppression Scar at \( z \approx 0.4 \): Evidence for a Decoherence Boundary in the FRACE+ Cosmology*  
  Focuses on the anomaly in Hubble residuals and presents the initial observational evidence.

- **Paper II:** *GRACE+: A Minimal Scalar Field Reconstruction for Resolving Late-Time Cosmological Tensions*  
  Introduces the theoretical mechanism and validates the scalar field dynamics through observational consistency checks.

---

## 📂 Repository Contents

- `fit_grace_stable_to_jla.py` – Scalar field fitting to Pantheon+ data  
- `deft_mcmc_analysis.py` – DEFT suppression feature MCMC analysis  
- `make_mock_jla.py` – Generator for mock SNe samples under ΛCDM  
- `jla_clean_processor.py` – Preprocessor for real SNe datasets  
- `plots/` – All residual and suppression feature visualizations

---

## 📊 Key Figures

Below are visual summaries from the core analyses included in this repository.

**1. χ²(z₀) Confidence Sweep**  
Best-fit rupture redshift \( z_0 = 0.01 \) with confidence intervals for Δχ² = 1, 4, 9  
![z0_chi2_confidence_sweep](Plots/z0_chi2_confidence_sweep.png)

**2. H₀ and ν₀ Posterior (Pantheon+ DEFT fit)**  
MCMC contours for Hubble constant and memory coupling  
![corner_deft_H0_nu0_fixedOm](Plots/corner_deft_H0_nu0_fixedOm.png)

**3. DEFT Residual Suppression Signature**  
Magnitude residuals from ΛCDM baseline showing suppression at \( z \approx 0.4 \)  
![deft_residual_suppression_feature_z0.4](Plots/deft_residual_suppression_feature_z0.4.png)

**4. Pantheon+ + DESI Joint Constraints on \( \chi(z) = \gamma \log(1+z) \)**  
Corner plot for scalar field amplitude γ and normalization χ₀  
![pantheon_desi_corner_chi0_gamma](Plots/pantheon_desi_corner_chi0_gamma.png)

---

## 🚀 Getting Started

Clone the repo and install dependencies:

```bash
git clone https://github.com/WhitWhitman/FRACE.git
cd FRACE
pip install -r requirements.txt
```

Run the main analysis pipeline:

```bash
python fit_grace_stable_to_jla.py
```

---

## 🧠 Suggested Readings

If you're new to this framework or cosmological phase transitions:

- [Pantheon+ Collaboration Data](https://github.com/PantheonPlusSH0ES)
- Verde, L., Treu, T., & Riess, A. (2019). Tensions between the early and late Universe.
- Carroll, S. M. (2019). *Spacetime and Geometry: An Introduction to General Relativity*.
- Original DEFT model documentation (available on request).

---

## 🧠 Contact

This project is led by Ken Whitman.  
For collaborations, corrections, or questions, please open an issue or contact [Ken via GitHub](https://github.com/WhitWhitman).
