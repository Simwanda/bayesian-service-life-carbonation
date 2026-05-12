# Bayesian updating code update for Section 4.4.2

This folder contains an updated MATLAB/UQLab-compatible workflow responding to reviewer comments R2.4 and R3.1.

## What changed

The previous workflow updated only the multiplicative model uncertainty factor `theta_xc`. The new workflow defines the probabilistic prior models from Table 2 and updates selected basic variables jointly with `theta_xc` using a random-walk Metropolis MCMC scheme.

Main new files:

- `table2_priors_carbonation.m` — defines all Table 2 stochastic priors for the fib Bulletin 34 and Malami models.
- `carbonation_xc_table2.m` — evaluates carbonation depth from the Table 2 basic variables.
- `bayes_update_table2_params.m` — generic MCMC updater for selected Table 2 variables.
- `run_bayes_update_table2_fib_malami.m` — main script to run the new reviewer-response update.
- `posterior_summary_table2.m` — posterior statistics for all updated variables.
- `fit_posterior_marginals_table2.m` — converts posterior samples to marginal distributions for UQLab.
- `buildInput_carbonation_posterior_table2.m` — creates a UQLab input object from fitted posterior marginals.
- `plot_xc_posterior_table2.m` — plots posterior predictive carbonation-depth curves.
- `section_4_4_2_replacement_text.tex` — draft replacement text for the manuscript.

The older theta-only scripts are kept for comparison, but the recommended script is:

```matlab
run_bayes_update_table2_fib_malami
```

## Important modelling note

With only one inspection campaign at 42 years, the posterior cannot uniquely identify all basic variables. Therefore, the update is intentionally prior-regularised: the prior distributions from Table 2 remain influential, while the inspection data shift the sensitive variables and the model discrepancy factor only to the extent supported by the likelihood.

If individual carbonation-depth measurements become available, replace the summary likelihood in `run_bayes_update_table2_fib_malami.m` with the raw vector `obs.y_mm`.

## Malami forward model update

The Malami option in `carbonation_xc_table2.m` has been updated to use the closed-form Malami implementation supplied by the authors, translated from Python to MATLAB. It now evaluates:

```matlab
D_c = sqrt(2*k_RH*k_T*k_c*CO2*R_NAC0_inv*t_days)*Wwet*1000;
```

with hydration-scaled CO2 binding capacity `a_CO2`, effective diffusivity `DeCO2`, `R_NAC0_inv = DeCO2/a_CO2`, humidity, temperature, curing, CO2 and time-of-wetness terms. The implementation still uses the Table 2 probabilistic variables as priors for Bayesian updating.
