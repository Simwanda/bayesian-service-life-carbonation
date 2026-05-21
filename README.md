# Bayesian Carbonation Updating with MATLAB/UQLab

MATLAB/UQLab scripts for Bayesian updating of carbonation-ingress models for reinforced concrete cooling towers using in-situ inspection data.

The repository supports two carbonation modelling frameworks:

- **fib Bulletin 34 resistance-based model**
- **Malami-type semi-mechanistic carbonation model**

The workflow was prepared to support manuscript revisions responding to reviewer comments on Bayesian posterior modelling. The current recommended workflow updates selected probabilistic variables from Table 2 jointly with a **multiplicative model-uncertainty factor** rather than updating only a single correction factor.

## Main features

- Prior distributions based on Table 2 of the related paper
- Bayesian updating using random-walk Metropolis MCMC
- Joint updating of selected material, exposure, and model-uncertainty variables
- Multiplicative model uncertainty through `theta_xc`
- Optional additive-discrepancy workflow through `delta_xc`
- Loop over individual carbonation-depth measurements
- Posterior predictive carbonation-depth curves
- Combined posterior samples from multiple individual measurement updates
- Prior-versus-posterior distribution comparison tables
- Histogram/PDF plots of posterior samples and fitted distributions
- Posterior Pearson and Spearman correlation matrices
- Empirical posterior CDF export for bounded variables such as `tc` and `pSR`
- MATLAB/UQLab-compatible posterior marginal definitions for reliability analysis
- Figure 8-style posterior carbonation-depth plots
- Beamer presentation material for reviewing/post-processing results

## Recommended workflow

The recommended workflow is the **multiplicative model-uncertainty workflow**:

```text
x_c,updated(t) = theta_xc * x_c,model(t, X)
```

where `theta_xc` is updated jointly with selected probabilistic variables from Table 2.

The additive discrepancy workflow is retained for comparison and for UQLab-style additive residual formulations, but the multiplicative workflow is more consistent with the carbonation model uncertainty used in the manuscript.

## Repository structure

| File or folder | Purpose |
|---|---|
| `run_bayes_update_table2_fib_malami.m` | Main single-update script for Bayesian updating of fib and Malami models using summary inspection data |
| `run_loop_posterior_all_measurements_multiplicative.m` | Recommended loop script: runs multiplicative Bayesian updating for each individual carbonation-depth measurement |
| `plot_combined_posterior_histograms_with_fits_multiplicative.m` | Plots histograms of combined posterior samples and overlays fitted Normal/Lognormal PDFs |
| `compare_prior_posterior_distribution_types_multiplicative.m` | Compares Table 2 prior PDFs with fitted posterior PDFs |
| `fit_empirical_distributions_pSR_tc_multiplicative.m` | Fits empirical/KDE posterior distributions for bounded variables `tc` and `pSR` |
| `export_empirical_CDF_pSR_tc_for_FReET.m` | Exports `x` and `F(x)` empirical CDF tables for `tc` and `pSR` for use in FReET |
| `table2_priors_carbonation.m` | Defines prior distributions from Table 2 |
| `carbonation_xc_table2.m` | Forward carbonation-depth model for fib and Malami options |
| `bayes_update_table2_params.m` | Random-walk Metropolis MCMC updater |
| `posterior_summary_table2.m` | Prints posterior summaries |
| `fit_posterior_marginals_table2.m` | Fits posterior marginals for UQLab input construction |
| `buildInput_carbonation_posterior_table2.m` | Creates UQLab posterior input object |
| `compare_prior_posterior_all_variables.m` | Plots prior vs posterior for all variables and exports correlation matrices |
| `produce_figure8_posterior_carbonation.m` | Produces Figure 8-style posterior carbonation-depth plot |
| `section_4_4_2_replacement_text.tex` | Draft manuscript replacement text for Section 4.4.2 |
| `MULTIPLICATIVE_LOOP/` or similar folder | Multiplicative-loop helper scripts, depending on the local repository version |
| `ADDITIVE_DISCREPANCY/` or similar folder | Optional additive-discrepancy helper scripts, depending on the local repository version |

Older theta-only scripts and additive-discrepancy scripts are retained for comparison. For the current manuscript response, the multiplicative Table 2 parameter-updating workflow is the preferred option.

## Requirements

- MATLAB
- UQLab
- Statistics and Machine Learning Toolbox recommended

The code uses MATLAB functions such as `fitdist`, `ksdensity`, `ecdf`, `corr`, `writetable`, and `exportgraphics`.

## Quick start: single posterior update

Run the main Bayesian updating workflow using summary inspection data:

```matlab
clear; clc; close all;
addpath(genpath(pwd));

run_bayes_update_table2_fib_malami
```

Compare prior and posterior distributions and produce posterior correlation matrices:

```matlab
compare_prior_posterior_all_variables
```

Produce the Figure 8-style posterior carbonation-depth plot:

```matlab
produce_figure8_posterior_carbonation
```

## Recommended run: multiplicative loop over individual measurements

To update the model using each individual carbonation-depth measurement and then combine posterior samples, run:

```matlab
clear; clc; close all;
addpath(genpath(pwd));

run_loop_posterior_all_measurements_multiplicative
```

This script uses the individual carbonation-depth observations:

```text
26, 19.1, 19.8, 9, 5.7, 18.2, 22.6, 18.6, 16.1, 17,
14.2, 11, 16.4, 11, 13.7, 12.2, 13.9, 12.9, 6.1, 15.6,
13.3, 12.6, 5, 21.2, 17.3, 15.7, 15.1, 12.6, 24.4, 26.5,
15.4, 13, 14.2, 21.5, 18.9, 15.5, 11.2, 21.5, 11.5, 21.2,
9.9, 22.1, 6.2, 7.3, 10.9, 8.6, 14.1, 7, 12, 14.4
```

The observation model used in the loop is typically:

```matlab
obs.y_std     = 0.0;
obs.N_eff     = 1;
obs.sigma_eps = 1.0;
```

Here, `sigma_eps = 1.0 mm` represents the uncertainty of the calibration target for each individual update. The spatial scatter of the full measurement set is retained separately for interpretation and plotting.

## Outputs from the multiplicative loop

The multiplicative loop saves outputs in:

```text
figs_multiplicative_loop/
```

Typical generated files include:

```text
figs_multiplicative_loop/combined_loop_results_multiplicative.mat
figs_multiplicative_loop/loop_posterior_predictive_summary_multiplicative.csv
figs_multiplicative_loop/all_posterior_mean_curves_fib_multiplicative.pdf
figs_multiplicative_loop/all_posterior_mean_curves_malami_multiplicative.pdf
figs_multiplicative_loop/combined_posterior_predictive_fib_multiplicative.pdf
figs_multiplicative_loop/combined_posterior_predictive_malami_multiplicative.pdf
figs_multiplicative_loop/combined_fit_summary_fib_multiplicative.csv
figs_multiplicative_loop/combined_fit_summary_malami_multiplicative.csv
figs_multiplicative_loop/combined_corr_fib_pearson_multiplicative.csv
figs_multiplicative_loop/combined_corr_fib_spearman_multiplicative.csv
figs_multiplicative_loop/combined_corr_malami_pearson_multiplicative.csv
figs_multiplicative_loop/combined_corr_malami_spearman_multiplicative.csv
```

## Posterior distribution fitting

After the multiplicative loop, fit and plot posterior marginal distributions using:

```matlab
plot_combined_posterior_histograms_with_fits_multiplicative
```

This creates histogram/PDF plots for each combined posterior variable and saves summary tables such as:

```text
figs_multiplicative_loop/posterior_histogram_fit_summary_fib_multiplicative.csv
figs_multiplicative_loop/posterior_histogram_fit_summary_malami_multiplicative.csv
```

The fitting script currently compares:

- Normal distribution
- Lognormal distribution, only for strictly positive samples

The preferred fit is selected using AIC.

## Prior-versus-posterior distribution comparison

To compare Table 2 prior distributions with the fitted posterior distributions, run:

```matlab
compare_prior_posterior_distribution_types_multiplicative
```

This exports:

```text
figs_multiplicative_loop/prior_vs_posterior_distribution_types_fib_multiplicative.csv
figs_multiplicative_loop/prior_vs_posterior_distribution_types_malami_multiplicative.csv
figs_multiplicative_loop/prior_vs_posterior_distribution_types_ALL_multiplicative.csv
```

These files are useful for transferring the posterior stochastic model into reliability software such as FReET.

## Empirical posterior distributions for `tc` and `pSR`

The variables `tc` and `pSR` are bounded:

```text
tc  in [1, 3] days
pSR in [0.05, 0.25]
```

Because these variables may not be well represented by Normal or Lognormal distributions, empirical posterior distributions can be used.

Run:

```matlab
fit_empirical_distributions_pSR_tc_multiplicative
```

This produces empirical PDF/CDF files and plots in:

```text
figs_multiplicative_loop/empirical_fit_pSR_tc/
```

Typical outputs include:

```text
fib_pSR_empirical_pdf_cdf.csv
fib_tc_empirical_pdf_cdf.csv
malami_pSR_empirical_pdf_cdf.csv
malami_tc_empirical_pdf_cdf.csv
empirical_fit_summary_ALL.csv
```

## Export empirical CDF files for FReET

If FReET supports empirical or user-defined distributions, export clean `x` and `F(x)` files using:

```matlab
export_empirical_CDF_pSR_tc_for_FReET
```

This creates:

```text
figs_multiplicative_loop/empirical_CDF_for_FReET/fib_pSR_CDF_for_FReET.csv
figs_multiplicative_loop/empirical_CDF_for_FReET/fib_tc_CDF_for_FReET.csv
figs_multiplicative_loop/empirical_CDF_for_FReET/malami_pSR_CDF_for_FReET.csv
figs_multiplicative_loop/empirical_CDF_for_FReET/malami_tc_CDF_for_FReET.csv
```

Each file contains:

```text
x,Fx
```

where `x` is the variable value and `Fx` is the empirical cumulative probability.

If direct empirical distributions are not supported in FReET, use one of the following approximations:

1. scaled beta distribution on the physical bounds;
2. truncated normal distribution on the physical bounds;
3. uniform distribution between posterior 2.5% and 97.5% quantiles.

For the current results, empirical distributions are mainly recommended for `tc` and `pSR`.

## Additive discrepancy workflow, optional

An additive-discrepancy version is retained for comparison with the standard UQLab observation model:

```text
x_c,updated(t) = x_c,model(t, X) + delta_xc
```

Typical additive scripts are:

```matlab
check_additive_delta_xc_setup
run_bayes_update_table2_fib_malami_additive
plot_xc_posterior_table2_additive
run_loop_posterior_all_measurements_additive
```

The additive discrepancy variable is `delta_xc` and must be present in `table2_priors_carbonation.m` for both `fib` and `malami`. This workflow is optional and should be treated as a computational comparison, not as the main manuscript workflow unless explicitly chosen.

## Modelling notes

### Malami forward model

The Malami option is implemented as a reduced closed-form semi-mechanistic forward model, not as a full finite-difference diffusion-reaction PDE solver. The implementation evaluates carbonation depth using humidity, temperature, curing, CO2 concentration, hydration-scaled CO2 binding capacity, effective CO2 diffusivity, inverse carbonation resistance, and time-of-wetness terms.

### Measurement error and model uncertainty

Two uncertainty components are distinguished:

- `theta_xc`: multiplicative model-uncertainty/model-discrepancy factor in the carbonation model;
- `sigma_eps`: observation-error term used only in the likelihood.

In the multiplicative workflow:

```text
x_c,updated(t) = theta_xc * x_c,model(t, X)
```

and the likelihood compares the predicted carbonation depth with the observed carbonation depth using `sigma_eps`.

The value `sigma_eps = 1.0 mm` is interpreted as uncertainty of the calibration target, while the full scatter of measured carbonation depths is interpreted as spatial variability across the structure.

### Identifiability

Because the available case-study information is concentrated at one inspection time, posterior updating is prior-regularised. The prior distributions remain influential, and the data update only those variables that can be meaningfully informed by the available inspection evidence. When the loop over individual measurements is used, posterior samples from all individual updates can be combined to represent the field-informed posterior stochastic model.

## Recommended `.gitignore`

Consider adding the following to `.gitignore` if you do not want to version generated files:

```text
*.asv
*.mat
figs/
figs_multiplicative_loop/
results/
```

Keep `.mat`, `figs/`, and `figs_multiplicative_loop/` only if you want to share posterior samples and generated figures.

## Citation

If you use this code, please cite the related paper. Until the final journal details are available, use the following provisional citation and update it after publication:

```bibtex
@unpublished{somodikova2026carbonation,
  title  = {Comparative assessment of carbonation-ingress modelling in reinforced concrete structures with application to a cooling tower},
  author = {Somod\'{i}kov\'{a}, Martina and Simwanda, Lenganji and S\'{y}kora, Miroslav and Lehk\'{y}, David},
  note   = {Manuscript under review / in revision},
  year   = {2026}
}
```

Please also cite UQLab when using the UQLab parts of this workflow:

```bibtex
@inproceedings{marelli2014uqlab,
  title     = {UQLab: A framework for uncertainty quantification in MATLAB},
  author    = {Marelli, Stefano and Sudret, Bruno},
  booktitle = {Vulnerability, Uncertainty, and Risk: Quantification, Mitigation, and Management},
  year      = {2014},
  publisher = {ASCE}
}
```

For the UQLab Bayesian inversion manual, cite:

```bibtex
@TechReport{UQdoc_20_113,
  author      = {Wagner, P.-R. and Nagel, J. and Marelli, S. and Sudret, B.},
  title       = {{UQLab user manual -- Bayesian inversion for model calibration and validation}},
  institution = {Chair of Risk, Safety and Uncertainty Quantification, ETH Zurich, Switzerland},
  year        = {2022},
  note        = {Report UQLab-V2.0-113}
}
```

## License

Add a license before making the repository public. For academic code, common choices are MIT, BSD-3-Clause, or GPL-3.0. If the code is connected to an unpublished manuscript or institutional project, confirm the preferred license with all co-authors before publication.
