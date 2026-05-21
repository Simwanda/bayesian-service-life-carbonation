# Bayesian Carbonation Updating with MATLAB/UQLab

MATLAB/UQLab scripts for Bayesian updating of carbonation-ingress models for reinforced concrete cooling towers using in-situ inspection data.

The repository supports two carbonation modelling frameworks:

- **fib Bulletin 34 resistance-based model**
- **Malami-type semi-mechanistic carbonation model**

The workflow was prepared to support manuscript revisions responding to reviewer comments on Bayesian posterior modelling, especially updating selected probabilistic input variables from Table 2 instead of updating only a single model-uncertainty factor.

## Main features

- Prior distributions based on Table 2 of the related paper
- Bayesian updating using random-walk Metropolis MCMC
- Joint updating of selected material, exposure, and model-uncertainty variables
- Posterior predictive carbonation-depth curves
- Prior-versus-posterior plots for all updated variables
- Posterior Pearson and Spearman correlation matrices
- MATLAB/UQLab-compatible posterior marginal definitions for reliability analysis
- Figure 8 reproduction script for posterior carbonation-depth evolution

## Repository structure

| File | Purpose |
|---|---|
| `run_bayes_update_table2_fib_malami.m` | Main script for Bayesian updating of fib and Malami models |
| `table2_priors_carbonation.m` | Defines prior distributions from Table 2 |
| `carbonation_xc_table2.m` | Forward carbonation-depth model for fib and Malami options |
| `bayes_update_table2_params.m` | Random-walk Metropolis MCMC updater |
| `posterior_summary_table2.m` | Prints posterior summaries |
| `fit_posterior_marginals_table2.m` | Fits posterior marginals for UQLab input construction |
| `buildInput_carbonation_posterior_table2.m` | Creates UQLab posterior input object |
| `compare_prior_posterior_all_variables.m` | Plots prior vs posterior for all variables and exports correlation matrices |
| `produce_figure8_posterior_carbonation.m` | Produces Figure 8-style posterior carbonation-depth plot |
| `section_4_4_2_replacement_text.tex` | Draft manuscript replacement text for Section 4.4.2 |

Older theta-only scripts are retained for comparison, but the recommended workflow is the Table 2 parameter-updating workflow.

## Requirements

- MATLAB
- UQLab
- Statistics and Machine Learning Toolbox recommended

## Quick start

Run the full Bayesian updating workflow:

```matlab
clear; clc;
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

## Outputs

The scripts save outputs in the `figs/` folder, including:

- prior/posterior distribution plots
- posterior predictive carbonation-depth plots
- posterior Pearson correlation matrices
- posterior Spearman correlation matrices
- posterior sample CSV files
- Figure 8-style PDF/PNG figures

Typical generated files include:

```text
figs/prior_posterior_fib_ALL.pdf
figs/prior_posterior_malami_ALL.pdf
figs/correlation_matrix_fib_pearson.csv
figs/correlation_matrix_fib_spearman.csv
figs/correlation_matrix_malami_pearson.csv
figs/correlation_matrix_malami_spearman.csv
figs/Figure8_posterior_carbonation_depth.pdf
figs/Figure8_posterior_carbonation_depth.png
```

## Modelling note

The Malami option is implemented as a reduced closed-form semi-mechanistic forward model, not as a full finite-difference diffusion-reaction PDE solver. The implementation evaluates carbonation depth using humidity, temperature, curing, CO2 concentration, hydration-scaled CO2 binding capacity, effective CO2 diffusivity, inverse carbonation resistance, and time-of-wetness terms.

Because the available case-study information contains one main inspection time at 42 years, posterior updating is prior-regularised. The prior distributions remain influential, and the data update only those variables that can be meaningfully informed by the available inspection evidence. If individual carbonation-depth observations become available, they can be inserted into the likelihood through `obs.y_mm` in `run_bayes_update_table2_fib_malami.m`.

## Recommended `.gitignore`

Consider adding the following to `.gitignore` if you do not want to version generated files:

```text
*.asv
*.mat
figs/
results/
```

Keep `.mat` and `figs/` only if you want to share posterior samples and generated figures.

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

## License

Add a license before making the repository public. For academic code, common choices are MIT, BSD-3-Clause, or GPL-3.0. If the code is connected to an unpublished manuscript or institutional project, confirm the preferred license with all co-authors before publication.

## Additive discrepancy and loop over individual measurements

For the UQLab-style additive discrepancy version, use:

```matlab
clear; clc; close all;
addpath(genpath(pwd));
check_additive_delta_xc_setup
run_bayes_update_table2_fib_malami_additive
plot_xc_posterior_table2_additive
```

To loop over the 50 individual carbonation-depth measurements and combine the posterior samples, run:

```matlab
clear; clc; close all;
addpath(genpath(pwd));
check_additive_delta_xc_setup
run_loop_posterior_all_measurements_additive
```

The additive discrepancy variable is `delta_xc` and must be present in `table2_priors_carbonation.m` for both `fib` and `malami`.
