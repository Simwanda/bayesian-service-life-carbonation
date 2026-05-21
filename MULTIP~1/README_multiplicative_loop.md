# Multiplicative Bayesian loop scripts

These scripts rerun the individual-measurement loop using the multiplicative model uncertainty factor `theta_xc`.

## Main script

```matlab
run_loop_posterior_all_measurements_multiplicative
```

This uses:

```matlab
obs.y_std = 0.0;
obs.sigma_eps = 1.0;
```

and updates `theta_xc` jointly with the selected Table 2 variables.

## Output folder

All outputs are saved in:

```text
figs_multiplicative_loop/
```

Important outputs:

```text
combined_loop_results_multiplicative.mat
loop_posterior_predictive_summary_multiplicative.csv
combined_fit_summary_fib_multiplicative.csv
combined_fit_summary_malami_multiplicative.csv
combined_corr_fib_spearman_multiplicative.csv
combined_corr_malami_spearman_multiplicative.csv
combined_corr_fib_pearson_multiplicative.csv
combined_corr_malami_pearson_multiplicative.csv
```

## Plot histogram + fitted PDFs

```matlab
plot_combined_posterior_histograms_with_fits_multiplicative
```

## Compare prior and posterior distribution types

```matlab
compare_prior_posterior_distribution_types_multiplicative
```

## Notes

This version does not use `delta_xc`. It returns to the multiplicative formulation:

```text
x_c,updated(t) = theta_xc * x_c,model(t, X)
```
