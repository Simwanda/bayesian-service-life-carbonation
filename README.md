# Bayesian Updating of Carbonation Ingress Models using MATLAB/UQLab

This repository contains MATLAB/UQLab scripts for Bayesian updating of carbonation ingress models for reinforced concrete cooling towers.

The code supports:

- fib Bulletin 34 carbonation model
- Malami-type semi-mechanistic carbonation model
- Prior distributions based on Table 2 of the paper
- Bayesian updating using MCMC
- Prior vs posterior comparison plots
- Posterior correlation matrices
- Export of posterior samples for further reliability analysis

## Main scripts

Run the Bayesian updating:

```matlab
run_bayes_update_table2_fib_malami
