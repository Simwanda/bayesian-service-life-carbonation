% ==========================================================
% run_bayes_update_fib_malami.m
% Bayesian updating of (A, n, theta_xc) using MCMC
% Observation model: y = theta * A * t^(0.5 - n) + eps
% ==========================================================
clear; clc; close all;

% -----------------------------
% USER: Provide observation data
% -----------------------------
% Option A (recommended): actual measurement vector
% t_obs = 42 * ones(Nobs,1);
% y_obs = [...]; % mm

% Option B: only mean & COV available -> create pseudo-sample
t0 = 42;                 % years
y_mean = 14.8;           % mm (example)
y_cov  = 0.36;           % example
Nobs   = 30;             % choose a reasonable count (e.g., number of tests)
rng(123,'twister');
y_obs  = max(0, y_mean + y_mean*y_cov*randn(Nobs,1)); % pseudo data (non-negative)
t_obs  = t0 * ones(Nobs,1);

% Measurement / residual sd (mm): if unknown, set to fraction of mean
sigma_eps = 1.5; % mm (tune; can also be y_mean*y_cov*0.3 etc.)

% -----------------------------
% Priors from your stochastic model
% -----------------------------
prior_common.n.mean   = 0.11;
prior_common.n.cov    = 0.706;     % lognormal in your model
prior_common.theta.mean = 1.0;
prior_common.theta.cov  = 0.20;    % often treated lognormal; we do lognormal to keep theta>0

prior_fib.A.mean    = 2.37;
prior_fib.A.cov     = 0.29;

prior_malami.A.mean = 2.35;
prior_malami.A.cov  = 0.27;

% -----------------------------
% MCMC settings
% -----------------------------
mcmc.Nburn = 20000;
mcmc.Nkeep = 40000;
mcmc.thin  = 5;

% proposal std in transformed space (log-params)
% tune these to get ~20-35% acceptance
mcmc.propStd = [0.08, 0.20, 0.10]; % [logA, logn, logtheta]

% -----------------------------
% Run updating for fib and Malami
% -----------------------------
post_fib = bayes_update_A_n_theta(prior_fib, prior_common, t_obs, y_obs, sigma_eps, mcmc);
post_mal = bayes_update_A_n_theta(prior_malami, prior_common, t_obs, y_obs, sigma_eps, mcmc);

% -----------------------------
% Summaries
% -----------------------------
fprintf('\nPosterior summaries (fib):\n');
disp(posterior_summary(post_fib));

fprintf('\nPosterior summaries (Malami):\n');
disp(posterior_summary(post_mal));

% -----------------------------
% OPTIONAL: Fit posterior moments for UQLab marginals
% (lognormal fit via mean/cov for A,n,theta)
% -----------------------------
fit_fib = fit_lognormal_moments(post_fib);
fit_mal = fit_lognormal_moments(post_mal);

save('posterior_samples_fib.mat','post_fib','fit_fib');
save('posterior_samples_malami.mat','post_mal','fit_mal');

% -----------------------------
% OPTIONAL: create posterior-updated UQLab Input objects
% (for reliability analysis afterwards)
% -----------------------------
% uqlab;
% myInputFib_post = buildInput_carbonation_posterior('fib','caseII',fit_fib);
% myInputMal_post = buildInput_carbonation_posterior('malami','caseII',fit_mal);
% % then run your MCS/FORM pipeline using these inputs

