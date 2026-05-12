% ==========================================================
% run_bayes_update_fib_malami.m
% Bayesian updating of theta_x ONLY
% Observation model:
%   y = theta_x * A0 * t^(0.5 - n0) + eps
% ==========================================================
clear; clc; close all;

% -----------------------------
% Observation data
% -----------------------------
t0 = 42;                 % years
y_mean = 14.8;           % mm
y_cov  = 0.36;
Nobs   = 30;

rng(123,'twister');
y_obs  = max(0, y_mean + y_mean*y_cov*randn(Nobs,1));
t_obs  = t0 * ones(Nobs,1);

sigma_eps = 5.3; % mm

% -----------------------------
% Fixed deterministic parameters
% -----------------------------
A_fib    = 2.37;
A_malami = 2.35;
n0       = 0.11;

% -----------------------------
% Prior for theta_x (Gaussian)
% -----------------------------
prior_common.theta.mean = 1.0;
prior_common.theta.cov  = 0.20;   % Gaussian

% -----------------------------
% MCMC settings
% -----------------------------
mcmc.Nburn = 20000;
mcmc.Nkeep = 40000;
mcmc.thin  = 5;
mcmc.propStd_theta = 0.05;  % RW proposal (tune ~20–30% acc.)

% -----------------------------
% Run Bayesian updating
% -----------------------------
post_fib = bayes_update_theta_x( ...
    prior_common, A_fib, n0, t_obs, y_obs, sigma_eps, mcmc);

post_mal = bayes_update_theta_x( ...
    prior_common, A_malami, n0, t_obs, y_obs, sigma_eps, mcmc);

% -----------------------------
% Summaries
% -----------------------------
disp('Posterior (fib):');
disp(posterior_summary(post_fib));

disp('Posterior (Malami):');
disp(posterior_summary(post_mal));

% -----------------------------
% Save
% -----------------------------
save('posterior_samples_fib.mat','post_fib');
save('posterior_samples_malami.mat','post_mal');
