%RUN_BAYES_UPDATE_TABLE2_FIB_MALAMI
% Bayesian updating of selected Table 2 probabilistic variables using MCMC.
% This is the reviewer-response version for Section 4.4.2:
%   - Table 2 stochastic variables are defined as prior distributions.
%   - The inspection data update selected basic variables and theta_xc jointly.
%   - The old theta-only script is kept in the package for comparison.
%
% NOTE on identifiability:
% With one inspection campaign, updating every uncertain variable is possible in
% a formal Bayesian sense but not all parameters are strongly identifiable. The
% posterior is therefore intentionally prior-regularised. The variable list below
% focuses on variables that are sensitive and have direct support in Table 2.

clear; clc; close all;

% -----------------------------
% Observation data from CT100/60
% -----------------------------
obs.t_years = 42;
obs.y_mean  = 14.78;  % mm
obs.y_std   = 5.3;    % mm, from CoV about 0.36
obs.N_eff   = 1;      % effective independent information; increase only if raw independent data are available
obs.sigma_eps = 0.0;  % additional model/measurement error on the mean, if required

% If actual individual measurements are available, replace the summary above by:
% obs.y_mm = [...];
% obs.t_years = 42*ones(numel(obs.y_mm),1);
% obs.sigma_eps = 1.0; % measurement noise [mm]

% -----------------------------
% MCMC settings
% -----------------------------
mcmc.Nburn = 20000;
mcmc.Nkeep = 40000;
mcmc.thin  = 5;
mcmc.rngSeed = 20260512;

% Proposal standard deviations in transformed variable space. These can be
% tuned to obtain an acceptance rate of about 0.20-0.50 for multi-parameter updates.

% -----------------------------
% fib Bulletin 34 update
% -----------------------------
update_fib = {'Racc_inv','kt','eps_t','tc','bc','CO2ppm','RH','tw','bw','pSR','theta_xc'};
mcmc.propStd = [0.07 0.07 0.06 0.08 0.05 0.05 0.05 0.07 0.06 0.08 0.08];
post_fib = bayes_update_table2_params('fib', obs, mcmc, update_fib);
fit_fib = fit_posterior_marginals_table2(post_fib);

% -----------------------------
% Malami update
% -----------------------------
update_malami = {'W_C','tc','bc','CO2ppm','RH','tw','bw','pSR','T','theta_xc'};
mcmc.propStd = [0.07 0.08 0.05 0.05 0.05 0.07 0.06 0.08 0.04 0.08];
post_malami = bayes_update_table2_params('malami', obs, mcmc, update_malami);
fit_malami = fit_posterior_marginals_table2(post_malami);

% -----------------------------
% Display and save
% -----------------------------
disp('Posterior summary: fib Bulletin 34'); disp(post_fib.summary);
disp('Posterior summary: Malami'); disp(post_malami.summary);

save('posterior_table2_fib.mat','post_fib','fit_fib');
save('posterior_table2_malami.mat','post_malami','fit_malami');

outDir = fullfile(pwd,'figs'); if ~exist(outDir,'dir'); mkdir(outDir); end
writetable([post_fib.summary; post_malami.summary], fullfile(outDir,'posterior_table2_parameters.csv'));

fprintf('\nAcceptance rates: fib = %.3f, Malami = %.3f\n', post_fib.accRate, post_malami.accRate);
fprintf('Saved posterior_table2_fib.mat, posterior_table2_malami.mat and figs/posterior_table2_parameters.csv\n');

% Optional after UQLab is initialized:
% uqlab;
% myInputFib_post = buildInput_carbonation_posterior_table2(fit_fib);
% myInputMalami_post = buildInput_carbonation_posterior_table2(fit_malami);
