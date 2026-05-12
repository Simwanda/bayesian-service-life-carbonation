function myInput = buildInput_carbonation_posterior(Aopt, Ccase, fit)
% Posterior-updated input for reliability analyses.
% Aopt : 'fib' or 'malami'
% Ccase: 'caseI' or 'caseII'
% fit  : struct from fit_lognormal_moments with fields:
%        fit.theta (mu_ln, sig_ln), fit.A0, fit.n0

Aopt  = char(Aopt);
Ccase = char(Ccase);

M = [];

% -------------------------------------------------
% Fixed (deterministic) parameters
% -------------------------------------------------
M(1).Name       = 'A';
M(1).Type       = 'Constant';
M(1).Parameters = fit.A0;     % REAL value

M(2).Name       = 'n';
M(2).Type       = 'Constant';
M(2).Parameters = fit.n0;     % REAL value

% -------------------------------------------------
% Posterior-updated model uncertainty factor
% -------------------------------------------------
M(3).Name       = 'theta_x';
M(3).Type       = 'Lognormal';
M(3).Parameters = [fit.theta.mu_ln, fit.theta.sig_ln];

% -------------------------------------------------
% Concrete cover
% -------------------------------------------------
CrDet  = 20.5;
CrMean = 20.5;
CrCov  = 0.30;

if strcmpi(Ccase,'caseII')
    [mu_ln_Cr, sig_ln_Cr] = logn_mu_sigma(CrMean, CrCov);
    M(4).Name       = 'C_r';
    M(4).Type       = 'Lognormal';
    M(4).Parameters = [mu_ln_Cr, sig_ln_Cr];
end

% -------------------------------------------------
% Create UQLab input
% -------------------------------------------------
InputOpts.Marginals = M;
InputOpts.Name = sprintf('InputPosterior_%s_%s', Aopt, Ccase);

myInput = uq_createInput(InputOpts);

% -------------------------------------------------
% Store deterministic case info
% -------------------------------------------------
myInput.Internal.CaseInfo.Ccase = lower(Ccase);
myInput.Internal.CaseInfo.CrDet = CrDet;

end

% =================================================
% Helper
% =================================================
function [mu_ln, sigma_ln] = logn_mu_sigma(m, cov)
sigma_ln = sqrt(log(1 + cov.^2));
mu_ln    = log(m) - 0.5*sigma_ln.^2;
end

