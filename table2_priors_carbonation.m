function priors = table2_priors_carbonation(modelName)
%TABLE2_PRIORS_CARBONATION Prior probabilistic models for CT100/60.
%   Priors follow Table 2 of the manuscript. The returned struct is used by
%   the MCMC routines and can also be converted to UQLab marginals.
%
%   modelName = 'fib' or 'malami'.
%
%   Distribution conventions used in this code:
%     Gaussian  : parameters are mean and standard deviation.
%     Lognormal : parameters stored as mean and CoV; internally transformed.
%     Uniform   : parameters are lower and upper bounds, derived from mean/CoV
%                 when the paper gives mean and CoV only.
%     Constant  : deterministic value, not updated by MCMC.

modelName = lower(char(modelName));

% -------------------------
% Constants common to models
% -------------------------
priors.modelName = modelName;
priors.t_ref_days = 28;
priors.t_ref_years = 28/365.25;
priors.A_prior_mean_fib_mm_y = 2.21;    % revised paper value, mm/year^(0.5-n)
priors.A_prior_mean_malami_mm_y = 2.17; % revised paper value, mm/year^(0.5-n)
priors.n_prior_mean = 0.11;

% -------------------------
% Table 2 stochastic variables
% -------------------------
P = empty_rv_array();

switch modelName
    case 'fib'
        P(end+1) = rv('Racc_inv',   'Lognormal', 9.80e-11, 0.42, '(m2/s)/(kg/m3)', 'Accelerated inverse resistance R_ACC,0^{-1}');
        P(end+1) = rv('kt',         'Lognormal', 1.25,     0.22, '-', 'Regression parameter for ACC-to-natural transfer');
        P(end+1) = rv('eps_t',      'Gaussian',  1.00e-11, 0.152, '(m2/s)/(kg/m3)', 'Error term of ACC test method; CoV interpreted relative to mean');
        P(end+1) = rv_uniform_cov('tc', 2.0, 0.289, 'days', 'Curing time, equivalent to Uniform[1,3] days');
        P(end+1) = rv('bc',         'Gaussian', -0.567,    0.042, '-', 'Rate constant exponent');
        P(end+1) = rv('CO2ppm',     'Gaussian', 392,       0.009, 'ppm', 'Long-term CO2 concentration');
        P(end+1) = rv('RH',         'Gaussian', 75.6,      0.035, '%', 'Relative humidity');
        P(end+1) = rv('tw',         'Gaussian', 0.156,     0.115, '-', 'Time of wetness');
        P(end+1) = rv('bw',         'Gaussian', 0.446,     0.365, '-', 'Time-of-wetness exponent');
        P(end+1) = rv_uniform_cov('pSR', 0.146, 0.40, '-', 'Probability of driving rain');
        P(end+1) = rv('theta_xc',   'Gaussian', 1.0,       0.20, '-', 'Multiplicative carbonation-depth model uncertainty');

    case 'malami'
        P(end+1) = rv('W_C',        'Gaussian', 0.55,      0.10, '-', 'Water-to-cement ratio');
        P(end+1) = rv('tc',         'Uniform',  1.0,       3.0, 'days', 'Curing time, Uniform[1,3] days');
        P(end+1) = rv('bc',         'Gaussian', -0.567,    0.042, '-', 'Rate constant exponent');
        P(end+1) = rv('CO2ppm',     'Gaussian', 392,       0.009, 'ppm', 'Long-term CO2 concentration');
        P(end+1) = rv('RH',         'Gaussian', 75.6,      0.035, '%', 'Relative humidity');
        P(end+1) = rv('tw',         'Gaussian', 0.156,     0.115, '-', 'Time of wetness');
        P(end+1) = rv('bw',         'Gaussian', 0.446,     0.365, '-', 'Time-of-wetness exponent');
        P(end+1) = rv_uniform_cov('pSR', 0.146, 0.40, '-', 'Probability of driving rain');
        P(end+1) = rv('T',          'Gaussian', 283.25,    0.002, 'K', 'Design temperature');
        P(end+1) = rv('theta_xc',   'Gaussian', 1.0,       0.20, '-', 'Multiplicative carbonation-depth model uncertainty');

        % Deterministic Table 2 material-composition inputs for Malami.
        priors.det.C = 375;        % kg/m3
        priors.det.P = 0;          % kg/m3
        priors.det.k = 1;          % -
        priors.det.phi_cl = 0.975; % -
        priors.det.CaO = 0.68;     % -
        priors.det.kurb = 1.0;     % -
        priors.det.rhoC = 3120;    % kg/m3
        priors.det.rhoW = 1000;    % kg/m3

    otherwise
        error('Unknown modelName: use ''fib'' or ''malami''.');
end

priors.rv = P;
priors.names = {P.Name};

% Store prior mean vector as a struct for reference calculations.
priors.meanStruct = prior_mean_struct(priors);

end

function P = empty_rv_array()
% Create a 0-by-1 struct array with exactly the same fields as rv().
% This avoids MATLAB versions that reject P(end+1)=struct(...) when P has no fields.
t = struct('Name', '', ...
           'Type', '', ...
           'Unit', '', ...
           'Comment', '', ...
           'Mean', NaN, ...
           'COV', NaN, ...
           'Std', NaN, ...
           'Lower', NaN, ...
           'Upper', NaN, ...
           'Update', true);
P = repmat(t,0,1);
end

function s = rv(name, type, mean_or_a, cov_or_b, unit, comment)
% For Gaussian and Lognormal, cov_or_b is CoV (not standard deviation).
% For Uniform, mean_or_a and cov_or_b are lower and upper bounds.
%
% IMPORTANT MATLAB COMPATIBILITY NOTE:
% All random-variable structs must contain the same fields before assignment
% with P(end+1)=..., otherwise MATLAB throws:
%   "Subscripted assignment between dissimilar structures."
% Therefore Lower/Upper/Std are initialized for every distribution type.
s = struct('Name', name, ...
           'Type', type, ...
           'Unit', unit, ...
           'Comment', comment, ...
           'Mean', NaN, ...
           'COV', NaN, ...
           'Std', NaN, ...
           'Lower', NaN, ...
           'Upper', NaN, ...
           'Update', true);

if strcmpi(type,'Uniform')
    s.Lower = mean_or_a;
    s.Upper = cov_or_b;
    s.Mean = 0.5*(s.Lower+s.Upper);
    s.Std = (s.Upper-s.Lower)/sqrt(12);
    s.COV = s.Std/abs(s.Mean);
else
    s.Mean = mean_or_a;
    s.COV = cov_or_b;
    s.Std = abs(mean_or_a)*cov_or_b;
end
end

function s = rv_uniform_cov(name, meanVal, covVal, unit, comment)
halfWidth = sqrt(3) * meanVal * covVal;
s = rv(name, 'Uniform', meanVal-halfWidth, meanVal+halfWidth, unit, comment);
s.Mean = meanVal;
s.COV = covVal;
end

function ms = prior_mean_struct(priors)
ms = struct();
for i = 1:numel(priors.rv)
    ms.(priors.rv(i).Name) = priors.rv(i).Mean;
end
if isfield(priors,'det')
    fn = fieldnames(priors.det);
    for k = 1:numel(fn)
        ms.(fn{k}) = priors.det.(fn{k});
    end
end
end
