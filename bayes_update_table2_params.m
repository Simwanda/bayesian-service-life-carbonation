function post = bayes_update_table2_params(modelName, obs, mcmc, updateNames)
%BAYES_UPDATE_TABLE2_PARAMS Bayesian updating of Table 2 basic variables.
%   post = bayes_update_table2_params(modelName, obs, mcmc, updateNames)
%
%   Responds to reviewer R2.4 by updating probabilistic model parameters
%   directly rather than updating theta_xc alone. The prior distributions are
%   those from Table 2. The likelihood can use either raw carbonation-depth
%   observations or the inspection mean/std summary.
%
%   obs fields:
%      t_years  : scalar or vector of inspection ages [years]
%      y_mm     : raw observations [mm], optional
%      y_mean   : mean carbonation depth [mm], used if y_mm is not given
%      y_std    : observed standard deviation [mm]
%      N_eff    : effective number of independent observations for the mean
%
%   mcmc fields:
%      Nburn, Nkeep, thin, propStd, rngSeed
%
%   updateNames optional cell array. Default: all stochastic Table 2 variables.

modelName = lower(char(modelName));
priors = table2_priors_carbonation(modelName);

if nargin < 4 || isempty(updateNames)
    updateNames = priors.names;
end
updateNames = updateNames(:)';

if ~isfield(mcmc,'rngSeed'); mcmc.rngSeed = 20260512; end
rng(mcmc.rngSeed,'twister');
if ~isfield(mcmc,'Nburn'); mcmc.Nburn = 20000; end
if ~isfield(mcmc,'Nkeep'); mcmc.Nkeep = 40000; end
if ~isfield(mcmc,'thin');  mcmc.thin  = 5; end
if ~isfield(mcmc,'propStd') || isempty(mcmc.propStd)
    mcmc.propStd = default_prop_std(priors, updateNames);
end

% Initial vector in unconstrained space
z = zeros(1,numel(updateNames));
for j = 1:numel(updateNames)
    rvj = get_rv(priors, updateNames{j});
    z(j) = x_to_z(rvj.Mean, rvj);
end

lp = logpost(z);
Ntot = mcmc.Nburn + mcmc.Nkeep*mcmc.thin;
Zkeep = zeros(mcmc.Nkeep,numel(updateNames));
acc = 0; kkeep = 0;

for it = 1:Ntot
    zprop = z + mcmc.propStd(:)'.*randn(size(z));
    lpprop = logpost(zprop);
    if log(rand) < (lpprop-lp)
        z = zprop; lp = lpprop; acc = acc + 1;
    end
    if it > mcmc.Nburn && mod(it-mcmc.Nburn,mcmc.thin)==0
        kkeep = kkeep + 1;
        Zkeep(kkeep,:) = z;
    end
end

% Convert saved samples to physical variables
samples = struct();
for j = 1:numel(updateNames)
    rvj = get_rv(priors, updateNames{j});
    samples.(updateNames{j}) = z_to_x(Zkeep(:,j), rvj);
end
% Add deterministic/not-updated variables at prior mean for traceability
for j = 1:numel(priors.names)
    nm = priors.names{j};
    if ~isfield(samples,nm)
        samples.(nm) = repmat(get_rv(priors,nm).Mean,mcmc.Nkeep,1);
    end
end
if isfield(priors,'det')
    fn = fieldnames(priors.det);
    for k = 1:numel(fn)
        samples.(fn{k}) = repmat(priors.det.(fn{k}),mcmc.Nkeep,1);
    end
end

[xc42, aux] = carbonation_xc_table2(modelName, samples, obs.t_years(1), priors);

post.modelName = modelName;
post.priors = priors;
post.updateNames = updateNames;
post.samples = samples;
post.aux.A = aux.A;
post.aux.n = aux.n;
post.aux.theta_xc = aux.theta_xc;
post.xc_at_obs = xc42;
post.accRate = acc/Ntot;
post.settings = mcmc;
post.obs = obs;
post.summary = posterior_summary_table2(post);

    function lpv = logpost(zz)
        X = priors.meanStruct;
        lpv = 0;
        for jj = 1:numel(updateNames)
            rv = get_rv(priors, updateNames{jj});
            x = z_to_x(zz(jj), rv);
            lpv = lpv + log_prior_z(zz(jj), rv); % prior + transform Jacobian
            X.(updateNames{jj}) = x;
        end
        if ~isfinite(lpv); return; end
        [muPred,~] = carbonation_xc_table2(modelName, X, obs.t_years, priors);
        muPred = muPred(:);
        lpv = lpv + log_likelihood_obs(muPred, obs);
    end
end

% -------------------------------------------------------------------------
function lp = log_likelihood_obs(muPred, obs)
if isfield(obs,'y_mm') && ~isempty(obs.y_mm)
    y = obs.y_mm(:);
    if numel(muPred)==1; muPred = repmat(muPred,numel(y),1); end
    if ~isfield(obs,'sigma_eps') || isempty(obs.sigma_eps)
        sigma = max(std(y), 1e-6);
    else
        sigma = obs.sigma_eps;
    end
    lp = sum(log_norm_pdf(y, muPred, sigma));
else
    % Summary-data likelihood: ybar | x ~ N(x, sqrt(y_std^2/N_eff + sigma_eps^2))
    if ~isfield(obs,'N_eff') || isempty(obs.N_eff); obs.N_eff = 1; end
    if ~isfield(obs,'sigma_eps') || isempty(obs.sigma_eps); obs.sigma_eps = 0; end
    sigmaMean = sqrt((obs.y_std.^2)./obs.N_eff + obs.sigma_eps.^2);
    lp = sum(log_norm_pdf(obs.y_mean(:), muPred, sigmaMean));
end
end

function v = default_prop_std(priors, updateNames)
v = zeros(1,numel(updateNames));
for i = 1:numel(updateNames)
    rv = get_rv(priors, updateNames{i});
    switch lower(rv.Type)
        case 'lognormal'; v(i)=0.08;
        case 'gaussian';  v(i)=0.12;
        case 'uniform';   v(i)=0.10;
        otherwise;        v(i)=0.10;
    end
end
end

function rv = get_rv(priors, name)
idx = find(strcmp({priors.rv.Name}, name),1);
if isempty(idx); error('Unknown random variable: %s', name); end
rv = priors.rv(idx);
end

function x = z_to_x(z, rv)
switch lower(rv.Type)
    case 'lognormal'
        [mu,sig] = logn_mom2par(rv.Mean, rv.COV);
        x = exp(mu + sig.*z); % z ~ N(0,1)
    case 'gaussian'
        x = rv.Mean + rv.Std.*z;
    case 'uniform'
        u = 1./(1+exp(-z));
        x = rv.Lower + (rv.Upper-rv.Lower).*u;
    otherwise
        error('Unsupported distribution: %s', rv.Type);
end
end

function z = x_to_z(x, rv)
switch lower(rv.Type)
    case 'lognormal'
        [mu,sig] = logn_mom2par(rv.Mean, rv.COV);
        z = (log(x)-mu)/sig;
    case 'gaussian'
        z = (x-rv.Mean)/rv.Std;
    case 'uniform'
        u = (x-rv.Lower)/(rv.Upper-rv.Lower);
        u = min(max(u,1e-9),1-1e-9);
        z = log(u/(1-u));
    otherwise
        error('Unsupported distribution: %s', rv.Type);
end
end

function lp = log_prior_z(z, rv)
% Prior density in the transformed z-space. For Gaussian/lognormal z is
% standard normal. For Uniform x, z uses logistic transform and includes dx/dz.
switch lower(rv.Type)
    case {'lognormal','gaussian'}
        lp = log_norm_pdf(z,0,1);
    case 'uniform'
        u = 1/(1+exp(-z));
        dx_dz = (rv.Upper-rv.Lower)*u*(1-u);
        if dx_dz <= 0
            lp = -Inf;
        else
            lp = -log(rv.Upper-rv.Lower) + log(dx_dz);
        end
    otherwise
        lp = -Inf;
end
end

function [mu,sig] = logn_mom2par(m,cov)
sig = sqrt(log(1+cov.^2));
mu = log(m)-0.5*sig.^2;
end

function lp = log_norm_pdf(x,mu,sig)
sig = max(sig, realmin);
lp = -0.5*log(2*pi*sig.^2) -0.5*((x-mu)./sig).^2;
end
