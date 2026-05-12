function post = bayes_update_A_n_theta(priorA, priorCommon, t_obs, y_obs, sigma_eps, mcmc)
% Bayesian updating of A, n, theta (all positive) using RW-Metropolis in log-space.
% priors: A lognormal from mean/cov; n lognormal from mean/cov; theta lognormal from mean/cov.

% --- convert priors (mean,cov) -> lognormal (mu_ln, sig_ln) ---
[muA, sigA] = logn_mu_sigma(priorA.A.mean, priorA.A.cov);
[mun, sign] = logn_mu_sigma(priorCommon.n.mean, priorCommon.n.cov);
[muth, sigth] = logn_mu_sigma(priorCommon.theta.mean, priorCommon.theta.cov);

% --- log-posterior in transformed variables: z = [logA logn logtheta] ---
    function lp = logpost(z)
        A = exp(z(1));
        n = exp(z(2));
        th = exp(z(3));

        % prior densities in original variables (lognormal) + Jacobian handled by using log(A) etc?
        % easiest: write log prior in terms of log(A), log(n), log(theta) directly:
        % if X ~ Lognormal(mu,sig), then log X ~ N(mu,sig)
        lp_prior = lognormpdf_log(A, muA, sigA) + ...
                   lognormpdf_log(n, mun, sign) + ...
                   lognormpdf_log(th, muth, sigth);

        % likelihood: y_i ~ N(th*A*t^(0.5-n), sigma_eps^2)
        mu_y = th .* A .* (t_obs .^ (0.5 - n));
        lp_like = sum(-0.5*log(2*pi*sigma_eps^2) - 0.5*((y_obs - mu_y)/sigma_eps).^2);

        lp = lp_prior + lp_like;
    end

% --- initial state at prior means ---
z = [log(priorA.A.mean), log(priorCommon.n.mean), log(priorCommon.theta.mean)];
lp = logpost(z);

Nburn = mcmc.Nburn;
Nkeep = mcmc.Nkeep;
thin  = mcmc.thin;
propStd = mcmc.propStd(:)';

% --- burn-in + sampling ---
Ntot = Nburn + Nkeep*thin;
Zkeep = zeros(Nkeep,3);

acc = 0;
kkeep = 0;

for it = 1:Ntot
    z_prop = z + propStd .* randn(1,3);
    lp_prop = logpost(z_prop);

    if log(rand) < (lp_prop - lp)
        z = z_prop;
        lp = lp_prop;
        acc = acc + 1;
    end

    if it > Nburn && mod(it-Nburn, thin)==0
        kkeep = kkeep + 1;
        Zkeep(kkeep,:) = z;
    end
end

accRate = acc / Ntot;

post.A     = exp(Zkeep(:,1));
post.n     = exp(Zkeep(:,2));
post.theta = exp(Zkeep(:,3));
post.accRate = accRate;
post.settings = mcmc;
post.sigma_eps = sigma_eps;

end

% ---- helpers ----
function [mu_ln, sigma_ln] = logn_mu_sigma(m, cov)
sigma_ln = sqrt(log(1 + cov.^2));
mu_ln    = log(m) - 0.5*sigma_ln.^2;
end

function lp = lognormpdf_log(x, mu_ln, sig_ln)
% log pdf of lognormal evaluated at x > 0
lp = -log(x) - log(sig_ln*sqrt(2*pi)) - 0.5*((log(x)-mu_ln)/sig_ln).^2;
end
