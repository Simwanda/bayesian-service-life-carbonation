function post = bayes_update_theta_x(prior, A0, n0, t_obs, y_obs, sigma_eps, mcmc)
% Bayesian updating of theta_x ONLY (Gaussian prior)
%
% y = theta_x * A0 * t^(0.5 - n0) + eps

% ---- prior parameters ----
mu_th = prior.theta.mean;
sig_th = prior.theta.mean * prior.theta.cov;

% ---- log-posterior ----
    function lp = logpost(theta)
        % enforce positivity (physical)
        if theta <= 0
            lp = -Inf;
            return;
        end

        % prior: Normal
        lp_prior = -0.5*log(2*pi*sig_th^2) ...
                   -0.5*((theta - mu_th)/sig_th)^2;

        % likelihood
        mu_y = theta * A0 .* (t_obs .^ (0.5 - n0));
        lp_like = sum( ...
            -0.5*log(2*pi*sigma_eps^2) ...
            -0.5*((y_obs - mu_y)/sigma_eps).^2 );

        lp = lp_prior + lp_like;
    end

% ---- MCMC init ----
theta = mu_th;
lp = logpost(theta);

Nburn = mcmc.Nburn;
Nkeep = mcmc.Nkeep;
thin  = mcmc.thin;
propStd = mcmc.propStd_theta;

Ntot = Nburn + Nkeep*thin;
ThetaKeep = zeros(Nkeep,1);

acc = 0;
k = 0;

% ---- RW-MH ----
for it = 1:Ntot
    theta_prop = theta + propStd*randn;
    lp_prop = logpost(theta_prop);

    if log(rand) < (lp_prop - lp)
        theta = theta_prop;
        lp = lp_prop;
        acc = acc + 1;
    end

    if it > Nburn && mod(it-Nburn,thin)==0
        k = k + 1;
        ThetaKeep(k) = theta;
    end
end

% ---- output ----
post.theta = ThetaKeep;
post.accRate = acc / Ntot;
post.settings = mcmc;
post.sigma_eps = sigma_eps;
post.A0 = A0;
post.n0 = n0;

end
