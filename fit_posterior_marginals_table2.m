function fit = fit_posterior_marginals_table2(post)
%FIT_POSTERIOR_MARGINALS_TABLE2 Fit posterior marginals for UQLab input.
% The posterior samples from MCMC are reduced to independent marginal models.
% This is consistent with the paper's stated independence assumption; if needed,
% posterior sample dependence can be kept separately for direct MCS.

fit = struct();
fit.modelName = post.modelName;
fit.updateNames = post.updateNames;
fit.rv = post.priors.rv;

for i = 1:numel(post.priors.rv)
    rv = post.priors.rv(i);
    name = rv.Name;
    x = post.samples.(name)(:);
    f.Name = name;
    f.OriginalType = rv.Type;
    f.Mean = mean(x);
    f.Std = std(x);
    f.COV = safe_cov(f.Mean, f.Std);
    f.Median = median(x);
    f.P2p5 = quantile(x,0.025);
    f.P97p5 = quantile(x,0.975);

    switch lower(rv.Type)
        case 'lognormal'
            [mu_ln,sig_ln] = logn_mom2par(max(f.Mean,realmin), f.COV);
            f.UQType = 'Lognormal';
            f.UQParameters = [mu_ln sig_ln];
        case 'gaussian'
            % Use Gaussian unless physical positivity is violated for theta_xc.
            f.UQType = 'Gaussian';
            f.UQParameters = [f.Mean f.Std];
        case 'uniform'
            f.UQType = 'Uniform';
            f.UQParameters = [f.P2p5 f.P97p5];
        otherwise
            f.UQType = 'Constant';
            f.UQParameters = f.Mean;
    end
    fit.marginals.(name) = f;
end

fit.derived.A.mean = mean(post.aux.A);
fit.derived.A.std = std(post.aux.A);
fit.derived.n.mean = mean(post.aux.n);
fit.derived.n.std = std(post.aux.n);
fit.derived.xc_obs.mean = mean(post.xc_at_obs);
fit.derived.xc_obs.std = std(post.xc_at_obs);
end

function [mu,sig] = logn_mom2par(m,cov)
sig = sqrt(log(1+cov.^2));
mu = log(m)-0.5*sig.^2;
end


function c = safe_cov(mu, sig)
if abs(mu) < 1e-12
    c = NaN;
else
    c = sig/abs(mu);
end
end
