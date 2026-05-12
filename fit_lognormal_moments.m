function fit = fit_lognormal_moments(post)
% Fits lognormal moment parameters for theta_x only
% A and n are deterministic and stored separately

fit = struct();

% -------------------------------------------------
% theta_x (lognormal)
% -------------------------------------------------
m   = mean(post.theta);
cov = std(post.theta) / m;

fit.theta = moment2logn(m, cov);

% -------------------------------------------------
% Store deterministic parameters (for traceability)
% -------------------------------------------------
fit.A0 = post.A0;
fit.n0 = post.n0;

end

% =================================================
% Helper
% =================================================
function p = moment2logn(m, cov)
sigma_ln = sqrt(log(1 + cov^2));
mu_ln    = log(m) - 0.5*sigma_ln^2;

p.mu_ln = mu_ln;
p.sig_ln = sigma_ln;
p.mean  = m;
p.cov   = cov;
end
