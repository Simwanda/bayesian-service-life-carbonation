function S = posterior_summary(post)

theta = post.theta(:);

S = table();
S.Parameter = {'theta_x'};
S.Mean    = mean(theta);
S.Std     = std(theta);
S.COV     = S.Std / S.Mean;
S.Median  = median(theta);
S.P2p5    = quantile(theta,0.025);
S.P97p5   = quantile(theta,0.975);
S.AccRate = post.accRate;

end
