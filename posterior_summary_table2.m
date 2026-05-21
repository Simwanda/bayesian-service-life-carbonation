function T = posterior_summary_table2(post)
%POSTERIOR_SUMMARY_TABLE2 Summary table for all updated Table 2 variables.

names = post.updateNames(:);
Model = repmat({post.modelName}, numel(names)+3, 1);
Parameter = [names; {'A_derived'; 'n_derived'; 'xc_at_obs'}];
Mean = nan(numel(Parameter),1);
Std = nan(numel(Parameter),1);
COV = nan(numel(Parameter),1);
Median = nan(numel(Parameter),1);
P2p5 = nan(numel(Parameter),1);
P97p5 = nan(numel(Parameter),1);
AccRate = nan(numel(Parameter),1);

for i = 1:numel(names)
    x = post.samples.(names{i})(:);
    Mean(i) = mean(x); Std(i) = std(x); COV(i) = safe_cov(Mean(i), Std(i));
    Median(i) = median(x); P2p5(i) = quantile(x,0.025); P97p5(i)=quantile(x,0.975);
    AccRate(i) = post.accRate;
end
extra = {post.aux.A(:), post.aux.n(:), post.xc_at_obs(:)};
for j = 1:3
    i = numel(names)+j;
    x = extra{j};
    Mean(i)=mean(x); Std(i)=std(x); COV(i)=safe_cov(Mean(i), Std(i));
    Median(i)=median(x); P2p5(i)=quantile(x,0.025); P97p5(i)=quantile(x,0.975);
end
T = table(Model, Parameter, Mean, Std, COV, Median, P2p5, P97p5, AccRate);
end


function c = safe_cov(mu, sig)
if abs(mu) < 1e-12
    c = NaN;
else
    c = sig/abs(mu);
end
end
