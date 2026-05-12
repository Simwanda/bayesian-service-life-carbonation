%EXPORT_POSTERIOR_TABLE2 Export posterior summaries for Table 2 Bayesian update.
clear; clc;
outDir = fullfile(pwd,'figs'); if ~exist(outDir,'dir'); mkdir(outDir); end
Sfib = load('posterior_table2_fib.mat');
Smal = load('posterior_table2_malami.mat');
T = [Sfib.post_fib.summary; Smal.post_malami.summary];
disp(T);
writetable(T, fullfile(outDir,'posterior_table2_parameters.csv'));
fprintf('Saved %s\n', fullfile(outDir,'posterior_table2_parameters.csv'));
