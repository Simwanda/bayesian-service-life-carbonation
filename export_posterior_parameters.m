% ==========================================================
% export_posterior_parameters.m
% Outputs posterior summaries for:
%   - theta_x (random)
%   - A0, n0   (deterministic)
% Saves CSV: figs/posterior_parameters.csv
% ==========================================================
clear; clc;

outDir = fullfile(pwd,'figs');
if ~exist(outDir,'dir'); mkdir(outDir); end
outCSV = fullfile(outDir,'posterior_parameters.csv');

% ---- load posterior files ----
Sfib = load('posterior_samples_fib.mat');
Smal = load('posterior_samples_malami.mat');

post_fib = Sfib.post_fib;

if isfield(Smal,'post_mal')
    post_mal = Smal.post_mal;
elseif isfield(Smal,'post_malami')
    post_mal = Smal.post_malami;
else
    error('Cannot find post_mal or post_malami inside posterior_samples_malami.mat');
end

% ---- build summary tables ----
Tfib = make_summary_table(post_fib, 'fib');
Tmal = make_summary_table(post_mal, 'Malami');

T = [Tfib; Tmal];

% ---- display ----
disp(T);

% ---- save ----
writetable(T, outCSV);
fprintf('\nSaved: %s\n', outCSV);

% ==========================================================
% Helper: summary table
% ==========================================================
function T = make_summary_table(post, modelName)

% ---------- theta_x (random) ----------
x = post.theta(:);

Model  = {modelName; modelName; modelName};
Param  = {'theta_x'; 'A'; 'n'};

Mean   = [
    mean(x);
    post.A0;
    post.n0
];

Std    = [
    std(x);
    NaN;
    NaN
];

COV    = [
    Std(1)/Mean(1);
    NaN;
    NaN
];

Median = [
    median(x);
    post.A0;
    post.n0
];

P2p5   = [
    quantile(x,0.025);
    NaN;
    NaN
];

P97p5  = [
    quantile(x,0.975);
    NaN;
    NaN
];

AccRate = [
    post.accRate;
    NaN;
    NaN
];

T = table(Model, Param, Mean, Std, COV, Median, P2p5, P97p5, AccRate);

end

