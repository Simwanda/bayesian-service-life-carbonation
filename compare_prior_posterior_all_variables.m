function compare_prior_posterior_all_variables()
%COMPARE_PRIOR_POSTERIOR_ALL_VARIABLES
% Compare prior and posterior distributions for all Table 2 variables and
% compute posterior correlation matrices for fib Bulletin 34 and Malami.
%
% Usage:
%   1) Run the Bayesian update first:
%        run_bayes_update_table2_fib_malami
%   2) Then run:
%        compare_prior_posterior_all_variables
%
% Outputs are saved in:
%   figs/prior_posterior_<model>_<variable>.png
%   figs/prior_posterior_<model>_ALL.pdf
%   figs/correlation_matrix_<model>_pearson.csv
%   figs/correlation_matrix_<model>_spearman.csv
%   figs/correlation_matrix_<model>_pearson.png
%   figs/correlation_matrix_<model>_spearman.png
%   figs/posterior_samples_<model>_with_aux.csv
%
% Notes:
%   - Prior PDFs are taken from table2_priors_carbonation.m.
%   - Posterior samples are taken from posterior_table2_fib.mat and
%     posterior_table2_malami.mat.
%   - Variables that were not updated are usually constant posterior samples;
%     they are plotted, but excluded from correlation matrices because their
%     variance is zero.

clc;
outDir = fullfile(pwd, 'figs');
if ~exist(outDir, 'dir'); mkdir(outDir); end

% Load posterior results. If files do not exist, tell the user what to run.
if ~exist('posterior_table2_fib.mat','file') || ~exist('posterior_table2_malami.mat','file')
    error(['Posterior .mat files were not found in the current folder. ', ...
           'Run run_bayes_update_table2_fib_malami first, then run this script.']);
end

Sfib = load('posterior_table2_fib.mat');
Smal = load('posterior_table2_malami.mat');

post_fib = Sfib.post_fib;
post_malami = Smal.post_malami;

make_prior_posterior_plots(post_fib, 'fib', outDir);
make_prior_posterior_plots(post_malami, 'malami', outDir);

make_correlation_outputs(post_fib, 'fib', outDir);
make_correlation_outputs(post_malami, 'malami', outDir);

fprintf('\nDone. Outputs saved in: %s\n', outDir);
end

% -------------------------------------------------------------------------
function make_prior_posterior_plots(post, modelName, outDir)
priors = post.priors;
vars = priors.rv;

pdfFile = fullfile(outDir, sprintf('prior_posterior_%s_ALL.pdf', modelName));
if exist(pdfFile,'file'); delete(pdfFile); end

for i = 1:numel(vars)
    rv = vars(i);
    nm = rv.Name;
    if ~isfield(post.samples, nm)
        warning('No posterior samples found for %s/%s. Skipping.', modelName, nm);
        continue;
    end

    xpost = post.samples.(nm)(:);
    xpost = xpost(isfinite(xpost));
    if isempty(xpost)
        warning('Posterior samples for %s/%s are empty. Skipping.', modelName, nm);
        continue;
    end

    % Define plotting range robustly using prior range and posterior quantiles.
    [xgrid, priorPDF] = prior_pdf_grid(rv, xpost);

    f = figure('Color','w', 'Visible','off');
    histogram(xpost, 45, 'Normalization','pdf', 'FaceAlpha',0.35, 'EdgeAlpha',0.15);
    hold on;
    plot(xgrid, priorPDF, 'LineWidth',2.0);

    % Posterior kernel or fitted smooth curve, if available.
    [xk, fk] = posterior_kernel_pdf(xpost);
    if ~isempty(xk)
        plot(xk, fk, '--', 'LineWidth',2.0);
        leg = {'Posterior samples','Prior PDF','Posterior KDE'};
    else
        leg = {'Posterior samples','Prior PDF'};
    end

    xlabel(label_with_unit(nm, rv.Unit), 'Interpreter','none');
    ylabel('Density');
    title(sprintf('%s: prior vs posterior for %s', model_label(modelName), nm), 'Interpreter','none');
    legend(leg, 'Location','best');
    grid on; box on;

    pngFile = fullfile(outDir, sprintf('prior_posterior_%s_%s.png', modelName, safe_name(nm)));
    exportgraphics(f, pngFile, 'Resolution', 300);
    exportgraphics(f, pdfFile, 'Append', true);
    close(f);

    % Save posterior sample summary for this variable as a small text report.
    fprintf('%s %-12s posterior mean = %.6g, std = %.6g, CoV = %.4g\n', ...
        modelName, nm, mean(xpost), std(xpost), std(xpost)/abs(mean(xpost)));
end
end

% -------------------------------------------------------------------------
function make_correlation_outputs(post, modelName, outDir)
% Use updated variables plus A and n auxiliary variables. Exclude constants.
baseNames = post.updateNames(:)';

samples = post.samples;
N = numel(samples.(baseNames{1}));

% Add auxiliary variables A and n because they are useful for explaining
% carbonation-depth updating and reliability calculations.
aux = struct();
if isfield(post,'aux')
    if isfield(post.aux,'A'); aux.A = post.aux.A(:); end
    if isfield(post.aux,'n'); aux.n = post.aux.n(:); end
end

allNames = baseNames;
X = zeros(N, numel(baseNames));
for j = 1:numel(baseNames)
    X(:,j) = samples.(baseNames{j})(:);
end

auxNames = fieldnames(aux);
for k = 1:numel(auxNames)
    v = aux.(auxNames{k});
    if numel(v)==N
        allNames{end+1} = auxNames{k}; %#ok<AGROW>
        X(:,end+1) = v(:); %#ok<AGROW>
    end
end

% Remove variables with zero or near-zero posterior variance.
sd = std(X,0,1,'omitnan');
keep = isfinite(sd) & sd > 1e-12 .* max(1, abs(mean(X,1,'omitnan')));
X = X(:,keep);
allNames = allNames(keep);

% Remove rows with any NaNs.
rowKeep = all(isfinite(X),2);
X = X(rowKeep,:);

Rpear = corrcoef(X);
Rspear = corrcoef(rank_columns(X));

Tpear = array2table(Rpear, 'VariableNames', matlab.lang.makeValidName(allNames), 'RowNames', allNames);
Tspear = array2table(Rspear, 'VariableNames', matlab.lang.makeValidName(allNames), 'RowNames', allNames);

writetable(Tpear, fullfile(outDir, sprintf('correlation_matrix_%s_pearson.csv', modelName)), 'WriteRowNames', true);
writetable(Tspear, fullfile(outDir, sprintf('correlation_matrix_%s_spearman.csv', modelName)), 'WriteRowNames', true);

% Also export posterior samples used for the correlations.
Tout = array2table(X, 'VariableNames', matlab.lang.makeValidName(allNames));
writetable(Tout, fullfile(outDir, sprintf('posterior_samples_%s_with_aux.csv', modelName)));

plot_corr_heatmap(Rpear, allNames, sprintf('%s posterior Pearson correlation', model_label(modelName)), ...
    fullfile(outDir, sprintf('correlation_matrix_%s_pearson.png', modelName)));
plot_corr_heatmap(Rspear, allNames, sprintf('%s posterior Spearman correlation', model_label(modelName)), ...
    fullfile(outDir, sprintf('correlation_matrix_%s_spearman.png', modelName)));

fprintf('\n%s correlation matrices saved. Variables included:\n', modelName);
fprintf('  %s\n', strjoin(allNames, ', '));
end

% -------------------------------------------------------------------------
function [xgrid, f] = prior_pdf_grid(rv, xpost)
qPost = quantile_local(xpost, [0.001 0.999]);

switch lower(rv.Type)
    case 'gaussian'
        lo = min(qPost(1), rv.Mean - 4*rv.Std);
        hi = max(qPost(2), rv.Mean + 4*rv.Std);
        if lo == hi; lo = rv.Mean - 1; hi = rv.Mean + 1; end
        xgrid = linspace(lo, hi, 500);
        f = norm_pdf_local(xgrid, rv.Mean, rv.Std);

    case 'lognormal'
        [mu, sig] = logn_mom2par(rv.Mean, rv.COV);
        priorLo = exp(mu - 4*sig);
        priorHi = exp(mu + 4*sig);
        lo = max(0, min(qPost(1), priorLo));
        hi = max(qPost(2), priorHi);
        if lo <= 0; lo = min(xpost(xpost>0)); end
        if isempty(lo) || ~isfinite(lo) || lo <= 0; lo = eps; end
        if hi <= lo; hi = lo*10; end
        xgrid = linspace(lo, hi, 500);
        f = logn_pdf_local(xgrid, mu, sig);

    case 'uniform'
        span = rv.Upper - rv.Lower;
        lo = min(qPost(1), rv.Lower - 0.10*span);
        hi = max(qPost(2), rv.Upper + 0.10*span);
        if hi <= lo; hi = lo + 1; end
        xgrid = linspace(lo, hi, 500);
        f = zeros(size(xgrid));
        idx = xgrid >= rv.Lower & xgrid <= rv.Upper;
        f(idx) = 1/span;

    otherwise
        error('Unsupported prior distribution: %s', rv.Type);
end
end

% -------------------------------------------------------------------------
function [xk, fk] = posterior_kernel_pdf(x)
% Kernel smoothing without requiring Statistics Toolbox ksdensity.
x = x(:);
x = x(isfinite(x));
if numel(unique(x)) < 5
    xk = []; fk = [];
    return;
end
q = quantile_local(x, [0.001 0.999]);
lo = q(1); hi = q(2);
if hi <= lo
    lo = min(x); hi = max(x);
end
if hi <= lo
    xk = []; fk = [];
    return;
end
pad = 0.1*(hi-lo);
xk = linspace(lo-pad, hi+pad, 400);

sx = std(x);
n = numel(x);
h = 1.06*sx*n^(-1/5); % Silverman's rule
if ~isfinite(h) || h <= 0
    xk = []; fk = [];
    return;
end

% Vectorized Gaussian KDE. For very large samples, thin for speed.
maxN = 8000;
if n > maxN
    idx = round(linspace(1,n,maxN));
    x = sort(x);
    x = x(idx);
    n = numel(x);
end
Z = (xk(:) - x(:)') ./ h;
fk = mean(exp(-0.5*Z.^2),2) ./ (h*sqrt(2*pi));
fk = fk(:)';
end

% -------------------------------------------------------------------------
function plot_corr_heatmap(R, names, ttl, outFile)
f = figure('Color','w', 'Visible','off');
imagesc(R, [-1 1]); axis equal tight;
colormap(parula); colorbar;
set(gca, 'XTick',1:numel(names), 'XTickLabel',names, ...
         'YTick',1:numel(names), 'YTickLabel',names, ...
         'TickLabelInterpreter','none', 'XTickLabelRotation',45);
title(ttl, 'Interpreter','none');

% Add numeric correlation labels.
for i = 1:size(R,1)
    for j = 1:size(R,2)
        text(j, i, sprintf('%.2f', R(i,j)), 'HorizontalAlignment','center', ...
            'FontSize', 8, 'Color','k');
    end
end

exportgraphics(f, outFile, 'Resolution', 300);
close(f);
end

% -------------------------------------------------------------------------
function R = rank_columns(X)
R = zeros(size(X));
for j = 1:size(X,2)
    R(:,j) = tied_rank_local(X(:,j));
end
end

function r = tied_rank_local(x)
% Simple tied ranks, independent of Statistics Toolbox.
[xs, idx] = sort(x(:));
r = zeros(size(xs));
n = numel(xs);
i = 1;
while i <= n
    j = i;
    while j < n && xs(j+1) == xs(i)
        j = j + 1;
    end
    r(i:j) = (i+j)/2;
    i = j + 1;
end
out = zeros(size(r));
out(idx) = r;
r = out;
end

% -------------------------------------------------------------------------
function q = quantile_local(x, p)
x = sort(x(:));
x = x(isfinite(x));
n = numel(x);
q = nan(size(p));
for i = 1:numel(p)
    if n == 0; q(i) = NaN; continue; end
    pos = 1 + (n-1)*p(i);
    lo = floor(pos); hi = ceil(pos);
    if lo == hi
        q(i) = x(lo);
    else
        q(i) = x(lo) + (pos-lo)*(x(hi)-x(lo));
    end
end
end

function f = norm_pdf_local(x, mu, sig)
sig = max(sig, realmin);
f = exp(-0.5*((x-mu)./sig).^2) ./ (sig*sqrt(2*pi));
end

function f = logn_pdf_local(x, mu, sig)
f = zeros(size(x));
idx = x > 0;
f(idx) = exp(-0.5*((log(x(idx))-mu)./sig).^2) ./ (x(idx).*sig*sqrt(2*pi));
end

function [mu,sig] = logn_mom2par(m,cov)
sig = sqrt(log(1+cov.^2));
mu = log(m)-0.5*sig.^2;
end

function s = label_with_unit(name, unit)
if isempty(unit) || strcmp(unit,'-')
    s = name;
else
    s = sprintf('%s [%s]', name, unit);
end
end

function s = model_label(modelName)
switch lower(modelName)
    case 'fib'
        s = 'fib Bulletin 34';
    case 'malami'
        s = 'Malami';
    otherwise
        s = modelName;
end
end

function s = safe_name(s)
s = regexprep(s, '[^A-Za-z0-9_]', '_');
end
