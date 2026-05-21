%PRODUCE_FIGURE8_POSTERIOR_CARBONATION
% Reproduces Fig. 8: posterior predictive carbonation depth x_c(t)
% for fib Bulletin 34 and Malami models after Table 2 Bayesian updating.
%
% Required workflow:
%   1) Run: run_bayes_update_table2_fib_malami
%   2) Run: produce_figure8_posterior_carbonation
%
% If posterior_table2_fib.mat or posterior_table2_malami.mat is missing,
% this script automatically runs run_bayes_update_table2_fib_malami first.
%
% Outputs:
%   figs/Figure8_posterior_carbonation_depth.pdf
%   figs/Figure8_posterior_carbonation_depth.png
%   figs/Figure8_posterior_carbonation_depth_data.csv

clear; clc; close all;

outDir = fullfile(pwd,'figs');
if ~exist(outDir,'dir'); mkdir(outDir); end

% -------------------------------------------------------------------------
% Ensure posterior samples exist
% -------------------------------------------------------------------------
if ~exist('posterior_table2_fib.mat','file') || ~exist('posterior_table2_malami.mat','file')
    fprintf('Posterior .mat files were not found. Running Bayesian update first...\n');
    run_bayes_update_table2_fib_malami;
end

Sfib = load('posterior_table2_fib.mat');
Smal = load('posterior_table2_malami.mat');

post_fib = Sfib.post_fib;
post_malami = Smal.post_malami;

% -------------------------------------------------------------------------
% Observation used in the paper
% -------------------------------------------------------------------------
t_obs = 42;          % years
y_obs = 14.78;       % mm
s_obs = 5.3;         % mm, approx. 14.78*0.36

% Time range for Fig. 8
t_years = linspace(0.1,100,300);

% -------------------------------------------------------------------------
% Posterior predictive model evaluations
% Each row is one posterior sample; each column is one time value.
% -------------------------------------------------------------------------
[x_fib, aux_fib] = carbonation_xc_table2('fib', post_fib.samples, t_years, post_fib.priors);
[x_mal, aux_mal] = carbonation_xc_table2('malami', post_malami.samples, t_years, post_malami.priors);

% Remove non-finite samples, just in case a tested distribution creates them.
x_fib = x_fib(all(isfinite(x_fib),2),:);
x_mal = x_mal(all(isfinite(x_mal),2),:);

% Mean and standard deviation bands used in the manuscript figure.
fib_mean = mean(x_fib,1,'omitnan');
fib_std  = std(x_fib,0,1,'omitnan');
mal_mean = mean(x_mal,1,'omitnan');
mal_std  = std(x_mal,0,1,'omitnan');

% Optional 95% credible intervals for checking. They are saved to CSV but not
% plotted by default, because Fig. 8 in the paper uses mean +/- standard deviation.
fib_q025 = quantile(x_fib,0.025,1);
fib_q975 = quantile(x_fib,0.975,1);
mal_q025 = quantile(x_mal,0.025,1);
mal_q975 = quantile(x_mal,0.975,1);

% -------------------------------------------------------------------------
% Save data behind the figure
% -------------------------------------------------------------------------
T = table(t_years(:), ...
    fib_mean(:), fib_std(:), fib_q025(:), fib_q975(:), ...
    mal_mean(:), mal_std(:), mal_q025(:), mal_q975(:), ...
    'VariableNames', {'t_years', ...
    'fib_mean_mm','fib_std_mm','fib_q025_mm','fib_q975_mm', ...
    'malami_mean_mm','malami_std_mm','malami_q025_mm','malami_q975_mm'});
writetable(T, fullfile(outDir,'Figure8_posterior_carbonation_depth_data.csv'));

% -------------------------------------------------------------------------
% Plot Fig. 8
% -------------------------------------------------------------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 8]);
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Common axes formatting
xLimits = [0 100];
yLimits = [0 30];
yTicks = 0:3:30;

% ---- fib Bulletin 34 panel
nexttile; hold on; box on;
plot(t_years, fib_mean, 'k-', 'LineWidth', 1.6);
plot(t_years, fib_mean + fib_std, 'k--', 'LineWidth', 1.0);
plot(t_years, max(fib_mean - fib_std,0), 'k--', 'LineWidth', 1.0, 'HandleVisibility','off');
errorbar(t_obs, y_obs, s_obs, 'ko', 'MarkerFaceColor','k', 'MarkerSize',4.5, 'CapSize',6, 'LineWidth',1.0);
xlabel('Time, t [years]');
ylabel('Carbonation depth, x_c(t) [mm]');
title('(a) fib Bulletin 34', 'FontWeight','normal');
xlim(xLimits); ylim(yLimits); yticks(yTicks);
legend({'Mean','Mean \pm Std','On-site data'}, 'Location','northwest', 'Box','off');
set(gca,'FontName','Times New Roman','FontSize',9,'LineWidth',0.8);

% ---- Malami panel
nexttile; hold on; box on;
plot(t_years, mal_mean, 'k-', 'LineWidth', 1.6);
plot(t_years, mal_mean + mal_std, 'k--', 'LineWidth', 1.0);
plot(t_years, max(mal_mean - mal_std,0), 'k--', 'LineWidth', 1.0, 'HandleVisibility','off');
errorbar(t_obs, y_obs, s_obs, 'ko', 'MarkerFaceColor','k', 'MarkerSize',4.5, 'CapSize',6, 'LineWidth',1.0);
xlabel('Time, t [years]');
ylabel('Carbonation depth, x_c(t) [mm]');
title('(b) Malami', 'FontWeight','normal');
xlim(xLimits); ylim(yLimits); yticks(yTicks);
legend({'Mean','Mean \pm Std','On-site data'}, 'Location','northwest', 'Box','off');
set(gca,'FontName','Times New Roman','FontSize',9,'LineWidth',0.8);

% -------------------------------------------------------------------------
% Export
% -------------------------------------------------------------------------
pdfName = fullfile(outDir,'Figure8_posterior_carbonation_depth.pdf');
pngName = fullfile(outDir,'Figure8_posterior_carbonation_depth.png');

try
    exportgraphics(fig, pdfName, 'ContentType','vector');
    exportgraphics(fig, pngName, 'Resolution', 300);
catch
    % Fallback for older MATLAB versions
    print(fig, pdfName, '-dpdf', '-painters');
    print(fig, pngName, '-dpng', '-r300');
end

fprintf('\nFigure 8 saved to:\n  %s\n  %s\n', pdfName, pngName);
fprintf('Figure data saved to:\n  %s\n', fullfile(outDir,'Figure8_posterior_carbonation_depth_data.csv'));

% Quick check at inspection time
[xfib42,~] = carbonation_xc_table2('fib', post_fib.samples, t_obs, post_fib.priors);
[xmal42,~] = carbonation_xc_table2('malami', post_malami.samples, t_obs, post_malami.priors);
fprintf('\nPosterior predictive at t = %.0f years:\n', t_obs);
fprintf('  fib Bulletin 34: mean = %.2f mm, std = %.2f mm, CoV = %.2f\n', ...
    mean(xfib42,'omitnan'), std(xfib42,0,'omitnan'), std(xfib42,0,'omitnan')/mean(xfib42,'omitnan'));
fprintf('  Malami:         mean = %.2f mm, std = %.2f mm, CoV = %.2f\n', ...
    mean(xmal42,'omitnan'), std(xmal42,0,'omitnan'), std(xmal42,0,'omitnan')/mean(xmal42,'omitnan'));
