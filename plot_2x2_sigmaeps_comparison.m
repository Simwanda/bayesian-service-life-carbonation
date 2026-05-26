% PLOT_2X2_SIGMAEPS_COMPARISON
% Creates a 2x2 panel figure:
% (a) fib,    sigma_eps = 1.0 mm
% (b) Malami, sigma_eps = 1.0 mm
% (c) fib,    sigma_eps = 2.5 mm
% (d) Malami, sigma_eps = 2.5 mm
%
% It uses the saved files:
%   figs_multiplicative_loop_1mm/combined_loop_results_multiplicative.mat
%   figs_multiplicative_loop_25mm/combined_loop_results_multiplicative.mat

clear; clc; close all;

% -------------------------------------------------------------------------
% Folders
% -------------------------------------------------------------------------
folder_1mm  = fullfile(pwd,'figs_multiplicative_loop_1mm');
folder_25mm = fullfile(pwd,'figs_multiplicative_loop_25mm');

file_1mm  = fullfile(folder_1mm ,'combined_loop_results_multiplicative.mat');
file_25mm = fullfile(folder_25mm,'combined_loop_results_multiplicative.mat');

if ~isfile(file_1mm)
    error('File not found: %s', file_1mm);
end

if ~isfile(file_25mm)
    error('File not found: %s', file_25mm);
end

% -------------------------------------------------------------------------
% Load data
% -------------------------------------------------------------------------
S1  = load(file_1mm);
S25 = load(file_25mm);

% Measured data
% These should be the same in both files, so we take from S1
y_values = S1.y_values(:);
t       = S1.t(:)';
t_obs   = 42;

y_mean = mean(y_values);
y_std  = std(y_values);

% -------------------------------------------------------------------------
% Combined posterior predictive summaries
% sigma_eps = 1.0 mm
% -------------------------------------------------------------------------
mean_fib_1  = mean(S1.meanCurves_fib,1);
std_fib_1   = std(S1.meanCurves_fib,0,1);

mean_mal_1  = mean(S1.meanCurves_mal,1);
std_mal_1   = std(S1.meanCurves_mal,0,1);

% -------------------------------------------------------------------------
% Combined posterior predictive summaries
% sigma_eps = 2.5 mm
% -------------------------------------------------------------------------
mean_fib_25 = mean(S25.meanCurves_fib,1);
std_fib_25  = std(S25.meanCurves_fib,0,1);

mean_mal_25 = mean(S25.meanCurves_mal,1);
std_mal_25  = std(S25.meanCurves_mal,0,1);

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 22 16]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% -------- (a) fib, sigma_eps = 1.0 mm --------
nexttile;
plot_one_panel(t, mean_fib_1, std_fib_1, t_obs, y_mean, y_std, ...
    '(a) fib Bulletin 34, \sigma_\epsilon = 1.0 mm');

% -------- (b) Malami, sigma_eps = 1.0 mm -----
nexttile;
plot_one_panel(t, mean_mal_1, std_mal_1, t_obs, y_mean, y_std, ...
    '(b) Malami, \sigma_\epsilon = 1.0 mm');

% -------- (c) fib, sigma_eps = 2.5 mm --------
nexttile;
plot_one_panel(t, mean_fib_25, std_fib_25, t_obs, y_mean, y_std, ...
    '(c) fib Bulletin 34, \sigma_\epsilon = 2.5 mm');

% -------- (d) Malami, sigma_eps = 2.5 mm -----
nexttile;
plot_one_panel(t, mean_mal_25, std_mal_25, t_obs, y_mean, y_std, ...
    '(d) Malami, \sigma_\epsilon = 2.5 mm');

% Optional overall title
%title(tl,'Comparison of combined posterior predictive carbonation depth','FontWeight','bold');

% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------
outDir = fullfile(pwd,'figs_sigmaeps_comparison');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

pngFile = fullfile(outDir,'combined_posterior_predictive_2x2_sigmaeps_comparison.png');
pdfFile = fullfile(outDir,'combined_posterior_predictive_2x2_sigmaeps_comparison.pdf');

exportgraphics(fig,pngFile,'Resolution',300);
exportgraphics(fig,pdfFile,'ContentType','vector');

fprintf('Saved:\n  %s\n  %s\n', pngFile, pdfFile);

% =========================================================================
% Local function
% =========================================================================
function plot_one_panel(t, m, s, t_obs, y_mean, y_std, titleStr)

    hold on; box on; grid on;

    % Mean
    plot(t, m, 'k-', 'LineWidth', 2.0);

    % Mean +/- std
    plot(t, m + s, 'k--', 'LineWidth', 1.2);
    plot(t, max(m - s, 0), 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');

    % Measured mean +/- std
    errorbar(t_obs, y_mean, y_std, 'ko', ...
        'MarkerFaceColor','k', 'CapSize',6, 'LineWidth',1.0);

    xlabel('Time, t [years]');
    ylabel('Carbonation depth, x_c(t) [mm]');
    title(titleStr,'Interpreter','tex');

    xlim([0 100]);
    ylim([0 35]);

    legend({'Mean of posterior means','Mean \pm Std across updates','Measured mean \pm std'}, ...
        'Location','northwest');
end