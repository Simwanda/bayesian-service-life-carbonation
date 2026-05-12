% ==========================================================
% plot_xc_mean_pm_std_posterior.m
% 2x2 figure (standard + log-log) for fib and Malami:
%   - Mean
%   - Mean ± Std
%   - On-site data at 42y (mean ± std)
%
% Uses:
%   - posterior samples for theta_x
%   - fixed A0 and n0
%   - additive residual scatter eps ~ N(0, sigma_eps^2)
%
% Output: figs/xc_time_mean_pm_std_posterior.pdf (vector)
% ==========================================================
clear; clc; close all;

% ------------ paths ------------
outDir = fullfile(pwd,'figs');
if ~exist(outDir,'dir'); mkdir(outDir); end
outFile = fullfile(outDir,'xc_time_mean_pm_std_posterior.pdf');

% ------------ load posterior samples ------------
Sfib = load('posterior_samples_fib.mat');    
post_fib = Sfib.post_fib;

Smal = load('posterior_samples_malami.mat'); 
if isfield(Smal,'post_mal')
    post_mal = Smal.post_mal;
else
    post_mal = Smal.post_malami;
end

% ------------ fixed model parameters ------------
A0_fib = post_fib.A0;
n0_fib = post_fib.n0;

A0_mal = post_mal.A0;
n0_mal = post_mal.n0;

% ------------ on-site measurement (mean ± 1 std) ------------
t_meas  = 42;          % years
xc_mean = 14.8;        % mm
xc_cov  = 0.36;        % -
xc_std  = xc_mean * xc_cov;

% ------------ predictive residual scatter (mm) ------------
sigma_eps = 2.5;

% ------------ time grids ------------
t_lin = linspace(1,100,400)';   % linear scale
t_log = logspace(0,2,250)';     % log scale (1..100)

% ------------ predictive mean/std ------------
[mFib_lin, sFib_lin] = meanStd_xc(post_fib, A0_fib, n0_fib, sigma_eps, t_lin);
[mMal_lin, sMal_lin] = meanStd_xc(post_mal, A0_mal, n0_mal, sigma_eps, t_lin);

[mFib_log, sFib_log] = meanStd_xc(post_fib, A0_fib, n0_fib, sigma_eps, t_log);
[mMal_log, sMal_log] = meanStd_xc(post_mal, A0_mal, n0_mal, sigma_eps, t_log);

% ------------ plot 2x2 figure ------------
fig = figure('Color','w');
set(fig,'Units','centimeters','Position',[2 2 18 14]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

fontName = 'Times New Roman';
fs = 8;

% ---------------- (1,1) fib standard ----------------
nexttile; hold on; box on; grid on;
set(gca,'FontName',fontName,'FontSize',fs,'LineWidth',1);

plot(t_lin, mFib_lin, 'k-', 'LineWidth',2);
plot(t_lin, mFib_lin + sFib_lin, 'k--', 'LineWidth',1.5);
plot(t_lin, mFib_lin - sFib_lin, 'k--', 'LineWidth',1.5, ...
    'HandleVisibility','off');

errorbar(t_meas, xc_mean, xc_std, 'ko', 'MarkerFaceColor','k', ...
    'LineWidth',1.2, 'CapSize',8);

xlabel('Time, $t$ [years]','Interpreter','latex');
ylabel('Carbonation depth, $x_c(t)$ [mm]','Interpreter','latex');
title('fib model','Interpreter','latex');

ylim([0 40]); xlim([0 100]);
legend({'Mean','Mean $\pm$ Std','On-site data'}, ...
    'Interpreter','latex','Location','northwest');

% ---------------- (1,2) Malami standard ----------------
nexttile; hold on; box on; grid on;
set(gca,'FontName',fontName,'FontSize',fs,'LineWidth',1);

plot(t_lin, mMal_lin, 'k-', 'LineWidth',2);
plot(t_lin, mMal_lin + sMal_lin, 'k--', 'LineWidth',1.5);
plot(t_lin, mMal_lin - sMal_lin, 'k--', 'LineWidth',1.5, ...
    'HandleVisibility','off');

errorbar(t_meas, xc_mean, xc_std, 'ko', 'MarkerFaceColor','k', ...
    'LineWidth',1.2, 'CapSize',8);

xlabel('Time, $t$ [years]','Interpreter','latex');
ylabel('Carbonation depth, $x_c(t)$ [mm]','Interpreter','latex');
title('Malami model','Interpreter','latex');

ylim([0 40]); xlim([0 100]);
legend({'Mean','Mean $\pm$ Std','On-site data'}, ...
    'Interpreter','latex','Location','northwest');

% ---------------- (2,1) fib log-log ----------------
nexttile; hold on; box on; grid on;
set(gca,'FontName',fontName,'FontSize',fs,'LineWidth',1);
set(gca,'XScale','log','YScale','log');

plot(t_log, mFib_log, 'k-', 'LineWidth',2);
plot(t_log, mFib_log + sFib_log, 'k--', 'LineWidth',1.5);
plot(t_log, max(mFib_log - sFib_log, eps), 'k--', 'LineWidth',1.5, ...
    'HandleVisibility','off');

errorbar(t_meas, xc_mean, xc_std, 'ko', 'MarkerFaceColor','k', ...
    'LineWidth',1.2, 'CapSize',8);

xlabel('Time, $t$ [years]','Interpreter','latex');
ylabel('Carbonation depth, $x_c(t)$ [mm]','Interpreter','latex');

xlim([1 100]); ylim([1 4e1]);
legend({'Mean','Mean $\pm$ Std','On-site data'}, ...
    'Interpreter','latex','Location','northwest');

% ---------------- (2,2) Malami log-log ----------------
nexttile; hold on; box on; grid on;
set(gca,'FontName',fontName,'FontSize',fs,'LineWidth',1);
set(gca,'XScale','log','YScale','log');

plot(t_log, mMal_log, 'k-', 'LineWidth',2);
plot(t_log, mMal_log + sMal_log, 'k--', 'LineWidth',1.5);
plot(t_log, max(mMal_log - sMal_log, eps), 'k--', 'LineWidth',1.5, ...
    'HandleVisibility','off');

errorbar(t_meas, xc_mean, xc_std, 'ko', 'MarkerFaceColor','k', ...
    'LineWidth',1.2, 'CapSize',8);

xlabel('Time, $t$ [years]','Interpreter','latex');
ylabel('Carbonation depth, $x_c(t)$ [mm]','Interpreter','latex');

xlim([1 100]); ylim([1 4e1]);
legend({'Mean','Mean $\pm$ Std','On-site data'}, ...
    'Interpreter','latex','Location','northwest');

% ------------ export (vector PDF) ------------
set(fig,'PaperPositionMode','auto');
exportgraphics(fig, outFile, ...
    'ContentType','vector', 'BackgroundColor','none');
close(fig);

fprintf('Saved: %s\n', outFile);

% ==========================================================
% Helper: predictive mean and std of x_c(t)
% ==========================================================
function [mxc, sxc] = meanStd_xc(post, A0, n0, sigma_eps, tvec)
% Predictive mean and std of carbonation depth
% x_c(t) = theta_x * A0 * t^(0.5 - n0) + eps
% eps is added once per trajectory (not per time point)

rng(1234,'twister');   % reproducible bands

theta = post.theta(:);
Ns = numel(theta);
Nt = numel(tvec);

T = tvec(:)';          % 1 x Nt
E = 0.5 - n0;          % scalar exponent

% deterministic trajectories
xc = (theta * A0) .* (T .^ E);   % Ns x Nt

% add residual scatter ONCE per sample
eps_i = sigma_eps * randn(Ns,1); % Ns x 1
xc = xc + eps_i;                 % broadcast over time

% physical constraint
xc(xc < 0) = 0;

% moments
mxc = mean(xc, 1).';
sxc = std(xc, 0, 1).';
end
