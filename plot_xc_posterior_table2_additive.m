%PLOT_XC_POSTERIOR_TABLE2_ADDITIVE
% Plot posterior predictive carbonation depth for the additive-discrepancy update.

clear; clc; close all;

outDir = fullfile(pwd,'figs');
if ~exist(outDir,'dir'); mkdir(outDir); end

Sfib = load('posterior_table2_fib_additive.mat');
Smal = load('posterior_table2_malami_additive.mat');
post_fib = Sfib.post_fib;
post_mal = Smal.post_malami;

t = linspace(0.1,100,250);
[xfib,~] = carbonation_xc_table2('fib', post_fib.samples, t, post_fib.priors);
[xmal,~] = carbonation_xc_table2('malami', post_mal.samples, t, post_mal.priors);

mfib = mean(xfib,1);
mmal = mean(xmal,1);

sfib_model = std(xfib,0,1);
smal_model = std(xmal,0,1);

% Optional: prior model scatter used for wider predictive band.
% Set usePriorScatter = false if you want only posterior parameter/discrepancy uncertainty.
usePriorScatter = false;
sigma_prior_fib = 4.5; % mm
sigma_prior_mal = 4.3; % mm

if usePriorScatter
    sfib = sqrt(sfib_model.^2 + sigma_prior_fib.^2);
    smal = sqrt(smal_model.^2 + sigma_prior_mal.^2);
else
    sfib = sfib_model;
    smal = smal_model;
end

% On-site data: mean and spatial scatter of individual measurements.
t_obs = 42;
xc_obs_mean = 14.78;
xc_obs_std = 5.3;

[~,idx42] = min(abs(t - t_obs));
fprintf('\nAdditive-discrepancy posterior at %.1f years:\n', t_obs);
fprintf('fib:    mean = %.2f mm, std = %.2f mm, CoV = %.3f\n', mfib(idx42), sfib(idx42), sfib(idx42)/mfib(idx42));
fprintf('Malami: mean = %.2f mm, std = %.2f mm, CoV = %.3f\n', mmal(idx42), smal(idx42), smal(idx42)/mmal(idx42));

fig = figure('Color','w','Units','centimeters','Position',[2 2 18 8]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

for k = 1:2
    nexttile; hold on; box on; grid on;
    if k == 1
        m = mfib; s = sfib; titleStr = 'fib Bulletin 34';
    else
        m = mmal; s = smal; titleStr = 'Malami';
    end
    plot(t,m,'k-','LineWidth',1.6);
    plot(t,m+s,'k--','LineWidth',1.0);
    plot(t,max(m-s,0),'k--','LineWidth',1.0,'HandleVisibility','off');
    errorbar(t_obs,xc_obs_mean,xc_obs_std,'ko','MarkerFaceColor','k','CapSize',6,'LineWidth',1.0);
    xlabel('Time, t [years]');
    ylabel('Carbonation depth, x_c(t) [mm]');
    title(titleStr);
    xlim([0 100]); ylim([0 35]);
    legend({'Posterior mean','Mean +/- Std','On-site data'},'Location','northwest');
end

pdfFile = fullfile(outDir,'xc_time_posterior_table2_additive.pdf');
pngFile = fullfile(outDir,'xc_time_posterior_table2_additive.png');
exportgraphics(fig,pdfFile,'ContentType','vector');
exportgraphics(fig,pngFile,'Resolution',300);

Tfig = table(t(:), mfib(:), sfib(:), mmal(:), smal(:), ...
    'VariableNames', {'t_years','fib_mean','fib_std','malami_mean','malami_std'});
writetable(Tfig,fullfile(outDir,'xc_time_posterior_table2_additive_data.csv'));

fprintf('\nSaved:\n  %s\n  %s\n', pdfFile, pngFile);
