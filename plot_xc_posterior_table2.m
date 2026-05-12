%PLOT_XC_POSTERIOR_TABLE2 Plot posterior predictive x_c(t) for Table 2 update.
clear; clc; close all;
outDir = fullfile(pwd,'figs'); if ~exist(outDir,'dir'); mkdir(outDir); end
Sfib = load('posterior_table2_fib.mat');
Smal = load('posterior_table2_malami.mat');
post_fib = Sfib.post_fib; post_mal = Smal.post_malami;

t = linspace(0.1,100,250);
[xfib,~] = carbonation_xc_table2('fib', post_fib.samples, t, post_fib.priors);
[xmal,~] = carbonation_xc_table2('malami', post_mal.samples, t, post_mal.priors);

mfib=mean(xfib,1); sfib=std(xfib,0,1);
mmal=mean(xmal,1); smal=std(xmal,0,1);

fig=figure('Color','w','Units','centimeters','Position',[2 2 18 8]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for k=1:2
    nexttile; hold on; box on; grid on;
    if k==1
        m=mfib; s=sfib; titleStr='fib Bulletin 34';
    else
        m=mmal; s=smal; titleStr='Malami';
    end
    plot(t,m,'k-','LineWidth',1.6);
    plot(t,m+s,'k--','LineWidth',1.0);
    plot(t,max(m-s,0),'k--','LineWidth',1.0,'HandleVisibility','off');
    errorbar(42,14.78,5.3,'ko','MarkerFaceColor','k','CapSize',6);
    xlabel('Time, t [years]'); ylabel('Carbonation depth, x_c(t) [mm]');
    title(titleStr); xlim([0 100]); ylim([0 35]);
    legend({'Mean','Mean +/- Std','On-site data'},'Location','northwest');
end
exportgraphics(fig, fullfile(outDir,'xc_time_posterior_table2.pdf'), 'ContentType','vector');
fprintf('Saved %s\n', fullfile(outDir,'xc_time_posterior_table2.pdf'));
