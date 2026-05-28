% RUN_LOOP_POSTERIOR_ALL_MEASUREMENTS_MULTIPLICATIVE
%
% Bayesian updating loop for individual carbonation-depth measurements
% using the MULTIPLICATIVE model uncertainty factor theta_xc.
%
% For each measured carbonation depth y_i:
%   - use y_i as obs.y_mean
%   - run Bayesian updating for fib and Malami
%   - update Table 2 variables jointly with theta_xc
%   - compute posterior predictive carbonation-depth curve
%
% Outputs are saved in the folder:
%   figs_multiplicative_loop/
%
% This script avoids chained dynamic-field indexing such as S.(name)(:),
% which causes "Invalid array indexing" in some MATLAB versions.

clear; clc; close all;

% ------------------------------------------------------------
% Output folder
% ------------------------------------------------------------
outDir = fullfile(pwd,'figs_multiplicative_loop');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% ------------------------------------------------------------
% Individual measured carbonation depths [mm]
% ------------------------------------------------------------
y_values = [ ...
    26, 19.1, 19.8, 9, 5.7, 18.2, 22.6, 18.6, 16.1, 17, ...
    14.2, 11, 16.4, 11, 13.7, 12.2, 13.9, 12.9, 6.1, 15.6, ...
    13.3, 12.6, 5, 21.2, 17.3, 15.7, 15.1, 12.6, 24.4, 26.5, ...
    15.4, 13, 14.2, 21.5, 18.9, 15.5, 11.2, 21.5, 11.5, 21.2, ...
    9.9, 22.1, 6.2, 7.3, 10.9, 8.6, 14.1, 7, 12, 14.4];

nObs = numel(y_values);

fprintf('Number of individual measurements = %d\n', nObs);
fprintf('Measured mean = %.2f mm\n', mean(y_values));
fprintf('Measured std  = %.2f mm\n', std(y_values));
fprintf('Measured CoV  = %.3f\n\n', std(y_values)/mean(y_values));

% ------------------------------------------------------------
% Observation model
% ------------------------------------------------------------
% Each y_i is used as one calibration target.
% sigma_eps controls how strongly each individual observation is fitted.
obs.t_years   = 42;
obs.y_std     = 0.0;
obs.N_eff     = 1;
obs.sigma_eps = 2.5;   % mm, uncertainty of individual calibration target

% ------------------------------------------------------------
% MCMC settings
% ------------------------------------------------------------
% Use quick settings first. Increase for final runs.
mcmc.Nburn   = 1000;
mcmc.Nkeep   = 3000;
mcmc.thin    = 5;
mcmc.rngSeed = 20260512;

% ------------------------------------------------------------
% Updated variables
% MULTIPLICATIVE model uncertainty:
% theta_xc is updated.
% delta_xc is NOT used.
% ------------------------------------------------------------
update_fib = {'Racc_inv','kt','eps_t','tc','bc', ...
              'CO2ppm','RH','tw','bw','pSR','theta_xc'};

update_malami = {'W_C','tc','bc','CO2ppm','RH', ...
                 'tw','bw','pSR','T','theta_xc'};

propStd_fib = [0.07 0.07 0.06 0.08 0.05 ...
               0.05 0.05 0.07 0.06 0.08 0.08];

propStd_malami = [0.07 0.08 0.05 0.05 0.05 ...
                  0.07 0.06 0.08 0.04 0.08];

% ------------------------------------------------------------
% Time vector for posterior prediction
% ------------------------------------------------------------
t = linspace(0.1,100,250);
[~,idx42] = min(abs(t - obs.t_years));

% Storage for posterior predictive curves
meanCurves_fib = zeros(nObs,numel(t));
stdCurves_fib  = zeros(nObs,numel(t));

meanCurves_mal = zeros(nObs,numel(t));
stdCurves_mal  = zeros(nObs,numel(t));

summaryRows = table();

% To avoid huge memory, keep only a subset of posterior samples per run.
maxCombinePerRun = 2000;

combined_fib = struct();
combined_mal = struct();

% ------------------------------------------------------------
% Main loop
% ------------------------------------------------------------
for i = 1:nObs

    fprintf('\n============================================================\n');
    fprintf('MULTIPLICATIVE update for measurement %d/%d: y = %.2f mm\n', i, nObs, y_values(i));
    fprintf('============================================================\n');

    obs.y_mean = y_values(i);

    % Use different seed for each measurement
    mcmc.rngSeed = 20260512 + i;

    % -----------------------------
    % fib update
    % -----------------------------
    mcmc.propStd = propStd_fib;
    post_fib_i = bayes_update_table2_params('fib', obs, mcmc, update_fib);

    [xfib_i,auxFib_i] = carbonation_xc_table2( ...
        'fib', post_fib_i.samples, t, post_fib_i.priors);

    meanCurves_fib(i,:) = mean(xfib_i,1);
    stdCurves_fib(i,:)  = std(xfib_i,0,1);

    post_fib_i.samples.xc42 = xfib_i(:,idx42);
    post_fib_i.samples = add_aux_if_available(post_fib_i.samples, auxFib_i);

    combined_fib = append_sample_subset(combined_fib, post_fib_i.samples, maxCombinePerRun);

    % -----------------------------
    % Malami update
    % -----------------------------
    mcmc.propStd = propStd_malami;
    post_mal_i = bayes_update_table2_params('malami', obs, mcmc, update_malami);

    [xmal_i,auxMal_i] = carbonation_xc_table2( ...
        'malami', post_mal_i.samples, t, post_mal_i.priors);

    meanCurves_mal(i,:) = mean(xmal_i,1);
    stdCurves_mal(i,:)  = std(xmal_i,0,1);

    post_mal_i.samples.xc42 = xmal_i(:,idx42);
    post_mal_i.samples = add_aux_if_available(post_mal_i.samples, auxMal_i);

    combined_mal = append_sample_subset(combined_mal, post_mal_i.samples, maxCombinePerRun);

    % -----------------------------
    % Store short summary
    % -----------------------------
    row = table( ...
        i, y_values(i), ...
        meanCurves_fib(i,idx42), stdCurves_fib(i,idx42), ...
        meanCurves_mal(i,idx42), stdCurves_mal(i,idx42), ...
        post_fib_i.accRate, post_mal_i.accRate, ...
        'VariableNames', { ...
            'ObsID','Measured_xc_mm', ...
            'fib_mean_xc42','fib_std_xc42', ...
            'malami_mean_xc42','malami_std_xc42', ...
            'fib_acceptance','malami_acceptance'});

    summaryRows = [summaryRows; row];

end

% ------------------------------------------------------------
% Save loop summary
% ------------------------------------------------------------
summaryFile = fullfile(outDir,'loop_posterior_predictive_summary_multiplicative.csv');
writetable(summaryRows,summaryFile);
fprintf('\nSaved loop summary: %s\n', summaryFile);

% ------------------------------------------------------------
% Plot all posterior mean curves: fib
% ------------------------------------------------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 14 9]);
hold on; box on; grid on;

for i = 1:nObs
    plot(t,meanCurves_fib(i,:),'Color',[0.55 0.55 0.55],'LineWidth',0.8);
end

plot(t,mean(meanCurves_fib,1),'k-','LineWidth',2.0);
errorbar(obs.t_years,mean(y_values),std(y_values),'ko', ...
    'MarkerFaceColor','k','CapSize',6,'LineWidth',1.0);

xlabel('Time, t [years]');
ylabel('Carbonation depth, x_c(t) [mm]');
title('fib Bulletin 34: multiplicative posterior mean curves');
legend({'Individual posterior means','Average posterior mean','Measured mean +/- std'}, ...
    'Location','northwest');

xlim([0 100]);
ylim([0 35]);

exportgraphics(fig, fullfile(outDir,'all_posterior_mean_curves_fib_multiplicative.pdf'), ...
    'ContentType','vector');
exportgraphics(fig, fullfile(outDir,'all_posterior_mean_curves_fib_multiplicative.png'), ...
    'Resolution',300);

% ------------------------------------------------------------
% Plot all posterior mean curves: Malami
% ------------------------------------------------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 14 9]);
hold on; box on; grid on;

for i = 1:nObs
    plot(t,meanCurves_mal(i,:),'Color',[0.55 0.55 0.55],'LineWidth',0.8);
end

plot(t,mean(meanCurves_mal,1),'k-','LineWidth',2.0);
errorbar(obs.t_years,mean(y_values),std(y_values),'ko', ...
    'MarkerFaceColor','k','CapSize',6,'LineWidth',1.0);

xlabel('Time, t [years]');
ylabel('Carbonation depth, x_c(t) [mm]');
title('Malami: multiplicative posterior mean curves');
legend({'Individual posterior means','Average posterior mean','Measured mean +/- std'}, ...
    'Location','northwest');

xlim([0 100]);
ylim([0 35]);

exportgraphics(fig, fullfile(outDir,'all_posterior_mean_curves_malami_multiplicative.pdf'), ...
    'ContentType','vector');
exportgraphics(fig, fullfile(outDir,'all_posterior_mean_curves_malami_multiplicative.png'), ...
    'Resolution',300);

% ------------------------------------------------------------
% Figure 8: Combined posterior predictive mean +/- std
% sigma_eps = 2.5 mm only
% Former panels (c) and (d) are now shown as panels (a) and (b).
% ------------------------------------------------------------

mean_fib_all = mean(meanCurves_fib,1);
std_fib_all  = std(meanCurves_fib,0,1);

mean_mal_all = mean(meanCurves_mal,1);
std_mal_all  = std(meanCurves_mal,0,1);

% Set default interpreters for clean labels
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fig = figure('Color','w','Units','centimeters','Position',[2 2 18 7.5]);
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Common limits
xLimits = [0 100];
yLimits = [0 35];

% -----------------------------
% (a) fib Bulletin 34
% -----------------------------
nexttile;
hold on; box on; grid on;

plot(t,mean_fib_all,'k-','LineWidth',2.0);
plot(t,mean_fib_all + std_fib_all,'k--','LineWidth',1.2);
plot(t,max(mean_fib_all - std_fib_all,0),'k--','LineWidth',1.2, ...
    'HandleVisibility','off');

errorbar(obs.t_years,mean(y_values),std(y_values),'ko', ...
    'MarkerFaceColor','k','MarkerSize',4.5, ...
    'CapSize',6,'LineWidth',1.0);

xlabel('Time, $t$ [years]');
ylabel('Carbonation depth, $x_c(t)$ [mm]');
title('(a) fib Bulletin 34','FontWeight','normal');

legend({'Mean of posterior means','Mean $\pm$ Std across updates', ...
    'Measured mean $\pm$ std'}, ...
    'Location','northwest','Box','on');

xlim(xLimits);
ylim(yLimits);
set(gca,'FontName','Times New Roman','FontSize',9,'LineWidth',0.8);

% -----------------------------
% (b) Malami
% -----------------------------
nexttile;
hold on; box on; grid on;

plot(t,mean_mal_all,'k-','LineWidth',2.0);
plot(t,mean_mal_all + std_mal_all,'k--','LineWidth',1.2);
plot(t,max(mean_mal_all - std_mal_all,0),'k--','LineWidth',1.2, ...
    'HandleVisibility','off');

errorbar(obs.t_years,mean(y_values),std(y_values),'ko', ...
    'MarkerFaceColor','k','MarkerSize',4.5, ...
    'CapSize',6,'LineWidth',1.0);

xlabel('Time, $t$ [years]');
ylabel('Carbonation depth, $x_c(t)$ [mm]');
title('(b) Malami','FontWeight','normal');

legend({'Mean of posterior means','Mean $\pm$ Std across updates', ...
    'Measured mean $\pm$ std'}, ...
    'Location','northwest','Box','on');

xlim(xLimits);
ylim(yLimits);
set(gca,'FontName','Times New Roman','FontSize',9,'LineWidth',0.8);

% Export new Figure 8
exportgraphics(fig, fullfile(outDir,'Figure8_posterior_carbonation_depth.pdf'), ...
    'ContentType','vector');
exportgraphics(fig, fullfile(outDir,'Figure8_posterior_carbonation_depth.png'), ...
    'Resolution',300);

% Also keep the older individual combined exports, if needed
exportgraphics(fig, fullfile(outDir,'combined_posterior_predictive_2panel_25mm.pdf'), ...
    'ContentType','vector');
exportgraphics(fig, fullfile(outDir,'combined_posterior_predictive_2panel_25mm.png'), ...
    'Resolution',300);

fprintf('\nSaved updated two-panel Figure 8:\n');
fprintf('  %s\n', fullfile(outDir,'Figure8_posterior_carbonation_depth.pdf'));
fprintf('  %s\n', fullfile(outDir,'Figure8_posterior_carbonation_depth.png'));

% ------------------------------------------------------------
% Combined posterior sample analysis
% Distribution fitting and correlation matrices
% ------------------------------------------------------------
fprintf('\nFitting combined posterior distributions: fib\n');
Tfit_fib = fit_combined_distributions(combined_fib,'fib');
writetable(Tfit_fib,fullfile(outDir,'combined_fit_summary_fib_multiplicative.csv'));

fprintf('Fitting combined posterior distributions: Malami\n');
Tfit_mal = fit_combined_distributions(combined_mal,'malami');
writetable(Tfit_mal,fullfile(outDir,'combined_fit_summary_malami_multiplicative.csv'));

% Correlation matrices
[Tcorr_fib_P,Tcorr_fib_S] = compute_corr_tables(combined_fib);
[Tcorr_mal_P,Tcorr_mal_S] = compute_corr_tables(combined_mal);

writetable(Tcorr_fib_P,fullfile(outDir,'combined_corr_fib_pearson_multiplicative.csv'), ...
    'WriteRowNames',true);
writetable(Tcorr_fib_S,fullfile(outDir,'combined_corr_fib_spearman_multiplicative.csv'), ...
    'WriteRowNames',true);

writetable(Tcorr_mal_P,fullfile(outDir,'combined_corr_malami_pearson_multiplicative.csv'), ...
    'WriteRowNames',true);
writetable(Tcorr_mal_S,fullfile(outDir,'combined_corr_malami_spearman_multiplicative.csv'), ...
    'WriteRowNames',true);

% Save MATLAB data
save(fullfile(outDir,'combined_loop_results_multiplicative.mat'), ...
    'y_values','t', ...
    'meanCurves_fib','stdCurves_fib', ...
    'meanCurves_mal','stdCurves_mal', ...
    'combined_fib','combined_mal', ...
    'summaryRows','Tfit_fib','Tfit_mal');

fprintf('\nAll multiplicative loop results saved in folder: %s\n', outDir);

% ========================================================================
% Local helper functions
% ========================================================================

function S = add_aux_if_available(S,aux)
% Add derived A and n to samples if carbonation_xc_table2 returns them.

    if isempty(aux)
        return;
    end

    if isstruct(aux)
        if isfield(aux,'A')
            A = aux.A;
            S.A = A(:);
        end
        if isfield(aux,'n')
            n = aux.n;
            S.n = n(:);
        end
    end
end

function Sout = append_sample_subset(Sout,Sin,maxN)
% Append a random subset of samples from Sin into Sout.
% Avoids chained indexing for MATLAB compatibility.

    names = fieldnames(Sin);
    firstName = names{1};
    xFirst = Sin.(firstName);
    n = numel(xFirst);

    if n > maxN
        idx = randperm(n,maxN);
    else
        idx = 1:n;
    end

    for j = 1:numel(names)
        name = names{j};
        x = Sin.(name);
        x = x(:);
        x = x(idx);

        if ~isfield(Sout,name)
            Sout.(name) = x;
        else
            old = Sout.(name);
            Sout.(name) = [old(:); x(:)];
        end
    end
end

function Tfit = fit_combined_distributions(S,modelName)
% Fit simple candidate distributions to each combined posterior variable.
%
% Candidate distributions:
%   - Normal
%   - Lognormal, only if all samples are positive
%
% The selected distribution is the one with lower AIC.

    names = fieldnames(S);

    Parameter = {};
    Model = {};
    BestPDF = {};
    Mean = [];
    Std = [];
    CoV = [];
    Median = [];
    P2p5 = [];
    P97p5 = [];
    AIC_Normal = [];
    AIC_Lognormal = [];

    for j = 1:numel(names)

        name = names{j};
        x = S.(name);
        x = x(:);
        x = x(isfinite(x));

        if numel(x) < 10
            continue;
        end

        if std(x) <= 0
            continue;
        end

        mu = mean(x);
        sig = std(x);
        covx = sig / max(abs(mu),eps);

        med = median(x);
        q = quantile(x,[0.025 0.975]);

        % Normal fit
        try
            pdN = fitdist(x,'Normal');
            pN = pdf(pdN,x);
            pN(pN < realmin) = realmin;
            logLN = sum(log(pN));
            aicN = 2*2 - 2*logLN;
        catch
            aicN = NaN;
        end

        % Lognormal fit, only for positive data
        if all(x > 0)
            try
                pdL = fitdist(x,'Lognormal');
                pL = pdf(pdL,x);
                pL(pL < realmin) = realmin;
                logLL = sum(log(pL));
                aicL = 2*2 - 2*logLL;
            catch
                aicL = NaN;
            end
        else
            aicL = NaN;
        end

        % Select best
        if isnan(aicL)
            best = 'Normal';
        elseif isnan(aicN)
            best = 'Lognormal';
        elseif aicL < aicN
            best = 'Lognormal';
        else
            best = 'Normal';
        end

        Model{end+1,1} = modelName;
        Parameter{end+1,1} = name;
        BestPDF{end+1,1} = best;
        Mean(end+1,1) = mu;
        Std(end+1,1) = sig;
        CoV(end+1,1) = covx;
        Median(end+1,1) = med;
        P2p5(end+1,1) = q(1);
        P97p5(end+1,1) = q(2);
        AIC_Normal(end+1,1) = aicN;
        AIC_Lognormal(end+1,1) = aicL;

    end

    Tfit = table(Model,Parameter,BestPDF,Mean,Std,CoV,Median,P2p5,P97p5, ...
        AIC_Normal,AIC_Lognormal);
end

function [TPearson,TSpearman] = compute_corr_tables(S)
% Compute Pearson and Spearman correlation matrices from combined samples.
% Compatible with older MATLAB versions.

    names0 = fieldnames(S);

    % Keep only numeric variables with finite values and non-zero variance
    keep = false(numel(names0),1);

    for j = 1:numel(names0)
        name = names0{j};
        x = S.(name);
        x = x(:);
        x = x(isfinite(x));

        if numel(x) > 2
            keep(j) = std(x,0) > 0;
        end
    end

    names = names0(keep);

    if isempty(names)
        error('No valid variables found for correlation analysis.');
    end

    % Make all variables the same length
    nMin = inf;
    for j = 1:numel(names)
        name = names{j};
        x = S.(name);
        x = x(:);
        x = x(isfinite(x));
        nMin = min(nMin,numel(x));
    end

    X = zeros(nMin,numel(names));

    for j = 1:numel(names)
        name = names{j};
        x = S.(name);
        x = x(:);
        x = x(isfinite(x));
        X(:,j) = x(1:nMin);
    end

    % Compute correlations
    Rpear = corr(X,'Type','Pearson','Rows','pairwise');
    Rspear = corr(X,'Type','Spearman','Rows','pairwise');

    % Convert variable names to valid MATLAB table names
    validNames = matlab.lang.makeValidName(names);

    TPearson = array2table(Rpear, ...
        'VariableNames',validNames, ...
        'RowNames',validNames);

    TSpearman = array2table(Rspear, ...
        'VariableNames',validNames, ...
        'RowNames',validNames);
end
