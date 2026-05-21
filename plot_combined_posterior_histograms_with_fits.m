% PLOT_COMBINED_POSTERIOR_HISTOGRAMS_WITH_FITS
%
% Plots posterior samples as histograms and overlays fitted PDFs.
%
% Requires:
%   figs/combined_loop_results.mat
%
% Produced by:
%   run_loop_posterior_all_measurements_additive.m
%
% Outputs:
%   figs/posterior_fits_fib/
%   figs/posterior_fits_malami/
%
% Each figure shows:
%   - posterior samples histogram
%   - fitted Normal PDF
%   - fitted Lognormal PDF, if samples are positive

clear; clc; close all;

% ------------------------------------------------------------
% Load combined posterior samples
% ------------------------------------------------------------
outDir = fullfile(pwd,'figs');
dataFile = fullfile(outDir,'combined_loop_results.mat');

if ~isfile(dataFile)
    error('File not found: %s. Run run_loop_posterior_all_measurements_additive first.', dataFile);
end

S = load(dataFile);

combined_fib = S.combined_fib;
combined_mal = S.combined_mal;

% ------------------------------------------------------------
% Output folders
% ------------------------------------------------------------
fibDir = fullfile(outDir,'posterior_fits_fib');
malDir = fullfile(outDir,'posterior_fits_malami');

if ~exist(fibDir,'dir')
    mkdir(fibDir);
end

if ~exist(malDir,'dir')
    mkdir(malDir);
end

% ------------------------------------------------------------
% Plot all variables
% ------------------------------------------------------------
Tfit_fib = plot_all_fits(combined_fib, 'fib', fibDir);
Tfit_mal = plot_all_fits(combined_mal, 'malami', malDir);

% Save fit summaries
writetable(Tfit_fib, fullfile(outDir,'posterior_histogram_fit_summary_fib.csv'));
writetable(Tfit_mal, fullfile(outDir,'posterior_histogram_fit_summary_malami.csv'));

fprintf('\nSaved histogram/PDF plots in:\n');
fprintf('  %s\n', fibDir);
fprintf('  %s\n', malDir);

fprintf('\nSaved fit summaries:\n');
fprintf('  %s\n', fullfile(outDir,'posterior_histogram_fit_summary_fib.csv'));
fprintf('  %s\n', fullfile(outDir,'posterior_histogram_fit_summary_malami.csv'));

% ========================================================================
% Local function
% ========================================================================

function Tfit = plot_all_fits(S, modelName, saveDir)
% Plot histogram and fitted PDFs for all variables in sample structure S.

    names = fieldnames(S);

    Model = {};
    Parameter = {};
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

        if numel(x) < 20
            fprintf('Skipping %s: too few samples.\n', name);
            continue;
        end

        if std(x) <= 0
            fprintf('Skipping %s: zero variance.\n', name);
            continue;
        end

        % Basic posterior statistics
        mu = mean(x);
        sig = std(x);
        covx = sig / max(abs(mu),eps);
        med = median(x);
        q = quantile(x,[0.025 0.975]);

        % Fit Normal
        pdN = [];
        aicN = NaN;

        try
            pdN = fitdist(x,'Normal');
            pN = pdf(pdN,x);
            pN(pN < realmin) = realmin;
            logLN = sum(log(pN));
            aicN = 2*2 - 2*logLN;
        catch ME
            warning('Normal fit failed for %s: %s', name, ME.message);
        end

        % Fit Lognormal only if positive
        pdL = [];
        aicL = NaN;

        if all(x > 0)
            try
                pdL = fitdist(x,'Lognormal');
                pL = pdf(pdL,x);
                pL(pL < realmin) = realmin;
                logLL = sum(log(pL));
                aicL = 2*2 - 2*logLL;
            catch ME
                warning('Lognormal fit failed for %s: %s', name, ME.message);
            end
        end

        % Select best PDF using AIC
        if isnan(aicL)
            best = 'Normal';
        elseif isnan(aicN)
            best = 'Lognormal';
        elseif aicL < aicN
            best = 'Lognormal';
        else
            best = 'Normal';
        end

        % ------------------------------------------------------------
        % Plot histogram and fitted PDFs
        % ------------------------------------------------------------
        fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
        hold on; box on; grid on;

        histogram(x,40,'Normalization','pdf', ...
            'FaceAlpha',0.35, ...
            'EdgeColor','none');

        xmin = min(x);
        xmax = max(x);

        if xmin == xmax
            xx = linspace(xmin-1,xmax+1,300);
        else
            pad = 0.05*(xmax - xmin);
            xx = linspace(xmin-pad,xmax+pad,400);
        end

        legendEntries = {'Posterior samples'};

        if ~isempty(pdN)
            plot(xx,pdf(pdN,xx),'k-','LineWidth',1.8);
            legendEntries{end+1} = sprintf('Normal, AIC = %.1f', aicN);
        end

        if ~isempty(pdL)
            % Lognormal pdf is defined only for positive xx
            xxL = xx;
            yyL = nan(size(xxL));
            idxPos = xxL > 0;
            yyL(idxPos) = pdf(pdL,xxL(idxPos));

            plot(xxL,yyL,'k--','LineWidth',1.8);
            legendEntries{end+1} = sprintf('Lognormal, AIC = %.1f', aicL);
        end

        xlabel(strrep(name,'_','\_'));
        ylabel('Probability density');

        title(sprintf('%s posterior: %s | Best fit: %s', ...
            modelName, strrep(name,'_','\_'), best), ...
            'Interpreter','tex');

        legend(legendEntries,'Location','best');

        % Add text box with statistics
        txt = sprintf(['Mean = %.4g\nStd = %.4g\nCoV = %.3g\n', ...
                       'Median = %.4g\n2.5%% = %.4g\n97.5%% = %.4g'], ...
                       mu, sig, covx, med, q(1), q(2));

        xText = xmin + 0.03*(xmax-xmin);
        yLim = ylim;
        yText = yLim(2) - 0.08*(yLim(2)-yLim(1));

        text(xText,yText,txt, ...
            'VerticalAlignment','top', ...
            'BackgroundColor','w', ...
            'EdgeColor',[0.4 0.4 0.4], ...
            'Margin',5);

        % Save
        safeName = matlab.lang.makeValidName(name);

        pdfFile = fullfile(saveDir, sprintf('%s_posterior_fit_%s.pdf', modelName, safeName));
        pngFile = fullfile(saveDir, sprintf('%s_posterior_fit_%s.png', modelName, safeName));

        exportgraphics(fig,pdfFile,'ContentType','vector');
        exportgraphics(fig,pngFile,'Resolution',300);

        close(fig);

        % Store summary
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

        fprintf('%s - %-12s: best = %s, mean = %.4g, std = %.4g, CoV = %.3g\n', ...
            modelName, name, best, mu, sig, covx);

    end

    Tfit = table(Model,Parameter,BestPDF,Mean,Std,CoV,Median,P2p5,P97p5, ...
        AIC_Normal,AIC_Lognormal);
end