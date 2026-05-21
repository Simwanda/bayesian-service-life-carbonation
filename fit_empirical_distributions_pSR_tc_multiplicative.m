% FIT_EMPIRICAL_DISTRIBUTIONS_PSR_TC_MULTIPLICATIVE
%
% Fits empirical/kernel posterior distributions for pSR and tc
% from the multiplicative Bayesian loop results.
%
% Requires:
%   figs_multiplicative_loop/combined_loop_results_multiplicative.mat
%
% Outputs:
%   figs_multiplicative_loop/empirical_fit_pSR_tc/
%       fib_tc_empirical_pdf_cdf.csv
%       fib_pSR_empirical_pdf_cdf.csv
%       malami_tc_empirical_pdf_cdf.csv
%       malami_pSR_empirical_pdf_cdf.csv
%       *_empirical_fit.png/pdf
%
% The empirical PDF is obtained using kernel density estimation.
% The empirical CDF is obtained directly from posterior samples.

clear; clc; close all;

% ------------------------------------------------------------
% Load combined posterior samples
% ------------------------------------------------------------
outDir = fullfile(pwd,'figs_multiplicative_loop');
dataFile = fullfile(outDir,'combined_loop_results_multiplicative.mat');

if ~isfile(dataFile)
    error('File not found: %s. Run run_loop_posterior_all_measurements_multiplicative first.', dataFile);
end

S = load(dataFile);

combined_fib = S.combined_fib;
combined_mal = S.combined_mal;

% ------------------------------------------------------------
% Output folder
% ------------------------------------------------------------
empDir = fullfile(outDir,'empirical_fit_pSR_tc');
if ~exist(empDir,'dir')
    mkdir(empDir);
end

% ------------------------------------------------------------
% Variables to fit empirically
% ------------------------------------------------------------
varsToFit = {'tc','pSR'};

% ------------------------------------------------------------
% Fit empirical distributions
% ------------------------------------------------------------
Tfib = fit_empirical_for_model(combined_fib,'fib',varsToFit,empDir);
Tmal = fit_empirical_for_model(combined_mal,'malami',varsToFit,empDir);

Tall = [Tfib; Tmal];

writetable(Tfib, fullfile(empDir,'empirical_fit_summary_fib.csv'));
writetable(Tmal, fullfile(empDir,'empirical_fit_summary_malami.csv'));
writetable(Tall, fullfile(empDir,'empirical_fit_summary_ALL.csv'));

disp(Tall)

fprintf('\nSaved empirical distribution results in:\n');
fprintf('  %s\n', empDir);

% ========================================================================
% Local functions
% ========================================================================

function Tsummary = fit_empirical_for_model(S,modelName,varsToFit,empDir)

    Model = {};
    Parameter = {};
    FitType = {};
    Mean = [];
    Std = [];
    CoV = [];
    Median = [];
    P2p5 = [];
    P97p5 = [];
    MinSample = [];
    MaxSample = [];
    Nsamples = [];

    for i = 1:numel(varsToFit)

        varName = varsToFit{i};

        if ~isfield(S,varName)
            warning('%s does not contain variable %s. Skipping.', modelName, varName);
            continue;
        end

        x = S.(varName);
        x = x(:);
        x = x(isfinite(x));

        if numel(x) < 20
            warning('%s %s has too few samples. Skipping.', modelName, varName);
            continue;
        end

        % ------------------------------------------------------------
        % Basic statistics
        % ------------------------------------------------------------
        mu = mean(x);
        sig = std(x);
        covx = sig / max(abs(mu),eps);
        med = median(x);
        q = quantile(x,[0.025 0.975]);

        % ------------------------------------------------------------
        % Kernel PDF
        % ------------------------------------------------------------
        xmin = min(x);
        xmax = max(x);
        pad = 0.05*(xmax - xmin);

        if pad == 0
            pad = 0.1;
        end

        xx = linspace(xmin-pad,xmax+pad,500);

        % KDE. Use support according to variable bounds.
        if strcmpi(varName,'tc')
            lowerBound = 1.0;
            upperBound = 3.0;
        elseif strcmpi(varName,'pSR')
            lowerBound = 0.05;
            upperBound = 0.25;
        else
            lowerBound = xmin;
            upperBound = xmax;
        end

        % Kernel density with bounded support
        try
            [pdf_kde,xx_kde] = ksdensity(x,xx, ...
                'Support',[lowerBound upperBound], ...
                'BoundaryCorrection','reflection');
        catch
            % Fallback for older MATLAB versions
            [pdf_kde,xx_kde] = ksdensity(x,xx);
        end
        
                % ------------------------------------------------------------
        % Empirical CDF
        % ------------------------------------------------------------
        try
            [Femp,Xemp] = ecdf(x);
        catch
            xs = sort(x);
            Xemp = xs;
            Femp = (1:numel(xs))' ./ numel(xs);
        end
        
        % Interpolate CDF onto the PDF grid for export
        % ecdf can return repeated X values, so make them unique first.
        [Xemp_unique, ia] = unique(Xemp, 'last');
        Femp_unique = Femp(ia);
        
        Fgrid = interp1(Xemp_unique, Femp_unique, xx_kde, 'previous', 'extrap');
        Fgrid(Fgrid < 0) = 0;
        Fgrid(Fgrid > 1) = 1;

        % ------------------------------------------------------------
        % Plot histogram + empirical PDF
        % ------------------------------------------------------------
        fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
        hold on; box on; grid on;

        histogram(x,35,'Normalization','pdf', ...
            'FaceAlpha',0.35, ...
            'EdgeColor','none');

        plot(xx_kde,pdf_kde,'k-','LineWidth',2.0);

        xlabel(strrep(varName,'_','\_'));
        ylabel('Probability density');

        title(sprintf('%s posterior empirical fit: %s', ...
            modelName, strrep(varName,'_','\_')), ...
            'Interpreter','tex');

        legend({'Posterior samples','Empirical KDE PDF'}, ...
            'Location','best');

        txt = sprintf(['Mean = %.4g\nStd = %.4g\nCoV = %.3g\n', ...
                       'Median = %.4g\n2.5%% = %.4g\n97.5%% = %.4g'], ...
                       mu, sig, covx, med, q(1), q(2));

        yLim = ylim;
        xLim = xlim;

        text(xLim(1) + 0.05*(xLim(2)-xLim(1)), ...
             yLim(2) - 0.08*(yLim(2)-yLim(1)), ...
             txt, ...
             'VerticalAlignment','top', ...
             'BackgroundColor','w', ...
             'EdgeColor',[0.4 0.4 0.4], ...
             'Margin',5);

        pngFile = fullfile(empDir,sprintf('%s_%s_empirical_fit.png',modelName,varName));
        pdfFile = fullfile(empDir,sprintf('%s_%s_empirical_fit.pdf',modelName,varName));

        exportgraphics(fig,pngFile,'Resolution',300);
        exportgraphics(fig,pdfFile,'ContentType','vector');

        close(fig);

        % ------------------------------------------------------------
        % Plot empirical CDF
        % ------------------------------------------------------------
        fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
        hold on; box on; grid on;

        stairs(Xemp,Femp,'k-','LineWidth',1.8);

        xlabel(strrep(varName,'_','\_'));
        ylabel('Empirical CDF');

        title(sprintf('%s posterior empirical CDF: %s', ...
            modelName, strrep(varName,'_','\_')), ...
            'Interpreter','tex');

        xlim([lowerBound upperBound]);
        ylim([0 1]);

        pngFile = fullfile(empDir,sprintf('%s_%s_empirical_cdf.png',modelName,varName));
        pdfFile = fullfile(empDir,sprintf('%s_%s_empirical_cdf.pdf',modelName,varName));

        exportgraphics(fig,pngFile,'Resolution',300);
        exportgraphics(fig,pdfFile,'ContentType','vector');

        close(fig);

        % ------------------------------------------------------------
        % Store summary
        % ------------------------------------------------------------
        Model{end+1,1} = modelName;
        Parameter{end+1,1} = varName;
        FitType{end+1,1} = 'Empirical_KDE';
        Mean(end+1,1) = mu;
        Std(end+1,1) = sig;
        CoV(end+1,1) = covx;
        Median(end+1,1) = med;
        P2p5(end+1,1) = q(1);
        P97p5(end+1,1) = q(2);
        MinSample(end+1,1) = min(x);
        MaxSample(end+1,1) = max(x);
        Nsamples(end+1,1) = numel(x);

        fprintf('%s - %s empirical KDE: mean = %.4g, std = %.4g, CoV = %.3g\n', ...
            modelName, varName, mu, sig, covx);

    end

    Tsummary = table(Model,Parameter,FitType,Mean,Std,CoV,Median, ...
        P2p5,P97p5,MinSample,MaxSample,Nsamples);

end