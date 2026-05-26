% FIT_BEST_CDF_PSR_TC_NO_BETA
%
% Fits simple CDF functions for bounded posterior variables tc and pSR.
% No scaled beta distribution is used.
%
% Candidate CDFs:
%   1. Normal
%   2. Lognormal, if samples are positive
%   3. Weibull, if samples are positive
%   4. Bounded sqrt CDF:
%        F(x) = sqrt((x-a)/(b-a))
%   5. Bounded power CDF:
%        F(x) = ((x-a)/(b-a))^gamma
%
% Input:
%   figs_multiplicative_loop_25mm/combined_loop_results_multiplicative.mat
%
% Output:
%   figs_multiplicative_loop_25mm/best_CDF_pSR_tc_no_beta/

clear; clc; close all;

% -------------------------------------------------------------------------
% Load posterior samples
% -------------------------------------------------------------------------
outDir = fullfile(pwd,'figs_multiplicative_loop_25mm');
dataFile = fullfile(outDir,'combined_loop_results_multiplicative.mat');

if ~isfile(dataFile)
    error('File not found: %s', dataFile);
end

S = load(dataFile);

saveDir = fullfile(outDir,'best_CDF_pSR_tc_no_beta');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% -------------------------------------------------------------------------
% Fit variables
% -------------------------------------------------------------------------
T1 = fit_one_variable(S.combined_fib,'fib','tc',1.0,3.0,saveDir);
T2 = fit_one_variable(S.combined_fib,'fib','pSR',0.05,0.25,saveDir);
T3 = fit_one_variable(S.combined_mal,'Malami','tc',1.0,3.0,saveDir);
T4 = fit_one_variable(S.combined_mal,'Malami','pSR',0.05,0.25,saveDir);

Tall = [T1; T2; T3; T4];

disp(Tall)

writetable(Tall,fullfile(saveDir,'best_CDF_summary_no_beta.csv'));

% -------------------------------------------------------------------------
% Write LaTeX notes
% -------------------------------------------------------------------------
latexFile = fullfile(saveDir,'best_CDF_latex_footnotes_no_beta.txt');
fid = fopen(latexFile,'w');

fprintf(fid,'%% Suggested LaTeX footnotes for fitted CDFs, no beta distribution\n\n');

fprintf(fid,['The bounded variables $t_c$ and $p_{SR}$ were represented using the best-fitting ', ...
             'CDF selected from normal, lognormal, Weibull, bounded square-root and bounded power forms. ', ...
             'The fitted CDF was selected by the smallest Kolmogorov--Smirnov distance to the empirical posterior CDF.\n\n']);

for i = 1:height(Tall)

    modelName = Tall.Model{i};
    varName   = Tall.Parameter{i};
    bestCDF   = Tall.BestCDF{i};

    if strcmpi(bestCDF,'Bounded power')
        fprintf(fid,'%s, $%s$: bounded power CDF $F(x)=((x-a)/(b-a))^{\\gamma}$ with $a=%.4g$, $b=%.4g$, and $\\gamma=%.4g$.\n', ...
            modelName, varName, Tall.Lower(i), Tall.Upper(i), Tall.Gamma(i));

    elseif strcmpi(bestCDF,'Bounded sqrt')
        fprintf(fid,'%s, $%s$: bounded square-root CDF $F(x)=\\sqrt{(x-a)/(b-a)}$ with $a=%.4g$ and $b=%.4g$.\n', ...
            modelName, varName, Tall.Lower(i), Tall.Upper(i));

    elseif strcmpi(bestCDF,'Lognormal')
        fprintf(fid,'%s, $%s$: lognormal CDF with fitted parameters $\\mu_{\\ln}=%.4g$ and $\\sigma_{\\ln}=%.4g$.\n', ...
            modelName, varName, Tall.Param1(i), Tall.Param2(i));

    elseif strcmpi(bestCDF,'Normal')
        fprintf(fid,'%s, $%s$: normal CDF with mean %.4g and standard deviation %.4g.\n', ...
            modelName, varName, Tall.Param1(i), Tall.Param2(i));

    elseif strcmpi(bestCDF,'Weibull')
        fprintf(fid,'%s, $%s$: Weibull CDF with scale %.4g and shape %.4g.\n', ...
            modelName, varName, Tall.Param1(i), Tall.Param2(i));
    end
end

fclose(fid);

fprintf('\nSaved:\n');
fprintf('  %s\n', fullfile(saveDir,'best_CDF_summary_no_beta.csv'));
fprintf('  %s\n', latexFile);

% =========================================================================
% Local function
% =========================================================================
function T = fit_one_variable(S,modelName,varName,a,b,saveDir)

    if ~isfield(S,varName)
        error('%s does not contain %s', modelName, varName);
    end

    x = S.(varName);
    x = x(:);
    x = x(isfinite(x));

    if numel(x) < 20
        error('Too few samples for %s %s', modelName, varName);
    end

    % Keep values inside physical bounds for bounded candidates
    xBound = min(max(x,a),b);

    % Empirical CDF
    xs = sort(xBound);
    n = numel(xs);
    Femp = (1:n)' ./ n;

    % Grid
    xx = linspace(a,b,600);

    % Storage
    candNames = {};
    KS = [];
    RMSE = [];
    Param1 = [];
    Param2 = [];
    GammaList = [];

    % ---------------------------------------------------------------------
    % Candidate 1: Normal
    % ---------------------------------------------------------------------
    try
        pdN = fitdist(x,'Normal');
        FN = cdf(pdN,xs);

        candNames{end+1,1} = 'Normal';
        KS(end+1,1) = max(abs(Femp - FN));
        RMSE(end+1,1) = sqrt(mean((Femp - FN).^2));
        Param1(end+1,1) = pdN.mu;
        Param2(end+1,1) = pdN.sigma;
        GammaList(end+1,1) = NaN;
    catch
    end

    % ---------------------------------------------------------------------
    % Candidate 2: Lognormal
    % ---------------------------------------------------------------------
    if all(x > 0)
        try
            pdL = fitdist(x,'Lognormal');
            FL = cdf(pdL,xs);

            candNames{end+1,1} = 'Lognormal';
            KS(end+1,1) = max(abs(Femp - FL));
            RMSE(end+1,1) = sqrt(mean((Femp - FL).^2));
            Param1(end+1,1) = pdL.mu;
            Param2(end+1,1) = pdL.sigma;
            GammaList(end+1,1) = NaN;
        catch
        end
    end

    % ---------------------------------------------------------------------
    % Candidate 3: Weibull
    % ---------------------------------------------------------------------
    if all(x > 0)
        try
            pdW = fitdist(x,'Weibull');
            FW = cdf(pdW,xs);

            candNames{end+1,1} = 'Weibull';
            KS(end+1,1) = max(abs(Femp - FW));
            RMSE(end+1,1) = sqrt(mean((Femp - FW).^2));
            Param1(end+1,1) = pdW.A;  % scale
            Param2(end+1,1) = pdW.B;  % shape
            GammaList(end+1,1) = NaN;
        catch
        end
    end

    % ---------------------------------------------------------------------
    % Candidate 4: Bounded sqrt CDF
    % F(x) = sqrt((x-a)/(b-a))
    % ---------------------------------------------------------------------
    z = (xs-a)/(b-a);
    z = min(max(z,0),1);

    Fsqrt = sqrt(z);

    candNames{end+1,1} = 'Bounded sqrt';
    KS(end+1,1) = max(abs(Femp - Fsqrt));
    RMSE(end+1,1) = sqrt(mean((Femp - Fsqrt).^2));
    Param1(end+1,1) = a;
    Param2(end+1,1) = b;
    GammaList(end+1,1) = 0.5;

    % ---------------------------------------------------------------------
    % Candidate 5: Bounded power CDF
    % F(x) = ((x-a)/(b-a))^gamma
    % gamma is fitted by minimizing CDF RMSE
    % ---------------------------------------------------------------------
    obj = @(g) sqrt(mean((Femp - z.^g).^2));

    try
        gammaHat = fminbnd(obj,0.05,5.0);
    catch
        gammaHat = 1.0;
    end

    Fpow = z.^gammaHat;

    candNames{end+1,1} = 'Bounded power';
    KS(end+1,1) = max(abs(Femp - Fpow));
    RMSE(end+1,1) = sqrt(mean((Femp - Fpow).^2));
    Param1(end+1,1) = a;
    Param2(end+1,1) = b;
    GammaList(end+1,1) = gammaHat;

    % ---------------------------------------------------------------------
    % Select best by KS distance
    % ---------------------------------------------------------------------
    [~,idxBest] = min(KS);
    bestName = candNames{idxBest};
    bestP1 = Param1(idxBest);
    bestP2 = Param2(idxBest);
    bestGamma = GammaList(idxBest);

    % ---------------------------------------------------------------------
    % Compute best CDF on grid
    % ---------------------------------------------------------------------
    Fbest = evaluate_cdf(bestName,xx,a,b,bestP1,bestP2,bestGamma);

    % ---------------------------------------------------------------------
    % Plot
    % ---------------------------------------------------------------------
    fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
    hold on; box on; grid on;

    stairs(xs,Femp,'Color',[0.45 0.45 0.45],'LineWidth',1.2);
    plot(xx,Fbest,'k-','LineWidth',2.0);

    xlabel(strrep(varName,'_','\_'));
    ylabel('CDF, F(x)');
    title(sprintf('%s %s: empirical CDF and best fitted CDF', ...
        modelName, strrep(varName,'_','\_')), ...
        'Interpreter','tex');

    legend({'Empirical CDF',sprintf('Best fit: %s',bestName)}, ...
        'Location','best');

    txt = sprintf('Best = %s\nKS = %.4g\nRMSE = %.4g', ...
        bestName, KS(idxBest), RMSE(idxBest));

    text(a+0.05*(b-a),0.92,txt, ...
        'VerticalAlignment','top', ...
        'BackgroundColor','w', ...
        'EdgeColor',[0.4 0.4 0.4], ...
        'Margin',5);

    pngFile = fullfile(saveDir,sprintf('%s_%s_best_CDF_no_beta.png',modelName,varName));
    pdfFile = fullfile(saveDir,sprintf('%s_%s_best_CDF_no_beta.pdf',modelName,varName));

    exportgraphics(fig,pngFile,'Resolution',300);
    exportgraphics(fig,pdfFile,'ContentType','vector');
    close(fig);

    % ---------------------------------------------------------------------
    % Save CDF grid
    % ---------------------------------------------------------------------
    Tcdf = table(xx(:),Fbest(:),'VariableNames',{'x','FittedCDF'});
    writetable(Tcdf,fullfile(saveDir,sprintf('%s_%s_best_CDF_values.csv',modelName,varName)));

    % ---------------------------------------------------------------------
    % Summary row
    % ---------------------------------------------------------------------
    Model = {modelName};
    Parameter = {varName};
    BestCDF = {bestName};
    Mean = mean(x);
    Std = std(x);
    CoV = Std / abs(Mean);
    Lower = a;
    Upper = b;
    Param1_out = bestP1;
    Param2_out = bestP2;
    Gamma = bestGamma;
    KS_distance = KS(idxBest);
    RMSE_CDF = RMSE(idxBest);

    T = table(Model,Parameter,BestCDF,Mean,Std,CoV,Lower,Upper, ...
        Param1_out,Param2_out,Gamma,KS_distance,RMSE_CDF, ...
        'VariableNames',{'Model','Parameter','BestCDF','Mean','Std','CoV', ...
        'Lower','Upper','Param1','Param2','Gamma','KS_distance','RMSE_CDF'});

end

function F = evaluate_cdf(name,x,a,b,p1,p2,gammaVal)

    switch lower(name)

        case 'normal'
            F = normcdf(x,p1,p2);

        case 'lognormal'
            F = logncdf(x,p1,p2);

        case 'weibull'
            F = wblcdf(x,p1,p2);

        case 'bounded sqrt'
            z = (x-a)/(b-a);
            z = min(max(z,0),1);
            F = sqrt(z);

        case 'bounded power'
            z = (x-a)/(b-a);
            z = min(max(z,0),1);
            F = z.^gammaVal;

        otherwise
            error('Unknown CDF type: %s', name);
    end

    F = min(max(F,0),1);

end