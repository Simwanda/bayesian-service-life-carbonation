% CONVERT_DELTA_XC_TO_THETA_XC
%
% Converts additive discrepancy delta_xc [mm] back to an equivalent
% multiplicative model uncertainty factor theta_xc.
%
% Additive model:
%   x_c_add(t) = x_c_base(t) + delta_xc
%
% Equivalent multiplicative model:
%   x_c_add(t) = theta_xc_eq(t) * x_c_base(t)
%
% Therefore:
%   theta_xc_eq(t) = 1 + delta_xc / x_c_base(t)
%
% The conversion is done at t_obs = 42 years.

clear; clc; close all;

outDir = fullfile(pwd,'figs');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

t_obs = 42;  % years

% ------------------------------------------------------------
% Load additive posterior results
% ------------------------------------------------------------
fibFile = 'posterior_table2_fib_additive.mat';
malFile = 'posterior_table2_malami_additive.mat';

if ~isfile(fibFile)
    error('File not found: %s', fibFile);
end

if ~isfile(malFile)
    error('File not found: %s', malFile);
end

Sfib = load(fibFile);
Smal = load(malFile);

post_fib = Sfib.post_fib;
post_mal = Smal.post_malami;

% ------------------------------------------------------------
% Convert fib delta_xc to equivalent theta_xc
% ------------------------------------------------------------
theta_fib = convert_one_model('fib', post_fib, t_obs);

% ------------------------------------------------------------
% Convert Malami delta_xc to equivalent theta_xc
% ------------------------------------------------------------
theta_mal = convert_one_model('malami', post_mal, t_obs);

% ------------------------------------------------------------
% Summaries
% ------------------------------------------------------------
T_fib = summarize_theta(theta_fib,'fib');
T_mal = summarize_theta(theta_mal,'malami');

T_all = [T_fib; T_mal];

disp(T_all)

% ------------------------------------------------------------
% Save equivalent theta samples and summary
% ------------------------------------------------------------
save(fullfile(outDir,'theta_xc_equivalent_from_delta.mat'), ...
    'theta_fib','theta_mal','T_all','t_obs');

writetable(T_all,fullfile(outDir,'theta_xc_equivalent_from_delta_summary.csv'));

% ------------------------------------------------------------
% Plot histograms and fitted PDFs
% ------------------------------------------------------------
plot_theta_fit(theta_fib,'fib',outDir);
plot_theta_fit(theta_mal,'malami',outDir);

fprintf('\nSaved:\n');
fprintf('  %s\n', fullfile(outDir,'theta_xc_equivalent_from_delta.mat'));
fprintf('  %s\n', fullfile(outDir,'theta_xc_equivalent_from_delta_summary.csv'));
fprintf('  %s\n', fullfile(outDir,'theta_xc_equivalent_from_delta_fib.png'));
fprintf('  %s\n', fullfile(outDir,'theta_xc_equivalent_from_delta_malami.png'));

% ========================================================================
% Local functions
% ========================================================================

function theta_eq = convert_one_model(modelName, post, t_obs)

    samples = post.samples;

    if ~isfield(samples,'delta_xc')
        error('The posterior samples do not contain delta_xc.');
    end

    delta = samples.delta_xc(:);

    % Create base samples with additive discrepancy removed
    samples_base = samples;
    samples_base.delta_xc = zeros(size(delta));

    % If theta_xc exists, keep it equal to 1 for the additive conversion
    if isfield(samples_base,'theta_xc')
        samples_base.theta_xc = ones(size(delta));
    end

    % Base carbonation prediction without additive discrepancy
    [x_base,~] = carbonation_xc_table2(modelName, samples_base, t_obs, post.priors);
    x_base = x_base(:);

    % Avoid division by zero
    idx = isfinite(x_base) & isfinite(delta) & abs(x_base) > eps;

    theta_eq = NaN(size(delta));
    theta_eq(idx) = 1 + delta(idx)./x_base(idx);

    theta_eq = theta_eq(isfinite(theta_eq));

end

function T = summarize_theta(theta_eq,modelName)

    theta_eq = theta_eq(:);
    theta_eq = theta_eq(isfinite(theta_eq));

    Mean = mean(theta_eq);
    Std = std(theta_eq);
    CoV = Std / abs(Mean);
    Median = median(theta_eq);
    q = quantile(theta_eq,[0.025 0.975]);

    % Fit Normal
    pdN = fitdist(theta_eq,'Normal');
    pN = pdf(pdN,theta_eq);
    pN(pN < realmin) = realmin;
    logLN = sum(log(pN));
    AIC_Normal = 2*2 - 2*logLN;

    % Fit Lognormal if positive
    if all(theta_eq > 0)
        pdL = fitdist(theta_eq,'Lognormal');
        pL = pdf(pdL,theta_eq);
        pL(pL < realmin) = realmin;
        logLL = sum(log(pL));
        AIC_Lognormal = 2*2 - 2*logLL;
    else
        AIC_Lognormal = NaN;
    end

    if isnan(AIC_Lognormal)
        BestPDF = {'Normal'};
    elseif AIC_Lognormal < AIC_Normal
        BestPDF = {'Lognormal'};
    else
        BestPDF = {'Normal'};
    end

    Model = {modelName};
    Parameter = {'theta_xc_equiv_from_delta'};
    P2p5 = q(1);
    P97p5 = q(2);

    T = table(Model,Parameter,BestPDF,Mean,Std,CoV,Median,P2p5,P97p5, ...
        AIC_Normal,AIC_Lognormal);

end

function plot_theta_fit(theta_eq,modelName,outDir)

    theta_eq = theta_eq(:);
    theta_eq = theta_eq(isfinite(theta_eq));

    fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
    hold on; box on; grid on;

    histogram(theta_eq,40,'Normalization','pdf', ...
        'FaceAlpha',0.35, ...
        'EdgeColor','none');

    xmin = min(theta_eq);
    xmax = max(theta_eq);
    pad = 0.05*(xmax-xmin);
    xx = linspace(xmin-pad,xmax+pad,400);

    legendEntries = {'Equivalent \theta_{xc} samples'};

    % Normal fit
    pdN = fitdist(theta_eq,'Normal');
    plot(xx,pdf(pdN,xx),'k-','LineWidth',1.8);
    legendEntries{end+1} = 'Normal fit';

    % Lognormal fit
    if all(theta_eq > 0)
        pdL = fitdist(theta_eq,'Lognormal');

        yyL = nan(size(xx));
        idxPos = xx > 0;
        yyL(idxPos) = pdf(pdL,xx(idxPos));

        plot(xx,yyL,'k--','LineWidth',1.8);
        legendEntries{end+1} = 'Lognormal fit';
    end

    xlabel('\theta_{x_c}^{eq}');
    ylabel('Probability density');
    title(sprintf('%s: equivalent multiplicative factor from \\delta_{xc}',modelName), ...
        'Interpreter','tex');

    legend(legendEntries,'Location','best');

    mu = mean(theta_eq);
    sig = std(theta_eq);
    covx = sig / abs(mu);
    q = quantile(theta_eq,[0.025 0.975]);

    txt = sprintf(['Mean = %.4f\nStd = %.4f\nCoV = %.3f\n', ...
                   '2.5%% = %.4f\n97.5%% = %.4f'], ...
                   mu, sig, covx, q(1), q(2));

    yLim = ylim;
    text(xmin + 0.05*(xmax-xmin), ...
         yLim(2) - 0.08*(yLim(2)-yLim(1)), ...
         txt, ...
         'VerticalAlignment','top', ...
         'BackgroundColor','w', ...
         'EdgeColor',[0.4 0.4 0.4], ...
         'Margin',5);

    pngFile = fullfile(outDir,sprintf('theta_xc_equivalent_from_delta_%s.png',modelName));
    pdfFile = fullfile(outDir,sprintf('theta_xc_equivalent_from_delta_%s.pdf',modelName));

    exportgraphics(fig,pngFile,'Resolution',300);
    exportgraphics(fig,pdfFile,'ContentType','vector');

    close(fig);

end