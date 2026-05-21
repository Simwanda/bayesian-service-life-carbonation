% COMPARE_PRIOR_POSTERIOR_DISTRIBUTION_TYPES_MULTIPLICATIVE
%
% Compares prior distribution types from Table 2 with fitted posterior
% distribution types from the multiplicative loop.

clear; clc;

outDir = fullfile(pwd,'figs_multiplicative_loop');

Tfib_post = readtable(fullfile(outDir,'combined_fit_summary_fib_multiplicative.csv'));
Tmal_post = readtable(fullfile(outDir,'combined_fit_summary_malami_multiplicative.csv'));

priors_fib = table2_priors_carbonation('fib');
priors_mal = table2_priors_carbonation('malami');

Tfib_prior = priors_to_table(priors_fib,'fib');
Tmal_prior = priors_to_table(priors_mal,'malami');

Tcompare_fib = build_compare_table(Tfib_prior,Tfib_post,'fib');
Tcompare_mal = build_compare_table(Tmal_prior,Tmal_post,'malami');

Tcompare_all = [Tcompare_fib; Tcompare_mal];

disp('Prior vs posterior distribution comparison: fib')
disp(Tcompare_fib)

disp('Prior vs posterior distribution comparison: Malami')
disp(Tcompare_mal)

writetable(Tcompare_fib,fullfile(outDir,'prior_vs_posterior_distribution_types_fib_multiplicative.csv'));
writetable(Tcompare_mal,fullfile(outDir,'prior_vs_posterior_distribution_types_malami_multiplicative.csv'));
writetable(Tcompare_all,fullfile(outDir,'prior_vs_posterior_distribution_types_ALL_multiplicative.csv'));

fprintf('\nSaved prior-vs-posterior comparison tables in:\n  %s\n', outDir);

% ========================================================================
% Local functions
% ========================================================================

function T = priors_to_table(priors,modelName)

    if isfield(priors,'rv')
        R = priors.rv;
    elseif isfield(priors,'Variables')
        R = priors.Variables;
    else
        error('Unknown prior structure. Expected priors.rv or priors.Variables.');
    end

    n = numel(R);

    Model = cell(n,1);
    Parameter = cell(n,1);
    PriorPDF = cell(n,1);
    PriorMean = nan(n,1);
    PriorStd = nan(n,1);
    PriorCoV = nan(n,1);

    for i = 1:n
        Model{i} = modelName;
        Parameter{i} = R(i).Name;
        PriorPDF{i} = R(i).Type;
        PriorMean(i) = R(i).Mean;
        PriorStd(i) = R(i).Std;
        PriorCoV(i) = R(i).COV;
    end

    T = table(Model,Parameter,PriorPDF,PriorMean,PriorStd,PriorCoV);
end

function Tcompare = build_compare_table(Tprior,Tpost,modelName)

    n = height(Tpost);

    Model = cell(n,1);
    Parameter = cell(n,1);

    PriorPDF = cell(n,1);
    PriorMean = nan(n,1);
    PriorStd = nan(n,1);
    PriorCoV = nan(n,1);

    PosteriorPDF = cell(n,1);
    PosteriorMean = nan(n,1);
    PosteriorStd = nan(n,1);
    PosteriorCoV = nan(n,1);

    PosteriorMedian = nan(n,1);
    PosteriorP2p5 = nan(n,1);
    PosteriorP97p5 = nan(n,1);

    AIC_Normal = nan(n,1);
    AIC_Lognormal = nan(n,1);

    for i = 1:n

        pName = char(Tpost.Parameter{i});

        Model{i} = modelName;
        Parameter{i} = pName;

        PosteriorPDF{i} = char(Tpost.BestPDF{i});
        PosteriorMean(i) = Tpost.Mean(i);
        PosteriorStd(i) = Tpost.Std(i);
        PosteriorCoV(i) = Tpost.CoV(i);
        PosteriorMedian(i) = Tpost.Median(i);
        PosteriorP2p5(i) = Tpost.P2p5(i);
        PosteriorP97p5(i) = Tpost.P97p5(i);

        if ismember('AIC_Normal',Tpost.Properties.VariableNames)
            AIC_Normal(i) = Tpost.AIC_Normal(i);
        end
        if ismember('AIC_Lognormal',Tpost.Properties.VariableNames)
            AIC_Lognormal(i) = Tpost.AIC_Lognormal(i);
        end

        idx = find(strcmp(Tprior.Parameter,pName),1);

        if ~isempty(idx)
            PriorPDF{i} = char(Tprior.PriorPDF{idx});
            PriorMean(i) = Tprior.PriorMean(idx);
            PriorStd(i) = Tprior.PriorStd(idx);
            PriorCoV(i) = Tprior.PriorCoV(idx);
        else
            PriorPDF{i} = 'Derived/not prior';
            PriorMean(i) = NaN;
            PriorStd(i) = NaN;
            PriorCoV(i) = NaN;
        end
    end

    Tcompare = table( ...
        Model,Parameter, ...
        PriorPDF,PriorMean,PriorStd,PriorCoV, ...
        PosteriorPDF,PosteriorMean,PosteriorStd,PosteriorCoV, ...
        PosteriorMedian,PosteriorP2p5,PosteriorP97p5, ...
        AIC_Normal,AIC_Lognormal);
end
