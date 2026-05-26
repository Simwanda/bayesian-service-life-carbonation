% EXPORT_POSTERIOR_SAMPLES_FOR_FREET_25MM
% Exports posterior samples for fib and Malami models to CSV.
%
% These CSV files should be given to Martina for FReET.
% IMPORTANT:
%   - x_c(t_obs) is NOT exported as an input variable.
%   - A and n are NOT exported as input variables.
%   - Only basic posterior variables are exported.
%   - For fib, T is exported as a deterministic constant only for CO2 conversion.

clear; clc;

% -------------------------------------------------------------------------
% Input folder
% -------------------------------------------------------------------------
inDir = fullfile(pwd, 'figs_multiplicative_loop_25mm');
matFile = fullfile(inDir, 'combined_loop_results_multiplicative.mat');

if ~isfile(matFile)
    error('File not found: %s', matFile);
end

S = load(matFile);

% -------------------------------------------------------------------------
% Output folder
% -------------------------------------------------------------------------
outDir = fullfile(inDir, 'posterior_samples_for_FReET');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% -------------------------------------------------------------------------
% Get posterior sample structures
% -------------------------------------------------------------------------
if isfield(S, 'combined_fib')
    fibS = S.combined_fib;
else
    error('combined_fib not found in %s', matFile);
end

if isfield(S, 'combined_mal')
    malS = S.combined_mal;
elseif isfield(S, 'combined_malami')
    malS = S.combined_malami;
else
    error('combined_mal or combined_malami not found in %s', matFile);
end

% -------------------------------------------------------------------------
% fib posterior samples
% -------------------------------------------------------------------------
T_fib_det = 283.25; % K, deterministic temperature for CO2 conversion

fibTable = table();

fibTable.Racc_inv = getcol(fibS, 'Racc_inv');
fibTable.kt       = getcol(fibS, 'kt');
fibTable.eps_t    = getcol(fibS, 'eps_t');
fibTable.tc       = getcol(fibS, 'tc');
fibTable.bc       = getcol(fibS, 'bc');
fibTable.CO2ppm   = getcol(fibS, 'CO2ppm');
fibTable.RH       = getcol(fibS, 'RH');
fibTable.tw       = getcol(fibS, 'tw');
fibTable.bw       = getcol(fibS, 'bw');
fibTable.pSR      = getcol(fibS, 'pSR');
fibTable.theta_xc = getcol(fibS, 'theta_xc');

% Deterministic temperature column for fib CO2 conversion
fibTable.T_K      = T_fib_det * ones(height(fibTable), 1);

fibFile = fullfile(outDir, 'posterior_samples_fib_25mm_for_FReET.csv');
writetable(fibTable, fibFile);

% -------------------------------------------------------------------------
% Malami posterior samples
% -------------------------------------------------------------------------
malTable = table();

malTable.W_C      = getcol(malS, 'W_C');
malTable.tc       = getcol(malS, 'tc');
malTable.bc       = getcol(malS, 'bc');
malTable.CO2ppm   = getcol(malS, 'CO2ppm');
malTable.RH       = getcol(malS, 'RH');
malTable.tw       = getcol(malS, 'tw');
malTable.bw       = getcol(malS, 'bw');
malTable.pSR      = getcol(malS, 'pSR');
malTable.T_K      = getcol(malS, 'T');
malTable.theta_xc = getcol(malS, 'theta_xc');

malFile = fullfile(outDir, 'posterior_samples_malami_25mm_for_FReET.csv');
writetable(malTable, malFile);

% -------------------------------------------------------------------------
% Also export quick summary for checking
% -------------------------------------------------------------------------
summaryFib = make_summary(fibTable);
summaryMal = make_summary(malTable);

writetable(summaryFib, fullfile(outDir, 'posterior_samples_fib_25mm_summary.csv'));
writetable(summaryMal, fullfile(outDir, 'posterior_samples_malami_25mm_summary.csv'));

fprintf('\nPosterior samples exported successfully:\n');
fprintf('  %s\n', fibFile);
fprintf('  %s\n', malFile);
fprintf('\nSummary files also saved in:\n');
fprintf('  %s\n', outDir);

% =========================================================================
% Local functions
% =========================================================================
function x = getcol(S, name)
    if ~isfield(S, name)
        error('Field "%s" not found in posterior sample structure.', name);
    end
    x = S.(name);
    x = x(:);
end

function Tout = make_summary(T)
    names = T.Properties.VariableNames(:);
    n = numel(names);

    Mean = zeros(n,1);
    Std  = zeros(n,1);
    CoV  = zeros(n,1);
    Q025 = zeros(n,1);
    Q975 = zeros(n,1);

    for i = 1:n
        x = T.(names{i});
        Mean(i) = mean(x, 'omitnan');
        Std(i)  = std(x, 0, 'omitnan');
        CoV(i)  = Std(i) / abs(Mean(i));
        q = quantile(x, [0.025 0.975]);
        Q025(i) = q(1);
        Q975(i) = q(2);
    end

    Parameter = names;
    Tout = table(Parameter, Mean, Std, CoV, Q025, Q975);
end