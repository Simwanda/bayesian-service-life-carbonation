% EXPORT_EMPIRICAL_CDF_PSR_TC_FOR_FREET
%
% Exports x and F(x) values for empirical posterior distributions
% of pSR and tc for use in FReET.
%
% Input:
%   figs_multiplicative_loop/combined_loop_results_multiplicative.mat
%
% Output:
%   figs_multiplicative_loop/empirical_CDF_for_FReET/
%       fib_pSR_CDF_for_FReET.csv
%       fib_tc_CDF_for_FReET.csv
%       malami_pSR_CDF_for_FReET.csv
%       malami_tc_CDF_for_FReET.csv
%
% Each CSV contains:
%   x      = variable value
%   Fx     = empirical cumulative probability

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
cdfDir = fullfile(outDir,'empirical_CDF_for_FReET');

if ~exist(cdfDir,'dir')
    mkdir(cdfDir);
end

% ------------------------------------------------------------
% Export variables
% ------------------------------------------------------------
export_one_model(combined_fib,'fib','pSR',cdfDir);
export_one_model(combined_fib,'fib','tc',cdfDir);

export_one_model(combined_mal,'malami','pSR',cdfDir);
export_one_model(combined_mal,'malami','tc',cdfDir);

fprintf('\nEmpirical CDF files for FReET saved in:\n');
fprintf('  %s\n', cdfDir);

% ========================================================================
% Local function
% ========================================================================

function export_one_model(S,modelName,varName,cdfDir)

    if ~isfield(S,varName)
        warning('%s does not contain variable %s. Skipping.', modelName, varName);
        return;
    end

    x = S.(varName);
    x = x(:);
    x = x(isfinite(x));

    if numel(x) < 10
        warning('%s %s has too few samples. Skipping.', modelName, varName);
        return;
    end

    % Sort posterior samples
    x = sort(x);

    % Empirical CDF values
    n = numel(x);
    Fx = (1:n)' ./ n;

    % Remove duplicate x-values because some software requires unique x
    [xUnique, ia] = unique(x,'last');
    FxUnique = Fx(ia);

    % Force the CDF to start close to 0 and end at 1
    FxUnique(1) = max(FxUnique(1), 1/n);
    FxUnique(end) = 1.0;

    % Create table
    T = table(xUnique(:),FxUnique(:), ...
        'VariableNames',{'x','Fx'});

    % Save CSV
    csvFile = fullfile(cdfDir, ...
        sprintf('%s_%s_CDF_for_FReET.csv',modelName,varName));

    writetable(T,csvFile);

    % Also save Excel file if supported
    xlsxFile = fullfile(cdfDir, ...
        sprintf('%s_%s_CDF_for_FReET.xlsx',modelName,varName));

    try
        writetable(T,xlsxFile);
    catch
        warning('Could not write Excel file. CSV file was saved.');
    end

    fprintf('Saved %s %s empirical CDF: %s\n', modelName, varName, csvFile);

end