function myInput = buildInput_carbonation_posterior_table2(fit)
%BUILDINPUT_CARBONATION_POSTERIOR_TABLE2 Create UQLab Input from posterior marginals.
%   Requires UQLab already on path and initialized with uqlab.

names = fit.updateNames;
M = struct([]);
for i = 1:numel(names)
    f = fit.marginals.(names{i});
    M(i).Name = f.Name;
    M(i).Type = f.UQType;
    M(i).Parameters = f.UQParameters;
end
InputOpts.Marginals = M;
InputOpts.Name = sprintf('Posterior_Table2_%s', fit.modelName);
myInput = uq_createInput(InputOpts);
end
