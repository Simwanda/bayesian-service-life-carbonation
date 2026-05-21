% CHECK_ADDITIVE_DELTA_XC_SETUP
% Quick check that delta_xc is defined in the active prior file.
clear; clc;
which table2_priors_carbonation -all
pFib = table2_priors_carbonation('fib');
pMal = table2_priors_carbonation('malami');
disp('fib variables:'); disp(pFib.names(:));
disp('Malami variables:'); disp(pMal.names(:));
assert(any(strcmp(pFib.names,'delta_xc')), 'delta_xc missing for fib');
assert(any(strcmp(pMal.names,'delta_xc')), 'delta_xc missing for Malami');
disp('OK: delta_xc exists in both fib and Malami priors.');
