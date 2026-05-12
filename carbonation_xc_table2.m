function [xc, aux] = carbonation_xc_table2(modelName, X, t_years, priors)
%CARBONATION_XC_TABLE2 Carbonation-depth forward model used in Bayesian updating.
%   [xc, aux] = carbonation_xc_table2(modelName, X, t_years, priors)
%   returns x_c in mm for each row/sample of X.
%
%   For modelName = 'fib', the code uses the reduced fib Bulletin 34 rate form.
%
%   For modelName = 'malami', the code now uses the Malami forward model exactly
%   in the closed-form implementation requested by the authors/reviewer response:
%       D_c = sqrt(2*k_RH*k_T*k_c*CO2*R_NAC0_inv*t_days)*Wwet*1000
%   with R_NAC0_inv = DeCO2/aCO2, hydration-scaled CO2 binding capacity,
%   Arrhenius temperature correction, humidity correction, curing correction, and
%   fib time-of-wetness exponent. The implemented equations correspond to the
%   Python function supplied by the user, translated to MATLAB, while still reading
%   Table 2 random variables from X/priors.

modelName = lower(char(modelName));
if nargin < 4 || isempty(priors)
    priors = table2_priors_carbonation(modelName);
end

% Make all fields column vectors of common length.
X = complete_with_deterministics(X, priors);
N = infer_N(X);
t_years = reshape(t_years,1,[]); % 1 x Nt

% Common exponent from fib weathering function, stored for diagnostics.
pSR = col(X.pSR,N);
tw  = col(X.tw,N);
bw  = max(col(X.bw,N), 1e-6);
n = 0.5 .* (pSR .* tw) .^ bw;
n(~isfinite(n)) = priors.n_prior_mean;
n = max(n, 0);

switch modelName
    case 'fib'
        A0mean = priors.A_prior_mean_fib_mm_y;
        ref = priors.meanStruct;

        RH = min(max(col(X.RH,N)./100, 0.01), 0.99);
        RHref = ref.RH/100;
        ke = ((1 - RH.^5) ./ (1 - RHref.^5)).^2.5;
        ke_ref = 1.0;

        tc = max(col(X.tc,N), 1e-6);
        bc = col(X.bc,N);
        kc = (tc./7).^bc;
        kc_ref = (ref.tc/7)^ref.bc;

        Rnat = max(col(X.kt,N).*col(X.Racc_inv,N) + col(X.eps_t,N), realmin);
        Rnat_ref = ref.kt*ref.Racc_inv + ref.eps_t;

        CO2 = max(col(X.CO2ppm,N), realmin);
        CO2_ref = ref.CO2ppm;

        scale = sqrt((ke./ke_ref) .* (kc./kc_ref) .* (Rnat./Rnat_ref) .* (CO2./CO2_ref));
        A = A0mean .* scale;

        theta = max(col(X.theta_xc,N), realmin);
        xc = (theta .* A) .* (t_years .^ (0.5 - n));

        aux.A = A;
        aux.n = n;
        aux.theta_xc = theta;
        aux.Rnat = Rnat;
        aux.kc = kc;
        aux.ke = ke;

    case 'malami'
        [Dc, auxM] = malami_model_forward_table2(X, t_years, priors, N);
        theta = max(col(X.theta_xc,N), realmin);
        xc = theta .* Dc;

        aux = auxM;
        aux.n = n;
        aux.theta_xc = theta;
        % Diagnostic A in units consistent with the supplied Python function:
        % D_c / t_days^(0.5-n), not used in the likelihood.
        t_days = t_years .* 365;
        aux.A = Dc ./ (t_days .^ (0.5 - n));

    otherwise
        error('Unknown modelName: use ''fib'' or ''malami''.');
end
end

function [D_c, aux] = malami_model_forward_table2(X, ty, priors, N)
%MALAMI_MODEL_FORWARD_TABLE2 MATLAB translation of the requested Python model.
%
% Python reference supplied by user:
% def malami_model(ty, W_C, C, Ps, T, t_c, k_urb, RH, ToW, Psr, phi_cl, CaO, k):
%     t = ty * 365
%     RH_ref = 0.65
%     k_RH = (RH / RH_ref) ** 2.6
%     W = W_C * C
%     Ceff = (C-Ps) + k * Ps
%     W_Ceff = W / Ceff
%     R = 8.314
%     Tref = 298.15
%     E_a = (24 - 8.7 * W_Ceff) * 1000
%     k_T = exp((E_a/R)*((1/Tref)-(1/T)))
%     b_c = -0.567
%     k_c = (t_c/7)**b_c
%     CO2_a = (44*424*101325)/(R*T*1e9)
%     CO2 = k_urb*CO2_a
%     alpha_inf = (1.031*W_Ceff)/(0.194+W_Ceff)
%     alpha_H = alpha_inf*(t/(2+t))
%     a_CO2 = 0.75*phi_cl*C*CaO*alpha_H*(44/56)
%     DeCO2 = 6.1e-6*porosity_term**3*(1-RH)**2.2
%     R_NAC0_inv = (DeCO2/a_CO2)*24*3600
%     n = 0.5*(Psr*ToW)**0.446
%     Wwet = (28/t)**n
%     D_c = sqrt(2*k_RH*k_T*k_c*CO2*R_NAC0_inv*t)*Wwet*1000

% Time in days, as in the supplied Python implementation.
t = ty .* 365; % 1 x Nt

% Table 2 variables. RH is stored as % in Table 2, converted to fraction here.
W_C    = max(col(X.W_C,N), 1e-6);
C      = max(col(X.C,N), realmin);
Ps     = max(col(X.P,N), 0);
T      = max(col(X.T,N), 1);
t_c    = max(col(X.tc,N), realmin);
k_urb  = max(col(X.kurb,N), realmin);
RH     = min(max(col(X.RH,N)./100, 1e-6), 0.999999);
ToW    = max(col(X.tw,N), 0);
Psr    = max(col(X.pSR,N), 0);
phi_cl = max(col(X.phi_cl,N), realmin);
CaO    = max(col(X.CaO,N), realmin);
k      = col(X.k,N);

% Keep the user-supplied constants, but allow Table 2 random bc/bw/CO2ppm to be
% used when present. At prior means this reduces to b_c=-0.567 and b_w=0.446.
if isfield(X,'bc'); b_c = col(X.bc,N); else; b_c = -0.567*ones(N,1); end
if isfield(X,'bw'); b_w = col(X.bw,N); else; b_w = 0.446*ones(N,1); end
if isfield(X,'CO2ppm'); CO2ppm = col(X.CO2ppm,N); else; CO2ppm = 424*ones(N,1); end

RH_ref = 0.65;
k_RH = (RH ./ RH_ref) .^ 2.6;

W = W_C .* C;
Ceff = (C - Ps) + k .* Ps;
Ceff = max(Ceff, realmin);
W_Ceff = W ./ Ceff;

R = 8.314;      % J/(mol.K)
Tref = 298.15;  % K
E_a = (24 - 8.7 .* W_Ceff) .* 1000; % J/mol
k_T = exp((E_a ./ R) .* ((1 ./ Tref) - (1 ./ T)));

k_c = (t_c ./ 7) .^ b_c;

M_CO2 = 44;
P_atm = 101325;
CO2_a = (M_CO2 .* CO2ppm .* P_atm) ./ (R .* T .* 1e9);
CO2 = k_urb .* CO2_a;

alpha_inf = (1.031 .* W_Ceff) ./ (0.194 + W_Ceff);
alpha_H = alpha_inf .* (t ./ (2 + t));

M_CaO = 56;
a_CO2 = 0.75 .* phi_cl .* C .* CaO .* alpha_H .* (M_CO2 ./ M_CaO);
a_CO2 = max(a_CO2, realmin);

rhoC = priors.det.rhoC;
rhoW = priors.det.rhoW;
porTerm = (((W - 0.267 .* Ceff) ./ 1000) ./ ((Ceff ./ rhoC) + (W ./ rhoW)));
porTerm = max(porTerm, realmin);
DeCO2 = 6.1e-6 .* porTerm .^ 3 .* (1 - RH) .^ 2.2;

R_NAC0 = a_CO2 ./ max(DeCO2, realmin);
R_NAC0_inv = (1 ./ R_NAC0) .* 24 .* 3600;

n = 0.5 .* (Psr .* ToW) .^ b_w;
Wwet = (28 ./ t) .^ n;

D_c = sqrt(2 .* k_RH .* k_T .* k_c .* CO2 .* R_NAC0_inv .* t) .* Wwet .* 1000;
D_c(~isfinite(D_c)) = NaN;

aux.k_RH = k_RH;
aux.k_T = k_T;
aux.k_c = k_c;
aux.CO2 = CO2;
aux.alpha_H = alpha_H;
aux.a_CO2 = a_CO2;
aux.DeCO2 = DeCO2;
aux.R_NAC0 = R_NAC0;
aux.R_NAC0_inv = R_NAC0_inv;
aux.n_malami = n;
aux.Wwet = Wwet;
aux.W_Ceff = W_Ceff;
end

function X = complete_with_deterministics(X, priors)
% Add deterministic Table 2 fields where missing.
if isfield(priors,'det')
    fn = fieldnames(priors.det);
    for i = 1:numel(fn)
        if ~isfield(X,fn{i}); X.(fn{i}) = priors.det.(fn{i}); end
    end
end
end

function N = infer_N(X)
fn = fieldnames(X); N = 1;
for i = 1:numel(fn)
    N = max(N, numel(X.(fn{i})));
end
end

function v = col(x,N)
if isscalar(x)
    v = repmat(x,N,1);
else
    v = x(:);
    if numel(v) ~= N
        error('Input vector length mismatch.');
    end
end
end
