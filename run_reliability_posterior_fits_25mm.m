% RUN_RELIABILITY_POSTERIOR_FITS_25MM
%
% Reliability analysis using posterior distribution fits from the
% sigma_epsilon = 2.5 mm Bayesian updating results.
%
% Outputs:
%   figs_reliability_posterior_25mm/Figure9_reliability_beta_time.pdf/png
%   figs_reliability_posterior_25mm/Figure10_FORM_importance_basic.pdf/png
%   figs_reliability_posterior_25mm/Figure11_FORM_importance_An.pdf/png
%   figs_reliability_posterior_25mm/Table5_reliability_results_42years.csv
%   figs_reliability_posterior_25mm/reliability_curves.csv
%   figs_reliability_posterior_25mm/FORM_importance_basic_42years.csv
%   figs_reliability_posterior_25mm/FORM_importance_An_42years.csv
%
% IMPORTANT:
%   1) This script uses the posterior distribution FITS for sigma_epsilon = 2.5 mm.
%   2) It includes the fitted bounded-power CDFs for t_c and p_SR.
%   3) sigma_epsilon is NOT added again in the reliability analysis.
%      It has already been used in the Bayesian likelihood to obtain the posterior distributions.
%   4) theta_xc is included in x_c(t).
%   5) Limit state: g(t) = C_r - x_c(t). Failure/depassivation if g(t) <= 0.
%
% Required in the MATLAB path:
%   - carbonation_xc_table2.m
%   - table2_priors_carbonation.m, if carbonation_xc_table2 requires priors
%
% Author: generated for Lenganji Simwanda

clear; clc; close all;

% Use LaTeX rendering in all figure text and tick labels.
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

rng(20260528);

% -------------------------------------------------------------------------
% User settings
% -------------------------------------------------------------------------
outDir = fullfile(pwd,'figs_reliability_posterior_25mm');
if ~exist(outDir,'dir'); mkdir(outDir); end

% Monte Carlo sample size. Increase to 1e6 for final paper if time permits.
Nmc = 300000;

% Time grid for reliability curves
tGrid = 0.5:1:70;

% Inspection time for Table 5 and FORM sensitivity factors
t_assess = 42;

% Concrete cover scenarios
Cr_mean = 20.5;       % mm
Cr_cov  = 0.30;       % only used in Case II
Cr_det  = 20.5;       % mm, Case I

% Target reliability level
beta_target = 1.5;

% -------------------------------------------------------------------------
% Define posterior random variables
% -------------------------------------------------------------------------
rv_fib = define_rvs_fib();
rv_mal = define_rvs_malami();

% Add deterministic constants that are required by carbonation_xc_table2
const_fib = struct();
const_fib.T = 283.25;     % K, deterministic temperature for CO2 conversion in fib
const_fib.kurb = 1.0;

const_mal = struct();
const_mal.C      = 375;
const_mal.P      = 0;
const_mal.k      = 1;
const_mal.phi_cl = 0.975;
const_mal.CaO    = 0.68;
const_mal.kurb   = 1.0;
const_mal.rhoC   = 3120;
const_mal.rhoW   = 1000;

% For aggregate A-n-theta reliability analysis
rv_fib_An = define_rvs_fib_An();
rv_mal_An = define_rvs_malami_An();

% -------------------------------------------------------------------------
% Monte Carlo reliability curves
% -------------------------------------------------------------------------
fprintf('\nRunning Monte Carlo reliability analysis with N = %d samples...\n',Nmc);

[Xfib_mc, names_fib] = sample_rvs(rv_fib,Nmc);
Xfib_mc = add_constants(Xfib_mc,const_fib,Nmc);

[Xmal_mc, names_mal] = sample_rvs(rv_mal,Nmc);
Xmal_mc = add_constants(Xmal_mc,const_mal,Nmc);

Cr_caseI = Cr_det * ones(Nmc,1);
Cr_caseII = sample_lognormal_mean_cov(Cr_mean,Cr_cov,Nmc);

[beta_fib_caseI_mcs, pf_fib_caseI_mcs] = mcs_curve('fib',Xfib_mc,Cr_caseI,tGrid);
[beta_fib_caseII_mcs,pf_fib_caseII_mcs] = mcs_curve('fib',Xfib_mc,Cr_caseII,tGrid);

[beta_mal_caseI_mcs, pf_mal_caseI_mcs] = mcs_curve('malami',Xmal_mc,Cr_caseI,tGrid);
[beta_mal_caseII_mcs,pf_mal_caseII_mcs] = mcs_curve('malami',Xmal_mc,Cr_caseII,tGrid);

% -------------------------------------------------------------------------
% FORM reliability curves
% -------------------------------------------------------------------------
fprintf('\nRunning FORM reliability curves...\n');

beta_fib_caseI_form  = zeros(size(tGrid));
beta_fib_caseII_form = zeros(size(tGrid));
beta_mal_caseI_form  = zeros(size(tGrid));
beta_mal_caseII_form = zeros(size(tGrid));

for i = 1:numel(tGrid)
    t = tGrid(i);

    beta_fib_caseI_form(i)  = form_beta_basic('fib',rv_fib,const_fib,t,'caseI',Cr_det,Cr_mean,Cr_cov);
    beta_fib_caseII_form(i) = form_beta_basic('fib',rv_fib,const_fib,t,'caseII',Cr_det,Cr_mean,Cr_cov);

    beta_mal_caseI_form(i)  = form_beta_basic('malami',rv_mal,const_mal,t,'caseI',Cr_det,Cr_mean,Cr_cov);
    beta_mal_caseII_form(i) = form_beta_basic('malami',rv_mal,const_mal,t,'caseII',Cr_det,Cr_mean,Cr_cov);

    if mod(i,10)==0
        fprintf('  FORM completed for %d/%d time points\n',i,numel(tGrid));
    end
end

% -------------------------------------------------------------------------
% Reliability summary at t = 42 years
% -------------------------------------------------------------------------
idx42 = find(abs(tGrid - t_assess)==min(abs(tGrid - t_assess)),1);

T42 = table( ...
    {'fib Bulletin 34 - Case I';'fib Bulletin 34 - Case II'; ...
     'Malami - Case I';'Malami - Case II'}, ...
    [pf_fib_caseI_mcs(idx42); pf_fib_caseII_mcs(idx42); pf_mal_caseI_mcs(idx42); pf_mal_caseII_mcs(idx42)], ...
    [beta_fib_caseI_mcs(idx42); beta_fib_caseII_mcs(idx42); beta_mal_caseI_mcs(idx42); beta_mal_caseII_mcs(idx42)], ...
    [normcdf(-beta_fib_caseI_form(idx42)); normcdf(-beta_fib_caseII_form(idx42)); ...
     normcdf(-beta_mal_caseI_form(idx42)); normcdf(-beta_mal_caseII_form(idx42))], ...
    [beta_fib_caseI_form(idx42); beta_fib_caseII_form(idx42); beta_mal_caseI_form(idx42); beta_mal_caseII_form(idx42)], ...
    'VariableNames',{'Model_case','Pf_MCS','Beta_MCS','Pf_FORM','Beta_FORM'});

disp(T42)
writetable(T42,fullfile(outDir,'Table5_reliability_results_42years.csv'));

% Save reliability curves
Tcurves = table(tGrid(:), ...
    pf_fib_caseI_mcs(:),  beta_fib_caseI_mcs(:),  beta_fib_caseI_form(:), ...
    pf_fib_caseII_mcs(:), beta_fib_caseII_mcs(:), beta_fib_caseII_form(:), ...
    pf_mal_caseI_mcs(:),  beta_mal_caseI_mcs(:),  beta_mal_caseI_form(:), ...
    pf_mal_caseII_mcs(:), beta_mal_caseII_mcs(:), beta_mal_caseII_form(:), ...
    'VariableNames',{'t_years', ...
    'pf_fib_caseI_mcs','beta_fib_caseI_mcs','beta_fib_caseI_form', ...
    'pf_fib_caseII_mcs','beta_fib_caseII_mcs','beta_fib_caseII_form', ...
    'pf_malami_caseI_mcs','beta_malami_caseI_mcs','beta_malami_caseI_form', ...
    'pf_malami_caseII_mcs','beta_malami_caseII_mcs','beta_malami_caseII_form'});
writetable(Tcurves,fullfile(outDir,'reliability_curves.csv'));

% -------------------------------------------------------------------------
% Plot Figure 9: time-dependent reliability
% -------------------------------------------------------------------------
plot_reliability_figure(outDir,tGrid,beta_target, ...
    beta_fib_caseI_mcs,beta_fib_caseI_form,beta_fib_caseII_mcs,beta_fib_caseII_form, ...
    beta_mal_caseI_mcs,beta_mal_caseI_form,beta_mal_caseII_mcs,beta_mal_caseII_form);

% -------------------------------------------------------------------------
% FORM sensitivity factors at 42 years: basic variables
% -------------------------------------------------------------------------
fprintf('\nComputing FORM sensitivity factors at t = %.1f years...\n',t_assess);

[~,alpha2_fib_caseI,names_fib_caseI]   = form_beta_basic('fib',rv_fib,const_fib,t_assess,'caseI',Cr_det,Cr_mean,Cr_cov);
[~,alpha2_fib_caseII,names_fib_caseII] = form_beta_basic('fib',rv_fib,const_fib,t_assess,'caseII',Cr_det,Cr_mean,Cr_cov);
[~,alpha2_mal_caseI,names_mal_caseI]   = form_beta_basic('malami',rv_mal,const_mal,t_assess,'caseI',Cr_det,Cr_mean,Cr_cov);
[~,alpha2_mal_caseII,names_mal_caseII] = form_beta_basic('malami',rv_mal,const_mal,t_assess,'caseII',Cr_det,Cr_mean,Cr_cov);

Timp_basic = make_importance_table( ...
    names_fib_caseI,alpha2_fib_caseI,names_fib_caseII,alpha2_fib_caseII, ...
    names_mal_caseI,alpha2_mal_caseI,names_mal_caseII,alpha2_mal_caseII);
writetable(Timp_basic,fullfile(outDir,'FORM_importance_basic_42years.csv'));

plot_importance_basic(outDir, ...
    names_fib_caseI,alpha2_fib_caseI,names_fib_caseII,alpha2_fib_caseII, ...
    names_mal_caseI,alpha2_mal_caseI,names_mal_caseII,alpha2_mal_caseII);

% -------------------------------------------------------------------------
% FORM sensitivity factors at 42 years: aggregate A, n, theta_xc, Cr
% -------------------------------------------------------------------------
[~,alpha2_fibAn_caseI,names_fibAn_caseI]   = form_beta_An(rv_fib_An,t_assess,'caseI',Cr_det,Cr_mean,Cr_cov);
[~,alpha2_fibAn_caseII,names_fibAn_caseII] = form_beta_An(rv_fib_An,t_assess,'caseII',Cr_det,Cr_mean,Cr_cov);
[~,alpha2_malAn_caseI,names_malAn_caseI]   = form_beta_An(rv_mal_An,t_assess,'caseI',Cr_det,Cr_mean,Cr_cov);
[~,alpha2_malAn_caseII,names_malAn_caseII] = form_beta_An(rv_mal_An,t_assess,'caseII',Cr_det,Cr_mean,Cr_cov);

Timp_An = make_importance_table( ...
    names_fibAn_caseI,alpha2_fibAn_caseI,names_fibAn_caseII,alpha2_fibAn_caseII, ...
    names_malAn_caseI,alpha2_malAn_caseI,names_malAn_caseII,alpha2_malAn_caseII);
writetable(Timp_An,fullfile(outDir,'FORM_importance_An_42years.csv'));

plot_importance_An(outDir, ...
    names_fibAn_caseI,alpha2_fibAn_caseI,names_fibAn_caseII,alpha2_fibAn_caseII, ...
    names_malAn_caseI,alpha2_malAn_caseI,names_malAn_caseII,alpha2_malAn_caseII);

fprintf('\nDONE. Results saved in:\n  %s\n',outDir);

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function rv = define_rvs_fib()
    % Posterior distribution fits, sigma_epsilon = 2.5 mm.
    rv = [];
    rv = add_rv(rv,'Racc_inv','lognormal',1.177241e-10,0.430867);
    rv = add_rv(rv,'kt','lognormal',1.324822,0.225232);
    rv = add_rv(rv,'eps_t','normal',9.966963e-12,0.152217);
    rv = add_rv(rv,'tc','bounded_power',1.829238,0.315874,1.0,3.0,0.7209);
    rv = add_rv(rv,'bc','normal',-0.5690245,0.041911);
    rv = add_rv(rv,'CO2ppm','lognormal',392.0933,0.008702);
    rv = add_rv(rv,'RH','normal',75.38409,0.034762);
    rv = add_rv(rv,'tw','normal',0.1550595,0.118126);
    rv = add_rv(rv,'bw','normal',0.5183326,0.287560);
    rv = add_rv(rv,'pSR','bounded_power',0.1355876,0.430757,0.05,0.25,0.7597);
    rv = add_rv(rv,'theta_xc','normal',1.086561,0.175945);
end

function rv = define_rvs_malami()
    rv = [];
    rv = add_rv(rv,'W_C','normal',0.567922,0.0930956);
    rv = add_rv(rv,'tc','bounded_power',1.899948,0.297335,1.0,3.0,0.8320);
    rv = add_rv(rv,'bc','normal',-0.567731,0.0415413);
    rv = add_rv(rv,'CO2ppm','lognormal',392.260152,0.0087984);
    rv = add_rv(rv,'RH','normal',75.277292,0.0348594);
    rv = add_rv(rv,'tw','normal',0.153987,0.119057);
    rv = add_rv(rv,'bw','normal',0.523812,0.271112);
    rv = add_rv(rv,'pSR','bounded_power',0.133728,0.430707,0.05,0.25,0.7395);
    rv = add_rv(rv,'T','lognormal',283.198264,0.00188433);
    rv = add_rv(rv,'theta_xc','normal',1.072589,0.180409);
end

function rv = define_rvs_fib_An()
    rv = [];
    rv = add_rv(rv,'A','lognormal',2.531550,0.252839);
    rv = add_rv(rv,'n','lognormal',0.07747914,0.648510);
    rv = add_rv(rv,'theta_xc','normal',1.086561,0.175945);
end

function rv = define_rvs_malami_An()
    rv = [];
    rv = add_rv(rv,'A','lognormal',0.213724,0.264053);
    rv = add_rv(rv,'n','lognormal',0.073877,0.610160);
    rv = add_rv(rv,'theta_xc','normal',1.072589,0.180409);
end

function rv = add_rv(rv,name,type,meanVal,covVal,a,b,gammaVal)
    r.name = name;
    r.type = type;
    r.mean = meanVal;
    r.cov  = covVal;
    r.std  = abs(meanVal)*covVal;
    if nargin >= 6
        r.a = a; r.b = b; r.gamma = gammaVal;
    else
        r.a = NaN; r.b = NaN; r.gamma = NaN;
    end
    rv = [rv; r];
end

function [X,names] = sample_rvs(rv,N)
    X = struct();
    names = cell(numel(rv),1);
    for i = 1:numel(rv)
        names{i} = rv(i).name;
        X.(rv(i).name) = sample_one_rv(rv(i),N);
    end
end

function x = sample_one_rv(r,N)
    switch lower(r.type)
        case 'normal'
            x = r.mean + r.std.*randn(N,1);
        case 'lognormal'
            [mu_ln,sig_ln] = lognormal_mu_sigma(r.mean,r.cov);
            x = exp(mu_ln + sig_ln.*randn(N,1));
        case 'bounded_power'
            u = rand(N,1);
            x = r.a + (r.b-r.a).*u.^(1/r.gamma);
        otherwise
            error('Unknown RV type: %s',r.type);
    end
end

function X = add_constants(X,const,N)
    fn = fieldnames(const);
    for i = 1:numel(fn)
        X.(fn{i}) = const.(fn{i})*ones(N,1);
    end
end

function x = sample_lognormal_mean_cov(meanVal,covVal,N)
    [mu_ln,sig_ln] = lognormal_mu_sigma(meanVal,covVal);
    x = exp(mu_ln + sig_ln.*randn(N,1));
end

function [mu_ln,sig_ln] = lognormal_mu_sigma(meanVal,covVal)
    sig_ln = sqrt(log(1+covVal.^2));
    mu_ln = log(meanVal) - 0.5*sig_ln.^2;
end

function [beta,pf] = mcs_curve(model,X,Cr,tGrid)
    beta = zeros(size(tGrid));
    pf   = zeros(size(tGrid));
    for i = 1:numel(tGrid)
        xc = eval_xc_model(model,X,tGrid(i));
        g = Cr - xc(:);
        pf(i) = mean(g <= 0);
        pf(i) = min(max(pf(i),1e-12),1-1e-12);
        beta(i) = -norminv(pf(i));
    end
end

function xc = eval_xc_model(model,X,t)
    % Uses the same forward function as the Bayesian updating workflow.
    try
        priors = table2_priors_carbonation();
    catch
        priors = [];
    end

    try
        [xc,~] = carbonation_xc_table2(model,X,t,priors);
    catch ME
        fprintf('\ncarbonation_xc_table2 failed. Check that this script is in the repository root.\n');
        rethrow(ME);
    end

    xc = xc(:);
end

function [beta,alpha2,names] = form_beta_basic(model,rv,const,t,caseName,Cr_det,Cr_mean,Cr_cov)
    rv_form = rv;
    if strcmpi(caseName,'caseII')
        rv_form = add_rv(rv_form,'Cr','lognormal',Cr_mean,Cr_cov);
    end

    names = {rv_form.name}';
    nvar = numel(rv_form);

    gfun = @(u) gfun_basic(u,model,rv_form,const,t,caseName,Cr_det);

    [beta,uStar,gradStar] = form_hlrf(gfun,nvar);

    if nargout > 1
        if norm(gradStar)==0
            alpha2 = zeros(nvar,1);
        else
            alpha2 = (gradStar(:)./norm(gradStar)).^2;
        end
    end
end

function g = gfun_basic(u,model,rv,const,t,caseName,Cr_det)
    X = struct();
    for i = 1:numel(rv)
        if strcmpi(rv(i).name,'Cr')
            continue;
        end
        X.(rv(i).name) = u_to_x(rv(i),u(i));
    end

    X = add_constants(X,const,1);

    if strcmpi(caseName,'caseII')
        idxCr = find(strcmp({rv.name},'Cr'),1);
        Cr = u_to_x(rv(idxCr),u(idxCr));
    else
        Cr = Cr_det;
    end

    xc = eval_xc_model(model,X,t);
    g = Cr - xc(1);
end

function [beta,alpha2,names] = form_beta_An(rv,t,caseName,Cr_det,Cr_mean,Cr_cov)
    rv_form = rv;
    if strcmpi(caseName,'caseII')
        rv_form = add_rv(rv_form,'Cr','lognormal',Cr_mean,Cr_cov);
    end

    names = {rv_form.name}';
    nvar = numel(rv_form);
    gfun = @(u) gfun_An(u,rv_form,t,caseName,Cr_det);

    [beta,~,gradStar] = form_hlrf(gfun,nvar);

    if nargout > 1
        if norm(gradStar)==0
            alpha2 = zeros(nvar,1);
        else
            alpha2 = (gradStar(:)./norm(gradStar)).^2;
        end
    end
end

function g = gfun_An(u,rv,t,caseName,Cr_det)
    A = NaN; n = NaN; theta = NaN; Cr = Cr_det;

    for i = 1:numel(rv)
        name = rv(i).name;
        val = u_to_x(rv(i),u(i));
        switch name
            case 'A'
                A = val;
            case 'n'
                n = val;
            case 'theta_xc'
                theta = val;
            case 'Cr'
                Cr = val;
        end
    end

    xc = theta .* A .* t.^(0.5-n);
    g = Cr - xc;
end

function x = u_to_x(r,u)
    p = normcdf(u);
    p = min(max(p,1e-12),1-1e-12);

    switch lower(r.type)
        case 'normal'
            x = r.mean + r.std.*u;
        case 'lognormal'
            [mu_ln,sig_ln] = lognormal_mu_sigma(r.mean,r.cov);
            x = exp(mu_ln + sig_ln.*u);
        case 'bounded_power'
            x = r.a + (r.b-r.a).*p.^(1/r.gamma);
        otherwise
            error('Unknown RV type: %s',r.type);
    end
end

function [beta,u,grad] = form_hlrf(gfun,nvar)
    % Basic HL-RF FORM implementation in independent standard-normal space.
    % Failure domain is g(u) <= 0.
    maxIter = 80;
    tolU = 1e-4;
    tolG = 1e-5;

    u = zeros(nvar,1);
    g0 = gfun(u);

    for iter = 1:maxIter
        [g,grad] = finite_diff_grad(gfun,u);

        ng = norm(grad);
        if ng < 1e-12
            warning('FORM gradient close to zero. Returning beta = NaN.');
            beta = NaN;
            return;
        end

        alpha = grad(:)/ng;
        beta_lin = (g - grad(:)'*u)/ng;
        uNew = -beta_lin * alpha;

        % Damping improves robustness for nonlinear models
        lambda = 1.0;
        gNew = gfun(uNew);
        while abs(gNew) > abs(g) && lambda > 0.05
            lambda = lambda/2;
            uTry = u + lambda*(uNew-u);
            gNew = gfun(uTry);
            uNew = uTry;
        end

        if norm(uNew-u) < tolU && abs(gNew) < max(tolG,1e-4*abs(g0))
            u = uNew;
            break;
        end

        u = uNew;
    end

    [~,grad] = finite_diff_grad(gfun,u);
    beta = norm(u);

    % Preserve sign if mean point is already in failure domain
    if g0 < 0
        beta = -beta;
    end
end

function [g,grad] = finite_diff_grad(gfun,u)
    g = gfun(u);
    n = numel(u);
    grad = zeros(n,1);

    h = 1e-4;
    for i = 1:n
        du = zeros(n,1);
        du(i) = h;
        gp = gfun(u+du);
        gm = gfun(u-du);
        grad(i) = (gp-gm)/(2*h);
    end
end

function plot_reliability_figure(outDir,t,beta_target, ...
    fibI_mcs,fibI_form,fibII_mcs,fibII_form, ...
    malI_mcs,malI_form,malII_mcs,malII_form)

    fig = figure('Color','w','Units','centimeters','Position',[2 2 18 7]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile; hold on; box on; grid on;
    plot(t,fibI_mcs,'k-','LineWidth',1.8);
    plot(t,fibI_form,'k--','LineWidth',1.4);
    plot(t,fibII_mcs,'b-','LineWidth',1.8);
    plot(t,fibII_form,'b--','LineWidth',1.4);
    yline(beta_target,'k-.','HandleVisibility','off');
    text(45, beta_target+0.08, '$\beta_t=1.5$', ...
    'Interpreter','latex', ...
    'FontSize',9, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'BackgroundColor','w', ...
    'Margin',1);

    xlabel('Time, $t$ [years]','Interpreter','latex');
    ylabel('Reliability index, $\beta(t)$ [-]','Interpreter','latex');
    title('(a) fib Bulletin 34','FontWeight','bold');
    xlim([0 70]); ylim([0 4]);
    legend({'MCS, Case I','FORM, Case I','MCS, Case II','FORM, Case II'},'Location','northeast');

    nexttile; hold on; box on; grid on;
    plot(t,malI_mcs,'k-','LineWidth',1.8);
    plot(t,malI_form,'k--','LineWidth',1.4);
    plot(t,malII_mcs,'b-','LineWidth',1.8);
    plot(t,malII_form,'b--','LineWidth',1.4);
    yline(beta_target,'k-.','HandleVisibility','off');
    text(45, beta_target+0.08, '$\beta_t=1.5$', ...
    'Interpreter','latex', ...
    'FontSize',9, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'BackgroundColor','w', ...
    'Margin',1);

    xlabel('Time, $t$ [years]','Interpreter','latex');
    ylabel('Reliability index, $\beta(t)$ [-]','Interpreter','latex');
    title('(b) Malami','FontWeight','bold');
    xlim([0 70]); ylim([0 4]);
    legend({'MCS, Case I','FORM, Case I','MCS, Case II','FORM, Case II'},'Location','northeast');

    exportgraphics(fig,fullfile(outDir,'Figure9_reliability_beta_time.pdf'),'ContentType','vector');
    exportgraphics(fig,fullfile(outDir,'Figure9_reliability_beta_time.png'),'Resolution',300);
end

function plot_importance_basic(outDir,nF1,aF1,nF2,aF2,nM1,aM1,nM2,aM2)
    fig = figure('Color','w','Units','centimeters','Position',[2 2 18 14]);
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile; plot_barh_importance(nF1,aF1,'(a) fib Bulletin 34: Case I');
    nexttile; plot_barh_importance(nM1,aM1,'(b) Malami: Case I');
    nexttile; plot_barh_importance(nF2,aF2,'fib Bulletin 34: Case II');
    nexttile; plot_barh_importance(nM2,aM2,'Malami: Case II');

    xlabel(tl,'FORM importance factors, $\alpha_i^2$','Interpreter','latex','FontWeight','bold');

    exportgraphics(fig,fullfile(outDir,'Figure10_FORM_importance_basic.pdf'),'ContentType','vector');
    exportgraphics(fig,fullfile(outDir,'Figure10_FORM_importance_basic.png'),'Resolution',300);
end

function plot_importance_An(outDir,nF1,aF1,nF2,aF2,nM1,aM1,nM2,aM2)
    fig = figure('Color','w','Units','centimeters','Position',[2 2 18 8]);
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile; plot_barh_importance(nF1,aF1,'(a) fib Bulletin 34: Case I');
    nexttile; plot_barh_importance(nM1,aM1,'(b) Malami: Case I');
    nexttile; plot_barh_importance(nF2,aF2,'fib Bulletin 34: Case II');
    nexttile; plot_barh_importance(nM2,aM2,'Malami: Case II');

    xlabel(tl,'FORM importance factors, $\alpha_i^2$','Interpreter','latex','FontWeight','bold');

    exportgraphics(fig,fullfile(outDir,'Figure11_FORM_importance_An.pdf'),'ContentType','vector');
    exportgraphics(fig,fullfile(outDir,'Figure11_FORM_importance_An.png'),'Resolution',300);
end

function plot_barh_importance(names,alpha2,titleStr)
    alpha2 = alpha2(:);
    names = names(:);
    [alpha2,idx] = sort(alpha2,'ascend');
    names = names(idx);
    labels = latex_names_for_plot(names);

    barh(alpha2,'FaceColor',[0.55 0.80 0.60],'EdgeColor','k');
    grid on; box on;

    yticks(1:numel(labels));
    yticklabels(labels);
    set(gca,'TickLabelInterpreter','latex','FontName','Times New Roman','FontSize',9,'LineWidth',0.8);

    xmax = max(0.55,ceil((max(alpha2)+0.05)*10)/10);
    xlim([0 xmax]);
    title(titleStr,'Interpreter','latex','FontWeight','bold');

    for i = 1:numel(alpha2)
        text(alpha2(i)+0.01,i,sprintf('%.2f',alpha2(i)), ...
            'VerticalAlignment','middle','FontSize',8,'FontWeight','bold', ...
            'Interpreter','latex');
    end
end

function labels = latex_names_for_plot(names)
    % Exact mapping from internal variable names to publication-quality labels.
    % This avoids broken labels such as $3 x c$ or $R-1 ACC,0$ in exported PDFs.
    labels = cell(numel(names),1);

    for i = 1:numel(names)
        if isstring(names(i))
            nm = char(names(i));
        elseif iscell(names)
            nm = char(names{i});
        else
            nm = char(names(i));
        end

        switch nm
            case 'Racc_inv'
                labels{i} = '$R^{-1}_{\mathrm{ACC},0}$';
            case 'kt'
                labels{i} = '$k_t$';
            case 'eps_t'
                labels{i} = '$\varepsilon_t$';
            case 'tc'
                labels{i} = '$t_c$';
            case 'bc'
                labels{i} = '$b_c$';
            case 'CO2ppm'
                labels{i} = '$[CO_2]_{\mathrm{ppm}}$';
            case 'RH'
                labels{i} = '$RH$';
            case 'tw'
                labels{i} = '$t_w$';
            case 'bw'
                labels{i} = '$b_w$';
            case 'pSR'
                labels{i} = '$p_{SR}$';
            case 'theta_xc'
                labels{i} = '$\theta_{x_c}$';
            case 'W_C'
                labels{i} = '$W/C$';
            case 'T'
                labels{i} = '$T$';
            case 'Cr'
                labels{i} = '$C_r$';
            case 'A'
                labels{i} = '$A$';
            case 'n'
                labels{i} = '$n$';
            otherwise
                labels{i} = ['$', strrep(nm,'_','\_'), '$'];
        end
    end
end

function T = make_importance_table(nF1,aF1,nF2,aF2,nM1,aM1,nM2,aM2)
    T = table();
    T = [T; make_one_importance_table('fib_caseI',nF1,aF1)];
    T = [T; make_one_importance_table('fib_caseII',nF2,aF2)];
    T = [T; make_one_importance_table('malami_caseI',nM1,aM1)];
    T = [T; make_one_importance_table('malami_caseII',nM2,aM2)];
end

function T = make_one_importance_table(caseName,names,alpha2)
    n = numel(names);
    Case = repmat({caseName},n,1);
    Parameter = names(:);
    Alpha2 = alpha2(:);
    T = table(Case,Parameter,Alpha2);
end
