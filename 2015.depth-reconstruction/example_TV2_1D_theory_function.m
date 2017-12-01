% Minimize total variations of first derivative of the signal, given signal measurements
% Author: Luca Carlone
% date: 2015-11-20
function [results] = example_TV2_1D_theory_function(settings)
  
if nargin < 1
    clc; close all
    settings.epsilon = 0.2; 
    settings.testL2 = 1;
    settings.nrCorners = 2;  % number of critical points in the piecewise linear function
    settings.percSamples = 0.05; % percentage of measured pixels (samples)
    settings.maxValue = 5;
    settings.addNoise = 1;
    settings.N = 1000; % total number of points
    settings.optMode = 'l1inf'; % l1, l1inf, l1reg
    settings.lambda_reg = NaN;
    settings.noiseMode = 'linf' % l2, linf
    settings.sampleMode = 'uniform'; % inSegments, all, onlyCorners, alsoCorners, uniform
    settings.doAddNeighbors = 0; % 1 or 2
    settings.doAddBoundary = 1; 
    settings.isDebug = 0;
    settings.tol = 1e-7;
    %% path include
    addpath('./lib')
    addpath(genpath('./NESTA_Luca')) 
    addpath('../myLib')
    addpath('./testTheory')
end
format long
try
    cvx_solver mosek 
    cvx_save_prefs
end
  
%% create TV2 matrix
N = settings.N;
tol = settings.tol;
epsilon = settings.epsilon;
isDebug = settings.isDebug;
[TV2] = createFiniteDiff2(N);

%% Create piecewise linear function x (true signal)
[zGT, xq, yq] = loadGTsignal_1D(settings,TV2);

%% create sample set
[samples, K, percSamples] = createSampleSet_1D(settings, TV2, zGT, xq);
fprintf('Percentage of samples: %g\n',percSamples)

%% create sparse sampling matrix
Rfull = eye(N);
R = Rfull(samples, :);

%% create measurements
epsilonK = epsilon * sqrt(K);
if settings.addNoise == 1
    switch settings.noiseMode
        case 'l2'
            noise = epsilon * (2*rand(K,1)-1);
            if(norm(noise) > epsilonK) error('wrong l2 noise'); end
            error('deprecated noise model l2')
        case 'linf'
            noise = epsilon * (2*rand(K,1)-1); % entry-wise in [-eps,+eps]
            if(norm(noise, Inf) > epsilon) error('wrong linf noise'); end
        otherwise
            error('wrong choice of noiseMode')
    end
else
    noise = zeros(K,1);
end
y = R * zGT + noise;
noiseNorm = norm(noise); 
noiseInfNorm = norm(noise, Inf);

%% l1 minimization
R = sparse(R);
[zlinf,cvx_time_linf,tvMin] = minimizeL1cvx(TV2,N,R,y,settings.optMode,epsilon,epsilonK,settings.lambda_reg);
% DUAL_PROBLEM_l1inf % uncomment to test it

[lowBound,upBound,violated] = computeErrorBounds1D(zlinf,samples,y,epsilon,zGT,isDebug);
sampledInSegments = includeTwoSamplesPerSegment(TV2,zGT,samples);
if sampledInSegments && ~isempty(violated) && epsilon > 0
   error('error bounds are violated') 
end
if sampledInSegments==0
   disp('not sample twice in each segment') 
end

%% 2-stage algorithm
if settings.doAddNeighbors > 0
    [zlinf_s2,cvx_time_linf_s2,~,TVsigns] = stage2_1D(TV2,N,y,R,samples,zlinf,tvMin,epsilon,epsilonK,settings.optMode);
else
    zlinf_s2 = zlinf; cvx_time_linf_s2 = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% duality check with *squared* norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isDebug && strcmp(settings.optMode,'l2')) check_dual; end
if(isDebug && strcmp(settings.optMode,'l1inf') && epsilon > 0) 
    [dstar, cost_zreg, ve] = check_dual_l1infSquared(TV2,N,R,y,epsilon,zlinf); end
if ~exist('dstar') dstar = NaN; end
if ~exist('cost_zreg') cost_zreg = NaN; end

%% calculate naiveTight using active constraints
if epsilon > 0
    % COMPUTE_NAIVE_TIGHT % uncomment to try
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check sufficient conditions for exact recovery [Nam, Davies, Elad, Gribonval, 2012]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(isDebug) check_nam_etal; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check conditions in [Vaiter et al., 2013]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isDebug && epsilon > 0) check_vaiter_etal; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check conditions in [Fadili et al., 2013]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_fadili_etal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% duality check with unsquared norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_dual_unsquared

%% test l2 constrained optimization
if settings.testL2
    if epsilon == 0
        zl2 = zlinf; cvx_time_l2 = cvx_time_linf;
        zl2_s2 = zlinf_s2; cvx_time_l2_s2 = cvx_time_linf_s2;
    else
        [zl2,cvx_time_l2,tvMin_l2] = minimizeL1cvx(TV2,N,R,y,'l1',epsilon,epsilonK,settings.lambda_reg);
        if settings.doAddNeighbors > 0
            [zl2_s2,cvx_time_l2_s2,~,~] = stage2_1D(TV2,N,y,R,samples,zl2,tvMin_l2,epsilon,epsilonK,'l1');
        else
            zl2_s2 = zl2; cvx_time_l2_s2 = -1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Naive Approach: connect the dots
disp('================================================')
disp('Naive Approach.')
naive_time_temp = tic;
zNaive = interp1(samples(:),y,[1:N]','linear','extrap');
naive_time = toc(naive_time_temp);

disp('================================================')
disp('Nesta with Naive initialization.')
%% Fast solvers
if epsilon == 0
    % zl2_s2 = l1_reweighted_1D(TV2,N,y,R,samples,epsilon);
    %% very fast, best solution in terms of accuracy and time
    [zNesta,timeNesta,iterNesta] = l2_reweighted_1D(TV2,y,R,samples);
    % zl2_s2 = NESTA_1D(TV2,N,y,R,samples,epsilon);
else
    %% this converges to a solution which is better than l1, but slightly suboptimal
    %% actually things get worse when the noise is large, in which case convergence is poor
    % [zIrls,timeNesta,iterNesta] = l2_reweighted_1D_linf_alm(TV2,y,R,samples,epsilon,1e-4,1e-4)
    %% the following does not converge to the l1 solution in the allocated number of function evaluations
    % [zIrls,timeNesta,iterNesta] = l2_reweighted_1D_linf_fmincon(TV2,y,R,samples,epsilon);
    %% the following does not converge at all
    % [zIrls,timeNesta,iterNesta] = projSubgradient_1D_linf(TV2,y,R,samples,epsilon);
    % [zIrls,timeNesta,iterNesta] = projSubgradient_1D_linf(TV2,y,R,samples,epsilon);
    % [zIrls,timeNesta,iterNesta] = l2_reweighted_1D_linf_nesterov(TV2,y,R,samples,epsilon);
    %% CPLEX is not much faster
    % [zcplex, cplex_time] = minimizeL1cplex_1D(y,samples,N)
     
    %% NESTA is actually slower here: no worries, we don't really need speed up in the 1D case
    % the reasons is that to get a good solution, we should set mu small
    % and small tolerances which slows down convergence
    mu = 0.0002;
    opts = [];
    opts.TOlVar = 1e-7;
    opts.verbose = 0; % show intermediate results
    opts.maxintiter = 10;
    opts.errFcn = @(x) norm( TV2 * x , 1 ); % only for debug: displayed at the end of each iteration
    opts.D = TV2;
    opts.stoptest = 1;
    opts.xplug = zNaive;
    opts.typemin = 'TV2';
    opts.typecon = 'Linf';
    opts.samples = samples;
    opts.maxiter = 10000;
    counter();
    Rhandle = @(x) R*x;
    Rt = R';
    Rthandle = @(x) Rt*x;
    Ac = @(z) counter(Rhandle,z);
    Atc = @(z) counter(Rthandle,z);
    [zNesta,iterNesta,~,~,timeNesta] = myNESTA(Ac,Atc,y,mu,epsilon,opts);
end

%% plot: comparison among techniques
figure; clf
plot(zGT,'-g', 'LineWidth', 3); hold on;
plot(zlinf_s2,'-b', 'LineWidth', 4);
plot(zlinf,'--c', 'LineWidth', 3);
plot(zl2_s2,'-k', 'LineWidth', 4);
plot(zl2,'--m', 'LineWidth', 3);
plot(zNaive,'-r', 'LineWidth', 3);
plot(zNesta,':b', 'LineWidth', 3);
plot(samples,y,'o','markersize',10)
legend('zGT','zinf-2stage','zinf','zl2-2stage','zl2','znaive','zNesta')
title('comparison of all approaches')

%% create figure for paper
% f = figure; clf
% dim = 24
% plot(zGT,'-g', 'LineWidth', 5); hold on;
% plot(zlinf,'-.m', 'LineWidth', 4);
% plot(zlinf_s2,'-b', 'LineWidth', 5);
% plot(zNaive,':r', 'LineWidth', 3);
% % plot(samples,y,'or','markersize',10)
% leg = legend('$z^\diamond$','naive','L1 ($\varepsilon=0.5$)','A1 ($\varepsilon=0.5$)'); 
% set(leg,'Interpreter','latex'); 
% xlabel('')
% ylabel('depth profile')
% set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
% set(ylabh, 'FontSize', dim); xlabh=get(gca,'xlabel');
% set(xlabh, 'FontSize', dim);
% ax = gca;
% ax.XTick = [];
% ax.YTick = [];
% saveas(f,'exp4-test2','epsc'); 
% box off
% axis off

%% plot: visualize sign assigment for 2-stage algorithm
figure; clf; hold on
plot(samples,y,'o','markersize',10)
plot(zGT,'-g', 'LineWidth', 4); hold on;
if exist('TVsigns','var')
  indNeg = find(TVsigns == -1);
  indPos = find(TVsigns == +1);
  plot(indNeg,zGT(indNeg),'Xr','markersize',4)
  plot(indPos,zGT(indPos),'+k','markersize',4)
end
title('ground truth and samples')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% statistics and optimality
if isDebug && exist('x','var') && exist('xls','var')
    plot(x,'-r', 'LineWidth', 3);
    plot(xls,'--m', 'LineWidth', 3);
    disp('=============================================')
    disp('relative differences in solutions (only when epsilon > 0)')
    fprintf('zlinf - x / zlinf = %g\n',norm(zlinf - x) / norm(zlinf));
    fprintf('zlinf - xls / zlinf = %g\n',norm(zlinf - xls) / norm(zlinf));
end

disp('=============================================')
disp('objective value')
cost_zGT = norm(TV2 * zGT,1); fprintf('cost_zGT = %g\n',cost_zGT);
cost_zlinf = norm(TV2 * zlinf,1); fprintf('cost_zlinf = %g\n',cost_zlinf);
cost_zlinf_s2 = norm(TV2 * zlinf_s2,1); fprintf('cost_zlinf_s2 = %g\n',cost_zlinf_s2);
cost_zl2 = norm(TV2 * zl2,1); fprintf('cost_zl2 = %g\n',cost_zl2);
cost_zl2_s2 = norm(TV2 * zl2_s2,1); fprintf('cost_zl2_s2 = %g\n',cost_zl2_s2);
cost_naive = norm(TV2 * zNaive,1); fprintf('cost_naive = %g\n',cost_naive);
cost_nesta = norm(TV2 * zNesta,1); fprintf('cost_nesta = %g\n',cost_nesta);

disp('=============================================')
disp('reconstruction error')
[e_zlinf , e_zlinf_l1, e_zlinf_sm, e_zlinf_linf] = computeError(zlinf, zGT); 
fprintf('e_zlinf = %g, e_zlinf_l1 = %g, e_zlinf_sm = %g, e_zlinf_linf = %g\n',e_zlinf,e_zlinf_l1,e_zlinf_sm,e_zlinf_linf);
[e_zlinf_s2 , e_zlinf_s2_l1, e_zlinf_s2_l1_sm, e_zlinf_s2_l1_linf] = computeError(zlinf_s2, zGT); 
fprintf('e_zlinf_s2 = %g, e_zlinf_s2_l1 = %g, e_zlinf_s2_l1_sm = %g, e_zlinf_s2_l1_linf = %g\n',e_zlinf_s2,e_zlinf_s2_l1,e_zlinf_s2_l1_sm,e_zlinf_s2_l1_linf);
[e_zl2 , e_zl2_l1, e_zl2_sm, e_zl2_linf] = computeError(zl2, zGT); 
fprintf('e_zl2 = %g, e_zl2_l1 = %g, e_zl2_sm = %g, e_zl2_linf = %g\n',e_zl2,e_zl2_l1,e_zl2_sm,e_zl2_linf);
[e_zl2_s2 , e_zl2_s2_l1, e_zl2_s2_sm, e_zl2_s2_linf] = computeError(zl2_s2, zGT); 
fprintf('e_zl2_s2 = %g, e_zl2_s2_l1 = %g, e_zl2_s2_sm = %g, e_zl2_s2_linf = %g\n',e_zl2_s2,e_zl2_s2_l1,e_zl2_s2_sm,e_zl2_s2_linf);
[e_naive , e_naive_l1, e_naive_sm, e_naive_inf] = computeError(zNaive, zGT); 
fprintf('e_naive = %g, e_naive_l1 = %g, e_naive_sm = %g, e_naive_inf = %g\n',e_naive,e_naive_l1,e_naive_sm,e_naive_inf);
[e_nesta , e_nesta_l1, e_nesta_sm, e_nesta_inf] = computeError(zNesta, zGT); 
fprintf('e_nesta = %g, e_nesta_l1 = %g, e_nesta_sm = %g, e_nesta_inf = %g\n',e_nesta , e_nesta_l1, e_nesta_sm, e_nesta_inf);

zlinf_l0 = getl0norm(TV2 * zlinf); 
zlinf_s2_l0 = getl0norm(TV2 * zlinf_s2);
zl2_l0 = getl0norm(TV2 * zl2);
zl2_s2_l0 = getl0norm(TV2 * zl2_s2);
zNaive_l0 = getl0norm(TV2 * zNaive);
zNesta_l0 = getl0norm(TV2 * zNesta);
zGT_l0 = getl0norm(TV2 * zGT); 

%% does the 2 stage approach work? 
if settings.addNoise == 0 && epsilon == 0 && e_zlinf_s2 > 1e-2
   e_zlinf_s2, warning('nonexact reconstruction') 
end

%% test optimality
if (isDebug)
    if  epsilon == 0 
        tol = 1e-6;
        nonSamples = setdiff([1:N],samples);
        [isOpt_z, ~, ujinf]= testOpt(TV2,zlinf,nonSamples,tol)
        if(isOpt_z > tol) error('optimality test does not work - zlinf'); end
        fprintf('optimality test - zlinf: norm(uj,Inf) = %g\n',ujinf)
        
        [isOpt_zGT, ~, ujinf] = testOpt(TV2,zGT,nonSamples,tol);
        if(isOpt_zGT > tol) warning('optimality test does not work - zGT'); end
        fprintf('optimality test - zGT: norm(uj,Inf) = %g\n',ujinf)
        
        [isOpt_rand, ~, ujinf] = testOpt(TV2,rand(N,1),nonSamples,tol);
        if(isOpt_rand < 1e-3) fprintf('isOpt = %g\n',isOpt_rand); error('optimality test is vacuous - rand'); end
        fprintf('optimality test - zrand: norm(uj,Inf) = %g\n',ujinf)
    else
        %% test optimality linf
        isOpt_z = testOptLinf(TV2,zlinf,samples,y,epsilon,tol)
        isOpt_zGT = testOptLinf(TV2,zGT,samples,y,epsilon,tol)
        isOpt_rand = testOptLinf(TV2,rand(N,1),samples,y,epsilon,tol)
    end
end

disp('=============================================')
disp('other statistics')
fprintf('noiseNorm = %g\n',noiseNorm);
fprintf('noiseInfNorm = %g\n',noiseInfNorm);
% The root-mean-square deviation (RMSD) or root-mean-square error (RMSE)
% RMSD = || y - y_hat || / sqrt(n) 
RMSD = epsilon / sqrt(K); fprintf('RMSD = %g\n',RMSD);
% signal to noise ratio: alternative definition of SNR (https://en.wikipedia.org/wiki/Signal-to-noise_ratio)
snr = ( max(abs(zGT)) / max(abs(noise)) )^2; fprintf('snr = %g\n',snr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zGT_norm2 = norm(zGT); zGT_norm1 = norm(zGT,1); zGT_normInf = norm(zGT,Inf);
results = struct('settings',settings,...
    'e_zlinf',e_zlinf,'e_zlinf_s2',e_zlinf_s2,'e_zl2',e_zl2,'e_zl2_s2',e_zl2_s2,'e_naive',e_naive,'e_nesta',e_nesta,...
    'e_zlinf_l1',e_zlinf_l1,'e_zlinf_s2_l1',e_zlinf_s2_l1,'e_zl2_l1',e_zl2_l1,'e_zl2_s2_l1',e_zl2_s2_l1,'e_naive_l1',e_naive_l1,'e_nesta_l1',e_nesta_l1,...
    'e_zlinf_sm',e_zlinf_sm,'e_zlinf_s2_l1_sm',e_zlinf_s2_l1_sm,'e_zl2_sm',e_zl2_sm,'e_zl2_s2_sm',e_zl2_s2_sm, 'e_naive_sm',e_naive_sm,'e_nesta_sm',e_nesta_sm,...
    'e_zlinf_linf',e_zlinf_linf,'e_zlinf_s2_l1_linf',e_zlinf_s2_l1_linf,'e_zl2_linf',e_zl2_linf,'e_zl2_s2_linf',e_zl2_s2_linf,'e_naive_inf',e_naive_inf,'e_nesta_inf',e_nesta_inf,...
    'dstar',dstar,'cost_zreg',cost_zreg,...
    'cost_zlinf',cost_zlinf,'cost_zlinf_s2',cost_zlinf_s2,'cost_zl2',cost_zl2,'cost_zl2_s2',cost_zl2_s2,'cost_naive',cost_naive,'cost_zGT',cost_zGT,'cost_nesta',cost_nesta, ...
    'cvx_time_linf',cvx_time_linf,'cvx_time_linf_s2',cvx_time_linf_s2,'cvx_time_l2',cvx_time_l2,'cvx_time_l2_s2',cvx_time_l2_s2,'naive_time',naive_time,...
    'timeNesta',timeNesta,'iterNesta',iterNesta,...,
    'zlinf_l0', zlinf_l0, 'zlinf_s2_l0', zlinf_s2_l0, 'zl2_l0', zl2_l0, 'zl2_s2_l0', zl2_s2_l0, 'zNaive_l0', zNaive_l0, 'zNesta_l0', zNesta_l0, 'zGT_l0', zGT_l0, ...
    'zGT_norm2',zGT_norm2,'zGT_norm1',zGT_norm1,'zGT_normInf',zGT_normInf,...
    'sampledInSegments',sampledInSegments,'K',K,'N',N,'percSamples',percSamples,'epsilonK',epsilonK,...
    'noiseNorm',noiseNorm,'noiseInfNorm',noiseInfNorm,'RMSD',RMSD,'snr',snr);
disp('- saved results')
disp('=============================================')


