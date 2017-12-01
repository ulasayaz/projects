% testing NESTA code

clc; close all
settings.epsilon = 0.2;
settings.testL2 = 1;
settings.nrCorners = 2;  % number of critical points in the piecewise linear function
settings.percSamples = 0.05; % percentage of measured pixels (samples)
settings.maxValue = 5;
settings.addNoise = 1;
settings.N = 2500; % total number of points
settings.optMode = 'l1'; % l1, l1inf, l1reg
settings.lambda_reg = 1;
settings.noiseMode = 'linf' % l2, linf
settings.sampleMode = 'uniform'; % inSegments, all, onlyCorners, alsoCorners, uniform
settings.doAddNeighbors = 0; % 1 or 2
settings.doAddBoundary = 1;
settings.isDebug = 0;
settings.tol = 1e-7;

%% path include
% addpath('../l1magic/Optimization')
% addpath('../grammSchmidt')
% addpath('../libOPL')
addpath('./algorithmsUA')
addpath(genpath('./NESTA')) 
addpath('../myLib')
addpath('./testTheory')
addpath('./lib')

format long
try
    cvx_solver mosek
    cvx_save_prefs
end
tryCPLEX = 0;

%% create TV2 matrix
N = settings.N;
tol = settings.tol;
epsilon = settings.epsilon;
isDebug = settings.isDebug;
%[TV2] = createFiniteDiff2(N);
[TV2] = createFiniteDiff1(N);

%% Create piecewise linear function x (true signal)
correct = false;
while correct == false
    % define a few query points (critical points)
    supp = randsample(2:N-1, settings.nrCorners);
    xq = [1, supp, N];
    yq = (rand(1, settings.nrCorners+2) - 0.5) * settings.maxValue;
    % generate linear interpolation
    zGT = interp1(xq,yq,1:N,'linear')';
    % let minimum be above zero
    minZGT = min(zGT);
    zGT = zGT - minZGT + 0.1;
    yq = yq - minZGT + 0.1;
    % if it did not end up creating a quasi straight line accept signal
    if length(find(abs(TV2*zGT)>1e-3)) == settings.nrCorners 
        correct = true;
    end
end

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
%[zlinf,cvx_time_linf,tvMin] = minimizeL1cvx(TV2,N,R,y,settings.optMode,epsilon,epsilonK,settings.lambda_reg);

%%%%%%%%%%%%% IMPORTANT TO REMEMBER ABOUT THIS PART
cvx_begin
cvx_quiet false
cvx_precision best
variable xTest(N,1);
minimize( norm(TV2 * xTest,1) )
% minimize( norm(xTest,1) )
subject to
norm(R * xTest - y) <= epsilon;
cvx_end

zlinf = xTest;
%%%%%%%%%%%%%%%%

if epsilon > 0
    %% plot: ground truth & samples & recovered

    fig1 = figure; clf
    plot(zGT,'-g', 'LineWidth', 3); hold on;
    plot(zlinf,'--k', 'LineWidth', 3);
    %plot(zlinf_s2,'-b', 'LineWidth', 3);
    %plot(zNaiveTight,'-b', 'LineWidth', 3);
    plot(samples,y,'o','markersize',3)
    errorbar(samples,y,epsilon*ones(size(y)),'.m')
    %     errorbar(tightUpSamples,y(tightUpSamplesIds),epsilon*ones(size(y(tightUpSamplesIds))),'.b','LineWidth', 2)
    %     errorbar(tightDownSamples,y(tightDownSamplesIds),epsilon*ones(size(y(tightDownSamplesIds))),'.b','LineWidth', 2)
    % plot(samples(tightSamples),zlinf(samples(tightSamples)), 'dr', 'MarkerFaceColor', 'r', 'markersize',8)
    legend('zGT','zlinf','zNaiveTight')
    grid minor
    
end

if epsilon == 0
    % zl2_s2 = l1_reweighted_1D(TV2,N,y,R,samples,epsilon);
    %% very fast, best solution in terms of accuracy and time
    [zIrls,timeIrls,iterIrls] = l2_reweighted_1D(TV2,y,R,samples);
    % zl2_s2 = NESTA_1D(TV2,N,y,R,samples,epsilon);
else
    %% this converges to a solution which is better than l1, but slightly suboptimal
    %% actually things get worse when the noise is large, in which case convergence is poor
    % [zIrls,timeIrls,iterIrls] = l2_reweighted_1D_linf_alm(TV2,y,R,samples,epsilon,1e-4,1e-4)
    %% the following does not converge to the l1 solution in the allocated number of function evaluations
    % [zIrls,timeIrls,iterIrls] = l2_reweighted_1D_linf_fmincon(TV2,y,R,samples,epsilon);
    %% the following does not converge at all
    % [zIrls,timeIrls,iterIrls] = projSubgradient_1D_linf(TV2,y,R,samples,epsilon);
    % [zIrls,timeIrls,iterIrls] = projSubgradient_1D_linf(TV2,y,R,samples,epsilon);
    %[zIrls,timeIrls,iterIrls] = l2_reweighted_1D_linf_nesterov(TV2,y,R,samples,epsilon);
    zl2_s2 = NESTA_1D(TV2,N,y,R,samples,epsilon);
    
    %% test the original nesta online
    fprintf('2) run NESTA\n\n');
    U = @(z) z;
    Ut = @(z) z;
    mu = 0.2; %--- can be chosen to be small
    opts = [];
    opts.maxintiter = 5;
    opts.TOlVar = 1e-5;
    opts.verbose = 0;
    opts.maxiter = 5000;
    opts.U = TV2;
    opts.Ut = TV2';
    opts.stoptest = 1;
    %opts.typemin = 'l1';
    opts.typemin = 'tv';
    % counter();
    % Ac = @(z) counter(A,z);
    % Atc = @(z) counter(At,z);

    [x_nesta,niter,resid,err] = NESTA(R,R',y,mu,epsilon,opts);
    
    disp('optimal cost with cvx')
    norm( TV2 * zlinf ,1 )
    disp('iterated cost with nesta online')
    norm( TV2 * zl2_s2 , 1)
    %norm( TV2 * zIrls , 1)
    disp('iterated cost with nesta online 2')
    norm( TV2 * x_nesta , 1)
    
    plot(zl2_s2,'color','r');
    
end


























