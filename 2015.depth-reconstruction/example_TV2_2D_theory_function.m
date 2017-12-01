% Minimize total variations of first derivative of the signal, given signal measurements
% Author: Luca Carlone
% date: 2015-11-20
function [results] = example_TV2_2D_theory_function(settings)

if nargin < 1
    clc; close all    
    settings.mode = 'pwLinear'; % pwLinear, toy, saddle, stripe
    settings.pwOption = 'diagonal'; % pyramid, diagonal
    settings.d = 100; % must be odd for testing
    settings.sampleMode = 'uniform';
    % edgesRGBrandom, edgesOnly, uniform, imgCornersOnly, cornersOnly, h_per_vline, h_per_vhline, h_per_area, grid, vlines, hvlines, hvgrid, randlines
    settings.h = 10; % depending on the previous setting, this is: h_per_area, size of grid cell, nr vlines, nothing
    settings.epsilon = 0.1;
    settings.addNoise = 1;
    settings.optMode = 'l1inf'; % l1, l1inf, l1reg
    settings.nrCorners = 3;
    settings.percSamples = 0.05; % percentage of measured pixels (samples)
    settings.maxValue = 5;
    settings.lambda_reg = NaN;
    settings.noiseMode = 'linf' % l2, linf
    settings.doAddNeighbors = 0; 
    settings.nrNeighbors = 1; 
    settings.includeAllDirs = 1; % 0 or 1
    settings.doAddBoundary = 0;
    settings.doAddRandom = 0;
    settings.isDebug = 0;
    settings.isBounded = 0;
    settings.testDiag = 1;
    settings.tol = 1e-5;
    settings.mu = 0.001;
    %% Only for real images:
    % settings.imageFilename = '/Users/fangchangma/Dropbox (MIT)/Research/sparse sensing/sparseSensing/datasets/gazebo_simulation/2016.01.12_house/depth/depth002.png';
    settings.dataset = 'gazebo';    % gazebo, zed, ortho_proj_gazebo
    % LUCA
    settings.path = '~/Dropbox (MIT)/sparseSensing/datasets/';
    % FANGCHANG
    % settings.path = '/Users/fangchangma/Dropbox (MIT)/Research/sparse sensing/sparseSensing/datasets/';
    % ULAS
    % settings.path = '/Users/ulasayaz/Desktop/Dropbox (MIT)/sparseSensing/';
    settings.thresVect = []; % using default threshold for canny algorithm
    settings.imgID = 7; % gazebo: 3, zed: 400
    settings.subSample = 0.1;
    settings.useInputImage = 0;
    settings.doSaveResults = 0;
    settings.resultsFolder = './resultsIROS/';
    %% path include
    addpath('../myLib')
    addpath('./lib')
    addpath('./testTheory')
    addpath('../libOPL')
    addpath('./testTheory')
    addpath(genpath('./NESTA_Luca')) 
end
if isfield(settings,'isBounded')==0
    settings.isBounded = 0;
    warning('isBounded undefined in original settings')
end
if isfield(settings,'mu')==0
    settings.mu = 0.001;
    warning('mu undefined in original settings')
end
format long
try
    cvx_solver mosek 
    cvx_save_prefs
end

%% create true signal
[ZGT, zGT, Nx, Ny, N, rgb] = loadGTsignal_2D(settings);
settings.Nx = Nx; settings.Ny = Ny; settings.N = N; tol = settings.tol; 
d = settings.d; epsilon = settings.epsilon; isDebug = settings.isDebug; mu = settings.mu;

%% Create TV matrices
[H,V] = createFiniteDiff2(Nx,Ny);
settings.H = H; settings.V = V;
[H1,V1] = createFiniteDiff1(Nx,Ny);
[VZ_GT, ZH_GT, normGrad_GT, VZ_uncut_GT, ZH_uncut_GT] = computeGrad(ZGT,V,H);
[areas,labels] = segmentBinaryImageAndGetAreas(normGrad_GT);

figure;
set(gcf,'name','original laplacians')
subplot(131); imshow(zeroOneTh(abs(VZ_uncut_GT),tol));   title('V * ZGT');
subplot(132); imshow(zeroOneTh(abs(ZH_uncut_GT),tol));   title('ZGT * H')
subplot(133); imshow(zeroOneTh(normGrad_GT,tol));  title('ZGT * grad')

%% create sample set
[samples, K, percSamples] = createSampleSet_2D(settings, ZGT, labels, rgb);
fprintf('Percentage of samples: %g\n',percSamples)

%% create sparse sampling matrix
Rfull = speye(N);
R = Rfull(samples, :);

%% create measurements
epsilonK = epsilon * sqrt(K);
if settings.addNoise == 1
    switch settings.noiseMode
        case 'l2'
            noise = epsilon * (2*rand(K,1)-1);
            if(norm(noise) > epsilonK) error('wrong l2 noise'); end
            error('deprecated choise of noiseMode: l2')
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

%% plot samples
% [X_sample, Y_sample] = ind2sub([Nx Ny],samples);
[Xq, Yq] = meshgrid(1:size(ZGT,2), 1:size(ZGT,1));
X_sample = Xq(samples)';
Y_sample = Yq(samples)';
if(isDebug)
    figure(101); clf
    if settings.useInputImage
        surf(ZGT); hold on; shading faceted % flat, axis equal;
    else
        surf(1:Nx, 1:Ny, ZGT); hold on; shading faceted % flat, axis equal;
    end 
    plot3(X_sample, Y_sample, y, 'or', 'markersize', 10);
    set(gcf,'name','original with samples');
end

%% create TV2 matrix for optimization
% vec(V * Z) = (I_Ny kron V) * vec(Z)
TV2_V = kron(speye(Ny), V);
if norm(TV2_V * zGT - vec(V * ZGT) ) > 1e-4
    error('vectorization 1 was wrong')
end
% vec(Z * H) = (H' kron I_Nx) * vec(Z)
TV2_H = kron(H',speye(Nx));
if norm(TV2_H * zGT - vec(ZGT * H) ) > 1e-4
    error('vectorization 2 was wrong')
end
% assemble matrix to be used for optimization
TV2 = sparse([TV2_V; TV2_H]); % not using diagonal terms

if size(TV2,1) ~=  2 * (N -Nx-Ny)
    error('did not understand dimension of TV2')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check conjecture and null space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~settings.useInputImage && isDebug) check_conjecture; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check sufficient conditions for exact recovery [Nam, Davies, Elad, Gribonval, 2012]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~settings.useInputImage && isDebug) check_nam_etal; end % takes as inputs zGT, and TV2

%% l1 minimization
disp('solving l1 minimization with linf constraints')
R = sparse(R);
[z,cvx_time] = minimizeL1cvx(TV2,N,R,y,settings.optMode,epsilon,epsilonK,settings.lambda_reg,settings.isBounded,settings.maxValue);
Z = reshape(full(z),Nx,Ny);
if abs( norm(TV2 * z,1) - ( norm( vec(Z * H), 1) + norm( vec(V * Z), 1) ) ) >1e-5
    error('mismatch between matrix and vectorized objective')
end

%% Naive Approach: connect the dots
disp('================================================')
disp('Naive Approach.')
naive_time_temp = tic;
Fun = scatteredInterpolant(X_sample, Y_sample, y, 'linear');
Znaive = Fun(Xq, Yq);
naive_time = toc(naive_time_temp);
if settings.isBounded
    Znaive(Znaive > settings.maxValue) = settings.maxValue; % saturate large values to maxValue
    Znaive(Znaive < 0) = 0; % saturate negative values to zero
end
zNaive = vec(Znaive);

%% FAST SOLVER (hopefully)
disp('================================================')
disp('Nesta with naive initialization')
%% very fast, best solution in terms of accuracy and time
% reweighted algorithm works great with eps = 0 but it is not explained in
% the paper, s we stick to NESTA for exact recovery too
%if epsilon == 0
%    [zFast,timeFast,iterFast] = l2_reweighted_1D(TV2,y,R,samples,1e-2);
%else
    % this is not much better than cvx
    % [zFast,timeFast,iterFast] = l2_reweighted_1D_linf(TV2,y,R,samples,epsilon,1e-2,1e-2);
    % %% the following do not improve a lot
    % [zFast,timeFast,iterFast] = fmincon_1D_linf(TV2,y,R,samples,epsilon);
    % [zFast,timeFast,iterFast] = projSubgradient_1D_linf(TV2,y,R,samples,epsilon);
    %[zFast,timeFast,iterFast] = l2_reweighted_1D_linf_lmdivide(TV2,y,R,samples,epsilon,1e-2,1e-2);
    %[zFast,timeFast,iterFast] = l2_reweighted_1D_linf_inverse(TV2,y,R,samples,epsilon,1e-2,1e-2);
    % [zFast,timeFast,iterFast] = l2_reweighted_1D_linf_quadProg(TV2,y,R,samples,epsilon,1e-3);
    %% Irresponsive:
    % [zFast,timeFast,iterFast] =  l2_reweighted_1D_linf_linProg(TV2,y,R,samples,epsilon);
    % [zFast,timeFast,iterFast] =  l2_reweighted_1D_linf_nesterov(TV2,y,R,samples,epsilon);
    % error('still looking for a fast solver here')
    %% CPLEX is not faster
    % [zcplex, cplex_time] = minimizeL1cplex(TV2,y,samples,N)
    [zFast,timeFast,iterFast] = solve_nesta_2D(TV2,y,R,samples,epsilon,'TV2',zNaive,mu);
%end
Zfast = reshape(full(zFast),Nx,Ny);
disp('================================================')

figure; set(gcf,'name','evaluation IRLS');
subplot(1,3,1)
surf(ZGT); hold on; shading faceted % flat, axis equal;
plot3(X_sample, Y_sample, y, 'or', 'markersize', 10); title('GT')
subplot(1,3,2)
surf(Zfast); title('Zfast')
subplot(1,3,3)
surf(Z); title('Z l1')

%% compute bounds
if strcmp(settings.sampleMode,'grid')
 [lowBound,upBound,violated] = computeErrorBounds2D(Z,samples,y,epsilon,ZGT);
end
warning('complete bounds')
% sampledInSegments = includeTwoSamplesPerSegment(TV2,zGT,samples);
% if sampledInSegments && ~isempty(violated) && epsilon > 0
%    error('error bounds are violated') 
% end
% if sampledInSegments==0
%    warning('not sample twice in each segment') 
% end

%% l1 minimization with diagonal terms in the objective
if settings.testDiag
   disp('****************** using diagonal terms! *******************');
   [Hd,Vd] = createFiniteDiffDiag(Nx,Ny);
   % vec(ABC) = (C' \kron A) vec(B) % Vd * X * Hd
   TV2_xy = kron(Hd',Vd);  
   % vec(ABC) = (C' \kron A) vec(B)
   % lambda * norm( vec(V1 * Z * H1), 1) 
   % TV2_xy = 2 * kron(H1',V1); 
   TV2_D = sparse([TV2_V; TV2_H; TV2_xy]);
   
   %% test with CVX
   [z_diag,cvx_time_diag] = minimizeL1cvx(TV2_D,N,R,y,settings.optMode,epsilon,epsilonK,settings.lambda_reg,settings.isBounded,settings.maxValue);
   Z_diag = reshape(z_diag,Nx,Ny);
   [e_z_diag,e_z_diag_l1,e_z_diag_sm,e_z_diag_inf] = computeError(Z_diag, ZGT); 
   cost_z_diag = norm(TV2_D * z_diag,1);
   z_diag_l0 = getl0norm(TV2_D * z_diag);
   
   %% test with NESTA
   [zFast_diag,timeFast_diag,iterFast_diag] = solve_nesta_2D(TV2_D,y,R,samples,epsilon,'TV2_diag',zNaive,mu);
   Zfast_diag = reshape(zFast_diag,Nx,Ny);
   [e_zFast_diag,e_zFast_diag_l1,e_zFast_diag_sm,e_zFast_diag_inf] = computeError(Zfast_diag, ZGT); 
   cost_zFast_diag = norm(TV2_D * zFast_diag,1);
   zFast_diag_l0 = getl0norm(TV2_D * zFast_diag); 
   
   % no longer important
   cost_zGT_diag = norm(TV2_D * zGT,1);
   zNaive_diag_l0 = getl0norm(TV2_D * zNaive);
   zGT_diag_l0 = getl0norm(TV2_D * zGT);
else
    cvx_time_diag = NaN; cost_z_diag = NaN; e_z_diag = NaN;
    e_z_diag_l1 = NaN; e_z_diag_sm = NaN; e_z_diag_inf = NaN; cost_zGT_diag = NaN;
    z_diag_l0 = NaN; zNaive_diag_l0 = NaN; zGT_diag_l0 = NaN;
    e_zFast_diag = NaN; e_zFast_diag_l1 = NaN; e_zFast_diag_sm = NaN; e_zFast_diag_inf= NaN;
    timeFast_diag = NaN; iterFast_diag = NaN;
    cost_zFast_diag = NaN; zFast_diag_l0 = NaN;
end

disp('=============================================')
disp('objective value')
cost_zGT = norm(TV2 * zGT,1); fprintf('cost_zGT = %g\n',cost_zGT);
cost_z = norm(TV2 * z,1); fprintf('cost_z = %g\n',cost_z);
cost_zNaive = norm(TV2 * zNaive,1); fprintf('cost_zNaive = %g\n',cost_zNaive);
cost_zFast = norm(TV2 * zFast,1); fprintf('cost_zFast = %g\n',cost_zFast);
if(settings.testDiag) 
    fprintf('cost_z_diag = %g\n',cost_z_diag); 
    fprintf('cost_zFast_diag = %g\n',cost_zFast_diag); 
end
cost_diff = cost_zGT - cost_z; fprintf('cost_diff = %g\n',cost_diff);

disp('=============================================')
disp('reconstruction error')
[e_z , e_z_l1, e_z_sm, e_z_inf] = computeError(Z, ZGT); 
fprintf('e_z = %g, e_z_l1 = %g, e_z_sm = %g, e_z_inf = %g\n',e_z,e_z_l1,e_z_sm,e_z_inf);
[e_naive , e_naive_l1, e_naive_sm, e_naive_inf] = computeError(Znaive, ZGT); 
fprintf('e_naive = %g, e_naive_l1 = %g, e_naive_sm = %g, e_naive_inf = %g\n',e_naive,e_naive_l1,e_naive_sm,e_naive_inf);
[e_zFast , e_zFast_l1, e_zFast_sm, e_zFast_inf] = computeError(Zfast, ZGT); 
fprintf('e_zFast = %g, e_zFast_l1 = %g, e_zFast_sm = %g, e_zFast_inf = %g\n',e_zFast,e_zFast_l1,e_zFast_sm,e_zFast_inf);
if(settings.testDiag) 
    fprintf('e_z_diag = %g, e_z_diag_l1 = %g, e_z_diag_sm = %g, e_z_diag_inf = %g\n',e_z_diag,e_z_diag_l1,e_z_diag_sm,e_z_diag_inf); 
    fprintf('e_zFast_diag = %g, e_zFast_diag_l1 = %g, e_zFast_diag_sm = %g, e_zFast_diag_inf = %g\n',e_zFast_diag,e_zFast_diag_l1,e_zFast_diag_sm,e_zFast_diag_inf); 
end

disp('=============================================')
disp('timing')
fprintf('cvx_time = %g\n',cvx_time);
fprintf('naive_time = %g\n',naive_time);
fprintf('cvx_time_diag = %g\n',cvx_time_diag);
fprintf('timeFast = %g\n',timeFast);
fprintf('iterFast = %g\n',iterFast);
if(settings.testDiag)
    fprintf('timeFast_diag = %g\n',timeFast_diag);
    fprintf('iterFast_diag = %g\n',iterFast_diag);
end

z_l0 = getl0norm(TV2 * z); 
zNaive_l0 = getl0norm(TV2 * zNaive);
zGT_l0 = getl0norm(TV2 * zGT); 
zFast_l0 = getl0norm(TV2 * zFast); 

if(isDebug && epsilon > 0)
    cvx_begin
    cvx_precision('best')
    variable Z(Nx,Ny);
    minimize( norm( vec(Z * H), 1) + norm( vec(V * Z), 1) )
    subject to
    norm(y - R * vec(Z)) <= epsilonK;
    cvx_end
    [errorTVl2_test] = computeError(Z,ZGT);
    if norm(errorTVl2_test - e_z) > 1e-3
        error('cost mismatch in debug mode')
    end
end

%% COMPUTE GRADIENTS
[VZ, ZH, normGrad] = computeGrad(Z,V,H);
[VZ_naive, ZH_naive, normGrad_naive] = computeGrad(Znaive,V,H);
[~, ~, normGradFast] = computeGrad(Zfast,V,H);

%% PLOTS
ss = get(0,'screensize'); %The screen size
width = ss(3);
height = ss(4);
vert = 400; %vertical pixels
horz = 400; %horizontal pixels
pos = 0; % arranges the figure positions

%% compute mask and visualize
index = zeros(N, 1);
index(samples) = 1;
mask = reshape(index, size(ZGT));
Fig = figure(10); set(Fig,'Position',[pos * horz, (height - pos * vert), horz, vert]); pos = pos+1;
hold on
%% column 1: gt and samples
subplot(3,4,1); imshow(zeroOne(ZGT)); title('img');
subplot(3,4,5); imshow(mask); hold on; title(['samples: ',num2str(percSamples)]);
plot(X_sample,Y_sample,'xb','markersize',5)
%% column 2: naive
subplot(3,4,2); imshow(zeroOne(Znaive)); title('naive')
subplot(3,4,6); imshow(zeroOne(Znaive-ZGT)); 
title(sprintf('naive error (l2:%2.2g)',e_naive))
subplot(3,4,10); imshow(zeroOneTh(normGrad_naive)); title('grad naive')
%% column 3: l1
subplot(3,4,3); imshow(zeroOne(Z)); title('l1')
subplot(3,4,7); imshow(zeroOne(Z-ZGT)); 
title(sprintf('l1 error (l2:%2.2g)',e_z)) 
subplot(3,4,11); imshow(zeroOneTh(normGrad)); title('grad l1')
%% column 3: lfast
subplot(3,4,4); imshow(zeroOne(Zfast)); title('Zfast')
subplot(3,4,8); imshow(zeroOne(Zfast-ZGT)); 
title(sprintf('fast error (l2:%2.2g)',e_zFast)) 
subplot(3,4,12); imshow(zeroOneTh(normGradFast)); title('grad fast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check signs z0 inside each patch:
if(isDebug)
    CHECK_SIGNS_INPATCH
end

%% test optimality
if (isDebug)
    tol = 1e-4;
    isOptNorm_z = testOpt(TV2, z, nonSamples, tol)
    if(isOptNorm_z > tol) error('optimality test failed - z'); end
    
    [isOptNorm_zGT,sigm_GT] = testOpt(TV2, zGT, nonSamples, tol);
    isOptNorm_zGT
    if(isOptNorm_zGT > tol) warning('optimality test failed - zGT'); end
    
    % isOptNorm_rand = testOpt(TV2,rand(N,1),nonSamples,tol);
    % if(isOptNorm_rand < 0.1) fprintf('isOpt = %g\n',isOptNorm_rand); error('optimality test is vacuous - rand'); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(settings.doSaveResults)
    if length(settings.thresVect) > 0 % we manually specify thresholds
        folderName = createFilename(settings,horzcat(settings.resultsFolder,'results_canny_',settings.dataset));
    else     
        folderName = createFilename(settings,horzcat(settings.resultsFolder,'results_',settings.dataset));
    end
    mkdir(folderName);
    imFolderName = horzcat(folderName,'/');
    mkdir(imFolderName);
    filename = horzcat(imFolderName,num2str(settings.imgID));
    
    if settings.maxValue ~= 10
       error('maxValue should be 10 with gazebo and zed') 
    end
    
    Fig = figure; 
    hold on
    subplot(231); imshow(ZGT/settings.maxValue); title(horzcat('full signal'));
    subplot(232); imshow(mask); title(sprintf('samples:%2.2g',percSamples));
    subplot(234); imshow(Znaive/settings.maxValue); title(sprintf('naive (error:%2.2g)',e_naive))
    subplot(235); imshow(Z/settings.maxValue); title(sprintf('l1 rec (error:%2.2g)',e_z))
    subplot(236); imshow(Z_diag/settings.maxValue); title(sprintf('l1 diag (error:%2.2g)',e_z_diag))
    saveas(Fig,horzcat(filename,'-overview'),'png'); 
    
    Fig = figure();
    imagesc(ZGT/settings.maxValue, [0 1]); 
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    saveas(Fig,horzcat(filename,'-original'),'png');
 
    Fig = figure();
    imagesc(mask, [0 1]);
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    saveas(Fig,horzcat(filename,'-samples'),'png');

    Fig = figure();
    imagesc(Z/settings.maxValue, [0 1]);
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    saveas(Fig,horzcat(filename,'-l1rec'),'png');

    Fig = figure();
    imagesc(Z_diag/settings.maxValue, [0 1]);
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    saveas(Fig,horzcat(filename,'-l1DiagRec'),'png');
    
    Fig = figure();
    imagesc(Znaive/settings.maxValue, [0 1]);
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    saveas(Fig,horzcat(filename,'-naive'),'png');
    
    close all
    save(horzcat(filename,'-results.mat'))
    
    %% create video:
    %         if( exist('settings','var')==1 && isfield(settings,'writerObjGrad') && isempty(settings.writerObjGrad)==0 )
    %         drawnow
    %         pause(3)
    %         writeVideo(settings.writerObjGrad, getframe(FigGrad) );
    %     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zGT_norm2 = norm(zGT); zGT_norm1 = norm(zGT,1); zGT_normInf = norm(zGT,Inf);
results = struct('settings',settings,...
    'e_z',e_z,'e_naive',e_naive,'e_z_diag',e_z_diag,'e_zFast',e_zFast,'e_zFast_diag',e_zFast_diag, ...
    'e_z_l1',e_z_l1,'e_naive_l1',e_naive_l1,'e_z_diag_l1',e_z_diag_l1,'e_zFast_l1',e_zFast_l1,'e_zFast_diag_l1',e_zFast_diag_l1, ...
    'e_z_inf',e_z_inf,'e_naive_inf',e_naive_inf,'e_z_diag_inf',e_z_diag_inf,'e_zFast_inf',e_zFast_inf,'e_zFast_diag_inf',e_zFast_diag_inf, ...
    'cost_z',cost_z,'cost_zFast',cost_zFast,'cost_zFast_diag',cost_zFast_diag,...
    'cost_zNaive',cost_zNaive,'cost_z_diag',cost_z_diag,'cost_zGT',cost_zGT,'cost_zGT_diag',cost_zGT_diag, ...
    'z_l0',z_l0,'zNaive_l0',zNaive_l0,'zGT_l0',zGT_l0,'zFast_l0',zFast_l0,'zFast_diag_l0',zFast_diag_l0,...
    'z_diag_l0',z_diag_l0, 'zNaive_diag_l0',zNaive_diag_l0,'zGT_diag_l0',zGT_diag_l0,...
    'zGT_norm2',zGT_norm2,'zGT_norm1',zGT_norm1,'zGT_normInf',zGT_normInf,...
    'cvx_time',cvx_time,'naive_time',naive_time,'cvx_time_diag',cvx_time_diag,...
    'timeFast',timeFast,'iterFast',iterFast,'timeFast_diag',timeFast_diag,'iterFast_diag',iterFast_diag, ...
    'K',K,'N',N,'percSamples',percSamples,'epsilonK',epsilonK,'noiseNorm',noiseNorm,'noiseInfNorm',noiseInfNorm);
disp('- saved results')
disp('=============================================')
 
   

