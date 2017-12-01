% Minimize total variations of first derivative of the signal, given signal measurements
% Author: Luca Carlone
% date: 2015-11-20
function [results] = example_TV2_2D_theory_temporal_function(settings,tempDataset)

format long
try
    cvx_solver mosek 
    cvx_save_prefs
end

%% create true signal
ZGT = tempDataset(end).ZGT; % current frame
zGT = tempDataset(end).zGT; % current frame
rgb = tempDataset(end).rgb; % current frame
Nx = settings.Nx;
Ny = settings.Ny;
N = settings.N;
tol = settings.tol;
epsilon = settings.epsilon;
isDebug = settings.isDebug;

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
samples = []; R = []; y = []; epsilonVec = [];
horizon = length(tempDataset);
for t=1:horizon % we collect samples from each frame in the time horizon
  samples = [samples tempDataset(t).samples];  
  R = [R; tempDataset(t).R];
  
  %% only for testing: consider the case of no odometry noise
  yt_test = tempDataset(t).y + (tempDataset(t).positionsTrue_z - tempDataset(horizon).positionsTrue_z);
  if norm(yt_test - zGT(tempDataset(t).samples),Inf)> epsilon
      error('inconsistent epsilon')
  end
  
  %% actual measurements: motion is noisy here
  yt = tempDataset(t).y + (tempDataset(t).positionsNoisy_z - tempDataset(horizon).positionsNoisy_z); 
  poseNoise = (horizon-t)*tempDataset(t).deltaTranEps; % each step adds up error
  % noise is sum of sensor noise (epsilon) and motion noise (poseNoise)
  if norm(yt - zGT(tempDataset(t).samples),Inf)> epsilon+poseNoise % if error is larger than predicted
      error('inconsistent epsilon + poseNoise')
  end  
  epsilonVec = [epsilonVec; (epsilon+poseNoise) * ones(size(yt))];   
  y = [y; yt];
end
for k=1:length(y)
    if abs(y(k) - R(k,:) * zGT) >  epsilonVec(k)
       y(k), R(k,:) * zGT
       error('zGT is infeasible - inconsistent noise model')
    end
end
K = length(y); 
percSamples = K/N; 
fprintf('Percentage of samples: %g\n',percSamples)
noiseNorm = norm(tempDataset(t).noise); 
noiseInfNorm = norm(tempDataset(t).noise, Inf);
epsilonK = -1;

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

%% l1 minimization
R = sparse(R);
[z,cvx_time1] = minimizeL1cvx_temp(TV2,N,R,y,epsilonVec);
Z = reshape(z,Nx,Ny);
if abs( norm(TV2 * z,1) - ( norm( vec(Z * H), 1) + norm( vec(V * Z), 1) ) ) >1e-5
    error('mismatch between matrix and vectorized objective')
end

%% Naive Approach: connect the dots
disp('Naive Approach.')
naive_time_temp = tic;
Fun = scatteredInterpolant(X_sample, Y_sample, y, 'linear');
Znaive = Fun(Xq, Yq);
naive_time = toc(naive_time_temp);
zNaive = vec(Znaive);

%% l1 minimization with diagonal terms in the objective
if settings.testDiag
   disp('****************** using diagonal terms! *******************');
   [Hd,Vd] = createFiniteDiffDiag(Nx,Ny);
   % vec(ABC) = (C' \kron A) vec(B)
   % Vd * X * Hd
   TV2_xy = kron(Hd',Vd); 
   TV2_D = sparse([TV2_V; TV2_H; TV2_xy]);
   [z_diag,cvx_time_diag] = minimizeL1cvx_temp(TV2_D,N,R,y,epsilonVec);
   Z_diag = reshape(z_diag,Nx,Ny);
   cost_z_diag = norm(TV2_D * z_diag,1);
   [e_z_diag,e_z_diag_l1,e_z_diag_sm,e_z_diag_inf] = computeError(Z_diag, ZGT); 
   cost_zGT_diag = norm(TV2_D * zGT,1);
   z_diag_l0 = getl0norm(TV2_D * z_diag);
   zNaive_diag_l0 = getl0norm(TV2_D * zNaive);
   zGT_diag_l0 = getl0norm(TV2_D * zGT);
else
    cvx_time_diag = NaN; cost_z_diag = NaN; e_z_diag = NaN;
    e_z_diag_l1 = NaN; e_z_diag_sm = NaN; e_z_diag_inf = NaN; cost_zGT_diag = NaN;
    z_diag_l0 = NaN; zNaive_diag_l0 = NaN; zGT_diag_l0 = NaN;
end

disp('=============================================')
disp('objective value')
cost_zGT = norm(TV2 * zGT,1); fprintf('cost_zGT = %g\n',cost_zGT);
cost_z = norm(TV2 * z,1); fprintf('cost_z = %g\n',cost_z);
cost_zn = norm(TV2 * zNaive,1); fprintf('cost_zNaive = %g\n',cost_zn);
if(settings.testDiag) fprintf('cost_z_diag = %g\n',cost_z_diag); end
cost_diff = cost_zGT - cost_z; fprintf('cost_diff = %g\n',cost_diff);

disp('=============================================')
disp('reconstruction error')
[e_z , e_z_l1, e_z_sm, e_z_inf] = computeError(Z, ZGT); fprintf('e_z = %g, e_z_l1 = %g, e_z_sm = %g, e_z_inf = %g\n',e_z,e_z_l1,e_z_sm,e_z_inf);
[e_n , e_n_l1, e_n_sm, e_n_inf] = computeError(Znaive, ZGT); fprintf('e_n = %g, e_n_l1 = %g, e_n_sm = %g, e_n_inf = %g\n',e_n,e_n_l1,e_n_sm,e_n_inf);
if(settings.testDiag) fprintf('e_z_diag = %g, e_z_diag_l1 = %g, e_z_diag_sm = %g, e_z_diag_inf = %g\n',e_z_diag,e_z_diag_l1,e_z_diag_sm,e_z_diag_inf); end
z_l0 = getl0norm(TV2 * z); 
zNaive_l0 = getl0norm(TV2 * zNaive);
zGT_l0 = getl0norm(TV2 * zGT); 

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
subplot(231); imshow(ZGT); title('img');
subplot(232); imshow(Z);  title('l1') 
subplot(233); 
if e_z > 1e-7
    imshow(zeroOne(Z-ZGT));
else
    imshow(zeros(size(ZGT)));
end
title(sprintf('l1 error (l2:%2.2g)',e_z)) 
subplot(234); imshow(mask); hold on; title(['samples: ',num2str(percSamples)]);
plot(X_sample,Y_sample,'xb','markersize',5)
subplot(235); imshow(Znaive); title('naive')
subplot(236); 
if e_n > 1e-7
    imshow(zeroOne(Znaive-ZGT)); 
else
    imshow(zeros(size(ZGT)));
end
title(sprintf('naive error (l2:%2.2g)',e_n))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(settings.doSaveResults)
    if length(settings.thresVect) > 0 % we manually specify thresholds
        folderName = createFilename(settings,horzcat(settings.resultsFolder,'results_canny_',settings.dataset));
    else
        
        folderName = createFilename(settings,horzcat(settings.resultsFolder,'results_',settings.dataset));
    end
    mkdir(folderName);
    imFolderName = horzcat(folderName,'/images/');
    mkdir(imFolderName);
    filename = horzcat(imFolderName,num2str(settings.imgID));
    
    Fig = figure; 
    hold on
    subplot(221); imshow(ZGT); title(horzcat('full signal'));
    subplot(222); imshow(mask); title(sprintf('samples:%2.2g',percSamples));
    subplot(223); imshow(Z); title(sprintf('l1 rec (error:%2.2g)',e_z))
    subplot(224); imshow(Znaive); title(sprintf('naive (error:%2.2g)',e_n))
    set(Fig, 'paperunits', 'points' )
    set(Fig,'papersize',[300 300]);
    set(Fig,'PaperPositionMode','Auto')
    saveas(Fig,horzcat(filename,'-overview'),'pdf');  
    
    Fig = figure();
    imagesc(ZGT, [0 1]); 
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    saveas(Fig,horzcat(filename,'-original'),'png');
 
    Fig = figure();
    imagesc(mask, [0 1]);
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    saveas(Fig,horzcat(filename,'-samples'),'png');

    Fig = figure();
    imagesc(Z, [0 1]);
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    saveas(Fig,horzcat(filename,'-l1rec'),'png');

    Fig = figure();
    imagesc(Znaive, [0 1]);
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    saveas(Fig,horzcat(filename,'-naive'),'png');
    
    % imwrite(ZGT,horzcat(filename,'-original.png'))
    % imwrite(mask,horzcat(filename,'-samples.png'))
    % imwrite(Z,horzcat(filename,'-l1rec.png'))
    % imwrite(Znaive,horzcat(filename,'-naive.png'))
    
    close all
    save(horzcat(filename,'-results.mat'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zGT_norm2 = norm(zGT); zGT_norm1 = norm(zGT,1); zGT_normInf = norm(zGT,Inf);
results = struct('settings',settings,...
    'e_z',e_z,'e_n',e_n,'e_z_diag',e_z_diag,...
    'e_z_l1',e_z_l1,'e_n_l1',e_n_l1,'e_z_diag_l1',e_z_diag_l1, ...
    'e_z_inf',e_z_inf,'e_n_inf',e_n_inf,'e_z_diag_inf',e_z_diag_inf, ...
    'cost_z',cost_z,'cost_zn',cost_zn,'cost_z_diag',cost_z_diag,'cost_zGT',cost_zGT,'cost_zGT_diag',cost_zGT_diag, ...
    'z_l0',z_l0,'zNaive_l0',zNaive_l0,'zGT_l0',zGT_l0,...
    'z_diag_l0',z_diag_l0, 'zNaive_diag_l0',zNaive_diag_l0,'zGT_diag_l0',zGT_diag_l0,...
    'zGT_norm2',zGT_norm2,'zGT_norm1',zGT_norm1,'zGT_normInf',zGT_normInf,...
    'cvx_time1',cvx_time1,'naive_time',naive_time,'cvx_time_diag',cvx_time_diag,...
    'K',K,'N',N,'percSamples',percSamples,'epsilonK',epsilonK,'noiseNorm',noiseNorm,'noiseInfNorm',noiseInfNorm);
disp('- saved results')
