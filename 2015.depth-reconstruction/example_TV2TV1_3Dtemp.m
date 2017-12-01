function [errorL1last,errorNaive] = example_TV2TV1_3Dtemp(settings)
% Minimize total variations of first derivative of the signal, given signal measurements

addpath('../data/traindata')
% addpath('~/Dropbox/3.FUTURE_WORK/sparseSensing/ZED_stereo_data_Nov17/depth_scaled')
addpath('~/Dropbox/3.FUTURE_WORK/sparseSensing/ZED_stereo_data_fixed_camera/depth_scaled/')

addpath('../myLib')
addpath('../plane_line_intersect')

try
    cvx_solver  mosek % sedumi %
    cvx_save_prefs
end

if nargin < 1
    clear all
    close all
    clc
    warning('using default arguments')
    percSamples = 0.05; % percentage of measured pixels (samples)
    subSample = 0.1;
    lambda = 0; % parameter to weigh the diagonal TV2 norms
    tau = 0.5; % parameter to weigh the temporal TV norm
    epsilon = 1e-2;
    addNoise = 0;
    mode = 'video' % load, planes, pwLinear
    optMode = 'standard'; % useEdges, standard, useTVsqrt
    msg = '0001';
    nrFrames = 4;
    doSave = 0;
    regenerateRandom = 1;
else
    percSamples = settings.percSamples; 
    subSample = settings.subSample; 
    lambda = settings.lambda; 
    theta = settings.theta; 
    tau = settings.tau;
    epsilon = settings.epsilon;
    addNoise = settings.addNoise; 
    mode = settings.mode; 
    optMode = settings.optMode; 
    msg = settings.msg;
    nrFrames = settings.nrFrames;
    regenerateRandom = settings.regenerateRandom;
    doSave = 0;
end
repeat = 0; % repeat = 1 uses previous image and sampling points, = 2 load samplesLog.mat

%% Load depth data
imName = '';
switch mode
    case 'video'
        startFrame = str2num(msg)
        [img, Nx, Ny, N] = loadImageSet(startFrame, subSample, nrFrames);      
    otherwise
        error('wrong choice of mode')
end

%% show the signal and temporal differences
dimg = diff(img,1,3); % take the differences in temporal direction

% show the difference for 1 frame
D = dimg(:,:,1); % we look at one of the frames
Ds = sort(abs(D(:)),'descend');
figure(8);
plot(Ds,'-o'); 

% show the *color bar* differences for all frames
figure(102)
spcount = nrFrames;
sd=factor_subplots(spcount); 
subplot(sd(1),sd(2),1);
imshow(img(:,:,1)); title('Original');
for i = 1: nrFrames-1
   subplot(sd(1),sd(2),i + 1); 
   imshow(dimg(:,:,i)); title(['difference of frames ',num2str(i),' and ',num2str(i+1)]);
   colormap('default');        
   colorbar;
end

%% set parameters
K = round(N * percSamples) % nr of measurements per frame

%% Take random measurements: R is the sampling matrix that takes K elements randomly
switch repeat 
    case 0 % create samples
        disp('Creating samples.')
        [M,samples] = generateSamplesSet(nrFrames,img,N,K,addNoise,epsilon,regenerateRandom);       
    case 1
        warning('repeat = 1: Using previous image and sampling points!')
    case 2
        for i = 1: nrFrames  
            R = speye(N);
            load('samplesLog.mat','samples')
            R = double(R(samples, :));
            img_i_vec = vec(img(:,:,i));
            y = img_i_vec(samples) + addNoise * epsilon * ( 2 * rand(K,1)-1);
            M(i).R = R;
            M(i).y = y;
        end
    otherwise
        error('repeat must be 0, 1, or 2')
end

%% l1 minimization
[H,V] = createFiniteDiff2(Nx,Ny); % Create TV matrices
[H1,V1] = createFiniteDiff1(Nx,Ny);
% [VZ_GT, ZH_GT, normGrad_GT, VZ_uncut_GT, ZH_uncut_GT] = computeGrad(img,V,H);
% VZ_GT_mask = 1-zeroOneTh(VZ_uncut_GT); % where there is no gradient, this matrix is 1, otherwise 0
% ZH_GT_mask = 1-zeroOneTh(ZH_uncut_GT);
% normGrad_GT_mask = zeroOneTh(normGrad_GT, 1e-4); % zero where no gradient
% Sh = speye(Nx); Sh = Sh(2:end-1,:);
% Sv = speye(Ny); Sv = Sv(:,2:end-1);

disp('starting CVX')
cvx_begin
variable Z(Nx,Ny,nrFrames);
% variable errSq(nrFrames);
variable norm1H(nrFrames);
variable norm1V(nrFrames);   
variable norm1diag(nrFrames);
minimize( tau * norm( vec( Z(:,:,2:nrFrames) - Z(:,:,1:nrFrames-1) ) , 1) + ... % temporal TV norm
          sum(norm1H) + sum(norm1V)  + ... % second derivatives
          lambda * sum(norm1diag) ) % diagonal derivative
subject to
    for t = 1:nrFrames
        norm(M(t).y - M(t).R * vec(Z(:,:,t))) <= epsilon; % errSq(t);
        norm( vec( Z(:,:,t) * H ) , 1) <= norm1H(t);
        norm( vec( V *  Z(:,:,t)) , 1) <= norm1V(t);
        norm( vec(V1 * Z(:,:,t) * H1), 1) <= norm1diag(t);
    end
    % sum( square ( errSq ) ) <= epsilon * nrFrames;
cvx_end

%% visualize results
figure(5)
sd=factor_subplots(nrFrames+1); 
subplot(sd(1),sd(2),1);
imshow(img(:,:,1)); title('Original');
for i = 1 : nrFrames
   subplot(sd(1),sd(2),i + 1);
   errorL1(i) = computeError(Z(:,:,i),img(:,:,i));
   imshow(Z(:,:,i)); title(sprintf('rec %d (e:%2.2g)',i,errorL1(i)))
end
errorL1last = errorL1(end)

%% Naive Approach: only for the latest image
disp('Naive Approach.')
[Xq, Yq] = meshgrid(1:size(img(:,:,end),2), 1:size(img(:,:,end),1));
X_sample = Xq(samples)';
Y_sample = Yq(samples)';
Fun = scatteredInterpolant(X_sample, Y_sample, M(nrFrames).y, 'linear');
naive_result = Fun(Xq, Yq);
errorNaive = computeError(naive_result,img(:,:,end));

%% compute mask and visualize: only for the latest image
index = zeros(N, 1);
index(samples) = 1;
mask = reshape(index, size(img(:,:,end)));
Fig = figure(10);
subplot(221); imshow(img(:,:,end)); title(horzcat('img-',num2str(str2num(msg)+nrFrames-1)))
subplot(222); imshow(mask); title(horzcat('%samples-',num2str(percSamples)));
subplot(223); imshow(Z(:,:,end));  title(sprintf('l1 rec (err:%2.2g)',errorL1(end)))
subplot(224); imshow(naive_result); title(sprintf('naive (err:%2.2g)',errorNaive))
myAxis = [0.5000   64.5000    0.5000   48.5000];
axis(myAxis)
hold on
if( exist('settings','var')==1 && isfield(settings,'writerObj') && isempty(settings.writerObj)==0 )
    drawnow
    pause(3)  
    set(Fig, 'paperunits', 'points' )
    set(Fig,'position',[0, 0, 300 250]); %
    set(Fig,'papersize',[300 250]);
    set(Fig,'PaperPositionMode','Auto')  
    fframe = getframe(Fig);
    writeVideo(settings.writerObj,fframe);
end

%% save results to file
filename = horzcat('./results/result-',mode,imName,'-sub',num2str(100*subSample),'-percSamples',num2str(100*percSamples),'-eps',num2str(log10(epsilon)),'-test',msg,'-comp')
set(Fig, 'paperunits', 'points' )
set(Fig,'position',[0, 0, 300 300]); %
set(Fig,'papersize',[300 300]);
set(Fig,'PaperPositionMode','Auto')
if(doSave) saveas(Fig,filename,'pdf'); end


