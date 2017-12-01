clear all
close all
clc

format long
addpath('../myLib')
% addpath('../plane_line_intersect')
try
    cvx_solver  mosek % sedumi %
    cvx_save_prefs
end

%% parameters
d = 20; % must be odd for testing
nrCorners = 2;  % select the number of corners
epsilon = 0;
addNoise = 0; % decides to add noise

recoveryMode = 'GAP'; % L1 , GAP
sampleMode = 'uniform'; % imgCornersOnly, cornersOnly, edgesOnly, h_per_area, uniform, square_patch
side = 40;
twin = 0;
offset = 20;
h = 3;
percSamples = 0.2;
doAddBoundary = 0;
doAddNeighbors = 0;

%% fixed parameters
mode = 'pwLinear';
pwOption = 'diagonal'; % pyramid, diagonal
if strcmp(mode,'pwLinear') == 0
    error('doSampleCorners ~= 1 with real images not implemented yet')
end
ZGT = createPWlinear(nrCorners, d, d);
%ZGT = createPyramid(d,pwOption);
%close;

Nx = size(ZGT,1);
Ny = size(ZGT,2);
N = Nx * Ny;
zGT = vec(ZGT); %vector version of the image
K = round(N * percSamples); % nr of measurements
fprintf('Nx=%d, Ny=%d, N=%d, K=%d\n',Nx,Ny,N,K)

%% Create TV matrices
[H,V] = createFiniteDiff2(Nx,Ny);
[H1,V1] = createFiniteDiff1(Nx,Ny);
[H2,V2] = createFiniteDiff1(Nx-1,Ny-1);
[H3,V3] = createFiniteDiff1(Nx-2,Ny-2);
[VZ_GT, ZH_GT, normGrad_GT, VZ_uncut_GT, ZH_uncut_GT] = computeGrad(ZGT,V,H);
[areas,labels] = segmentBinaryImageAndGetAreas(normGrad_GT)
close;
sparsity = sum(zeroOneTh(normGrad_GT(:)))

tol = 1e-7;

% show laplacians of the original 
figure;
set(gcf,'name','original laplacians')
subplot(131); imshow(zeroOneTh(abs(VZ_uncut_GT),tol));   title('V * ZGT');
subplot(132); imshow(zeroOneTh(abs(ZH_uncut_GT),tol));   title('ZGT * H')
subplot(133); imshow(zeroOneTh(normGrad_GT,tol));  title('ZGT * grad')
%% Take samples
disp('Creating samples.')
Rfull = speye(N);
Id = speye(N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch sampleMode
    case 'imgCornersOnly' % 4 image corners only
        samples = [1  Nx  Nx*Ny  Nx*(Ny-1)+1];
    case 'cornersOnly'
        pwOption
        if strcmp(pwOption,'diagonal')
            samples = [1  Nx  Nx*Ny  Nx*(Ny-1)+1]; % corners
        else
            if (rem(Nx,2)>0) warning('this test requires Nx=Ny odd'); end
            samples = [1  Nx  Nx*Ny  Nx*(Ny-1)+1  Nx*floor(Ny/2) + ceil(Nx/2)]; % corners
        end
    case 'edgesOnly'
        zNx = zeros(Nx,1);
        HDiff = [zNx ZGT * H zNx];
        hcorners = find(abs(vec(HDiff)) > tol);
        zNy = zeros(Ny,1);
        VDiff = [zNy'; V * ZGT; zNy'];
        vcorners = find(abs(vec(VDiff)) > tol);
        samples = union(hcorners, vcorners);
        samples = samples(:)'; % row vector
    case 'h_per_area'
        disp('h_per_area')
        samples = distributeSamples2D(Nx, Ny, h,labels)
    case 'square_patch'
        disp('square_patch')
        samples = createPatches(Nx,Ny,side,offset ,twin);
    case 'uniform'
        warning('uniform sampling')
        samples = randperm(N, K);
    otherwise
        error('sampleMode selection was wrong')
end
if(doAddBoundary==1) samples = addBoundary(samples,Nx,Ny); end
if(doAddNeighbors==1) samples = addNeighbors(samples,Nx,Ny,1,1); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if K ~= length(samples) 
    percSamples = -1;
    K = length(samples);
end
y = zGT(samples) + addNoise * epsilon * ( 2 * rand(K,1)-1);
R = double(Rfull(samples, :));

%% plot samples on the original 3d signal
figure(101); clf
% [X_sample, Y_sample] = ind2sub([Nx Ny],samples);
[Xq, Yq] = meshgrid(1:size(ZGT,2), 1:size(ZGT,1));
X_sample = Xq(samples)';
Y_sample = Yq(samples)';
surf(1:Nx, 1:Ny, ZGT); hold on; shading faceted % flat, axis equal;
plot3(X_sample, Y_sample, y, 'or', 'markersize', 10)
set(gcf,'name','original with samples')

%% check theory
tol = 1e-5;
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
TV2 = [TV2_V; TV2_H];

% get sparsity
[Ix,Jx,DIx,DJx,sIx] = getVaiterMat(TV2,zGT,tol);
s = numel(Ix); % sparsity
nnz = getl0norm(ZGT*H) + getl0norm(V*ZGT);

p = size(TV2,1);


%% recovery
rec_tmp = tic;
switch recoveryMode
    case 'L1'
        disp('starting CVX')
        cvx_begin
        cvx_precision best
        variable Z(Nx,Ny);
        minimize( norm( vec(Z * H), 1) + norm( vec(V * Z), 1) )
        subject to
        norm(y - R * vec(Z)) <= epsilon;
        cvx_end
        [errorTVl2,errorTVl1] = computeError(Z,ZGT);
        z = vec(Z);
        
        if abs( norm(TV2 * z,1) - ( norm( vec(Z * H), 1) + norm( vec(V * Z), 1) ) ) > 1e-6
            error('mismatch between matrix and vectorized objective')
        end
        
    case 'GAP'
        k = 0;
        t = 1; % in (0,1]
        J = [1:p]; % initial cosupport
        sigma = 1e-3; % stopping precision
               
        z = solveLeastSquares(R,TV2(J,:),y); % initial solution
        z0 = z;
        % lambda = 1e-5; % in (1e-4,1e-6)
        % z = solveLeastSquares(R,TV2(J,:),y,lambda);
        maxIter = N;
        % residual = norm(y - R * z);
        for i = 1:maxIter
            if mod(i,10) == 0
                i
            end
            alpha = TV2 * z;
            Jdelete = sort( find( abs(alpha) >= t* max(abs(alpha)) ) ); % find largest entries
            J = setdiff(J,Jdelete); % update cosupport
            z = solveLeastSquares(R,TV2(J,:),y); %update solution
            % residual = norm( y - R * z);
            if norm(z - z0) < sigma
                %disp('breaking the while loop: solution is static.')
                % break;
            end
            z0 = z;     
        end  
        fprintf('nr iterations: %d\n',i)
        Z = reshape(z,Nx,Ny);
        [errorTVl2,errorTVl1] = computeError(Z,ZGT);
end
rec_time = toc(rec_tmp)

figure
plot(z,'-or'); hold on
plot(z0,'-ok');
plot(zGT,'-*b');

%% Naive Approach: connect the dots
disp('Naive Approach.')
[Xq, Yq] = meshgrid(1:size(ZGT,2), 1:size(ZGT,1));
X_sample = Xq(samples)';
Y_sample = Yq(samples)';
Fun = scatteredInterpolant(X_sample, Y_sample, y, 'linear');
Znaive = Fun(Xq, Yq);
[errorNaivel2,errorNaivel1] = computeError(Znaive,ZGT);

[X_sample, Y_sample] = ind2sub([Nx Ny],samples);

%% display statistics
fprintf('Error l1 sol: %g (l2), %g (l1)\n',errorTVl2,errorTVl1)
fprintf('Error naive sol: %g (l2), %g (l1)\n',errorNaivel2,errorNaivel1)
fprintf('TV2 l1 sol: %g\n', norm( vec(Z * H), 1) + norm( vec(V * Z), 1) )
fprintf('TV2 GT sol: %g\n', norm( vec(ZGT * H), 1) + norm( vec(V * ZGT), 1) )
fprintf('TV2 naive sol: %g\n', norm( vec(Znaive * H), 1) + norm( vec(V * Znaive), 1) )

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
subplot(221); imshow(ZGT); title(horzcat('img'));
subplot(222); imshow(mask); title(['samples with ',num2str(percSamples)]); hold on
plot(X_sample,Y_sample,'*b', 'MarkerSize', 20)
subplot(223); imshow(Z);  title(sprintf('TV2 rec (error l2:%2.2g, l1:%2.2g)',errorTVl2,errorTVl1))
subplot(224); imshow(Znaive); title(sprintf('naive (error l2:%2.2g, l1:%2.2g)',errorNaivel2, errorNaivel1))

%% compute mask and visualize
Fig = figure(11); clf
colormap hsv
cm = colormap;
hold on; set(gcf,'name','original VS reconstruction')
mesh(1:Nx, 1:Ny, ZGT); % black
mesh(1:Nx, 1:Ny, Z);
view(3)

%% Check the gradients
[VZ, ZH, normGrad, VZ_uncut, ZH_uncut] = computeGrad(Z,V,H);
[VZ_naive, ZH_naive, normGrad_naive] = computeGrad(Znaive,V,H);
Fig = figure(4); set(Fig,'Position',[pos * horz, (height- vert), horz, vert]); pos = pos+1; hold on
subplot(331); imshow(zeroOneTh(abs(VZ_GT)));  title('V * img');
subplot(332); imshow(zeroOneTh(abs(ZH_GT)));  title('img * H')
subplot(333); imshow(zeroOneTh(normGrad_GT)); title('img * grad')
subplot(334); imshow(zeroOneTh(abs(VZ)));  title('V * Z')
subplot(335); imshow(zeroOneTh(abs(ZH)));  title('Z * H')
subplot(336); imshow(zeroOneTh(normGrad));  title('est. grad')
subplot(337); imshow(zeroOneTh(abs(VZ_naive)));  title('V * naive')
subplot(338); imshow(zeroOneTh(abs(ZH_naive)));  title('naive * H')
subplot(339); imshow(zeroOneTh(normGrad_naive));  title('naive grad')

%% test optimality conjecture
close all;
figure
imshow(mask,'InitialMagnification','fit'); hold on;
plot(X_sample,Y_sample,'*b', 'MarkerSize', 10)

figure
subplot(121);
imshow(zeroOneTh(normGrad_GT),'InitialMagnification','fit');  title('original grad');
hold on;
plot(X_sample,Y_sample,'*b', 'MarkerSize', 10)
subplot(122);
imshow(zeroOneTh(normGrad),'InitialMagnification','fit');  title('est. grad'); hold on;
plot(X_sample,Y_sample,'*b', 'MarkerSize', 10)

figure
subplot(121)
surf(1:Nx, 1:Ny, ZGT); hold on; shading faceted % flat, axis equal;
plot3(X_sample, Y_sample, y, 'or', 'markersize', 10)
set(gcf,'name','original with samples')
subplot(122)
surf(1:Nx, 1:Ny, Z); hold on; shading faceted % flat, axis equal;
plot3(X_sample, Y_sample, y, 'or', 'markersize', 10)
set(gcf,'name','recovered with samples')

% show directional laplacians
figure
subplot(221);
surf(ZH); hold; plot(X_sample,Y_sample,'*b', 'MarkerSize', 10); title('recovered ZH');
subplot(222);
surf(VZ); hold; plot(X_sample,Y_sample,'*b', 'MarkerSize', 10); title('recovered VZ');
subplot(223)
surf(ZH_GT); hold; plot(X_sample,Y_sample,'*b', 'MarkerSize', 10); title('original ZH_GT');
subplot(224)
surf(VZ_GT); hold; plot(X_sample,Y_sample,'*b', 'MarkerSize', 10); title('original VZ_GT');

%% surf figures
Fig = figure(8); set(Fig,'Position',[pos * horz, (height-vert), horz, vert]); pos = pos+1;
[Xq, Yq] = meshgrid(1:size(ZGT,2), 1:size(ZGT,1));
subplot(131); surf(Xq,Yq,ZGT); title('original');
subplot(132); surf(Xq,Yq,Z); title('Z');
set(gcf,'name','original with samples')
subplot(133); ERZ = zeroOne(abs(Z-ZGT)); imshow(ERZ);

%% surf figures
Fig = figure; 
subplot(121); surf(Xq,Yq,ZGT); title('original');
subplot(122); surf(Xq,Yq,Z); hold on; title('Z');
plot3(X_sample, Y_sample, Z(samples), 'or', 'markersize', 20)

%% test optimality
tol = 1e-7;
isOpt_z = testOpt(TV2,z,zGT,nonSamples,tol)
if(isOpt_z > tol) error('optimality test does not work - z'); end

[isOpt_zGT,sigm_GT,g_GT] = testOpt(TV2,zGT,zGT,nonSamples,tol);
isOpt_zGT
if(isOpt_zGT > tol) warning('optimality test does not work - zGT'); end

isOpt_rand = testOpt(TV2,rand(N,1),zGT,nonSamples,tol);
if(isOpt_rand < 1) fprintf('isOpt = %g\n',isOpt_rand); error('optimality test is vacuous - rand'); end


