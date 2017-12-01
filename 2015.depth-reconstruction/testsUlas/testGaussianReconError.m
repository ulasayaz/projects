% tests reconstruction error vs noise level with Gaussian measurements


try
    cvx_solver  mosek % sedumi %
    cvx_save_prefs
end


close all
clc
%percSamples = 0.03; % percentage of measured pixels (samples)

lambda = 0; % parameter to weigh the diagonal TV2 norms
theta = 0; % parameter to weigh the diagonal TV3 norms
beta = 0; % parameter to weigh the TV norm

addNoise = 1; % decides to add noise
mode = 'pwLinear' % load_(png,jpg,pgm), planes, pwLinear
optMode = 'standard'; % useEdges, standard, denoising, useTVsqrt, matrixCompletion
optMeasure = 'Gaussian'; % uniform, Gaussian
nrCorners = 4 ;  % select the number of corners
repeat_img = 0; % repeat_img = 1 uses previous image and sampling points
repeat_samples = 0;

d = 30;
img = createPWlinear(nrCorners, d, d);

eps = [0:0.01:0.12];
T = numel(eps);

trial = 20;

result = zeros(T,1);


tol = 1e-5;


%% set parameters
Nx = size(img,1);
Ny = size(img,2);
N = Nx * Ny;
zGT = vec(img); %vector version of the image
K = round(N * percSamples); % nr of measurements
fprintf('Nx=%d, Ny=%d, N=%d, K=%d\n',Nx,Ny,N,K)

surf(1:Nx, 1:Ny, img); hold on; shading faceted % flat, axis equal;

%% Create TV matrices
[H,V] = createFiniteDiff2(Nx,Ny);

[VZ_GT, ZH_GT, normGrad_GT, VZ_uncut_GT, ZH_uncut_GT] = computeGrad(img,V,H);

figure
subplot(211); imshow(zeroOneTh(abs(VZ_GT))); title('V * img');
subplot(212); imshow(zeroOneTh(abs(ZH_GT)));        title('img * H')

%imshow(zeroOneTh(normGrad_GT));  title('img * grad');

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

[Ix,Jx,DIx,DJx,sIx] = getVaiterMat(TV2,zGT,tol);

s = numel(Ix) % sparsity

%nnz = getl0norm(img*H) + getl0norm(V*img)

%% kappa (theory)

kappa = 1 - (3/16/pi) * (2 - 2 * (Nx+Ny)/N - s/N)^2

percSamples = 0.4; % percentage of measured pixels (samples)
K = round(N * percSamples); % nr of measurements

for i = 1:T
    i
    
    epsilon = eps(i);
    total = 0;
    for j = 1: trial
        if mod(j,5) == 0
            j
        end
        %% Gaussian measurements
        error = randn(K,1);
        error = error/norm(error);
        
        R = randn(K,N) / sqrt(K);
        y = R * zGT + addNoise * epsilon * error;
        
        %% L1 optimization
        
        cvx_begin
        variable Z(Nx,Ny);
        minimize( norm( vec(Z * H), 1) + norm( vec(V * Z), 1) )
        %lambda * norm( vec(V1 * Z * H1), 1)  + ...
        %theta * (norm( vec(Z * H1 * H2 * H3), 1) + norm( vec(V3 * V2 * V1 * Z), 1)) + ...
        %beta * ( norm( vec(Z * H1), 1) + norm( vec(V1 * Z), 1) ) )
        subject to
        % y == R * z;
        norm(y - R * vec(Z)) <= epsilon;
        cvx_end
        
        total = total + norm(img - Z);
    end
    
    result(i) = total/trial;
    
end

figure(1)
plot(eps,result,'-o')

















