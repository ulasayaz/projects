clear all
close all
clc

ADD_LIBRARIES

%% settings l1-inf
settings.percSamples = 1; 
settings.optMode = 'l1inf'; 
settings.lambda_reg = -1;
settings.noiseMode = 'linf';
settings.doAddNeighbors = 1;

%% common settings
settings.nrCorners = 3;
settings.maxValue = 5;
settings.addNoise = 0;
settings.mode = 'pwLinear'; % pwLinear, toy, saddle, stripe
settings.pwOption = 'diagonal'; % pyramid, diagonal
settings.d = 20; % must be odd for testing
settings.sampleMode = 'edgesOnly';
settings.h = -1;
settings.testDiag = 1;
settings.isDebug = 1;
settings.doAddBoundary = 0;
settings.doAddRandom = 0;
settings.nrNeighbors = 1;
settings.includeAllDirs = 0; % 0 or 1
settings.tol = 1e-5;
settings.dataset = '';
settings.path = '';
settings.imgID = -1;
settings.subSample = -1;
settings.useInputImage = 0;
settings.doSaveResults = 0;
settings.epsilon = 0;

if isfield(settings,'isBounded')==0
    settings.isBounded = 0;
    warning('isBounded undefined in original settings')
end
if isfield(settings,'mu')==0
    settings.mu = 0.001;
    warning('mu undefined in original settings')
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
%[areas,labels] = segmentBinaryImageAndGetAreas(normGrad_GT);
labels = 0;

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

%% check conditions in [KabanavaRauhut Cosparsity Survey and Vaiter et. al. "Decomposable priors" 2013]

zGTdiff = TV2 * zGT;

codim = size(TV2,1);

supportSet = sort(find(abs(zGTdiff) > tol));
cosupportSet = setdiff(1:codim,supportSet);

s = numel(supportSet);

% Msupp = zeros(s,codim); % restricts a vector to supportSet
% Msupp(:,supportSet) = eye(s);
% Mcosupp = zeros(codim-s);


suppSigns = sign(zGTdiff(supportSet));

% if sort(supportSet') ~= sort(supp-1)
%     error('support sets do not match!');
% end

D = TV2';
samplesComp = setdiff(1:N,samples);

Dsamp = D(samples,:);
DsampCo = D(samplesComp,:);

DstarJ = TV2(cosupportSet,:);
% next line takes a while if d is large
DJrow = orth(full(DstarJ'));

DstarJnull = null(full(DstarJ));

% calculate constants 
H = full(R * DstarJnull);
sing1 = sqrt(eig( H' * H));
Ca = min(sing1);

F = full(DstarJ * DJrow);
sing2 = sqrt(eig( F' * F));
Cdj = min(sing2);


if isempty(null(full(DsampCo)))
    error('Null space of DsampCo is empty!');
end

%cvx_solver  mosek % sedumi %
%cvx_save_prefs

cvx_begin
cvx_quiet true
cvx_precision best

variable u(codim,1);
%minimize( norm( Dsamp * u ) )
minimize norm(u(cosupportSet),Inf)

subject to
norm(DsampCo * u) <= tol;
u(supportSet) == suppSigns;
norm(u(cosupportSet),Inf) <= 1

cvx_end

kappa = norm(u(cosupportSet),Inf);
eta = Dsamp * u;

figure(2)
plot(abs(u))

onesInU = sort(find(abs(u) > 1 - tol))

if kappa > 1
    warning('dual vector fails, i.e., kappa > 1')
end

fprintf('norm of ||eta|| = %g \n', norm( eta ))
fprintf('norm of A = %g \n', norm( full(R) ))
fprintf('kappa = %g \n', kappa)
fprintf('Ca = %g \n', Ca)
fprintf('Cdj = %g \n', Cdj)

% constant in the error term in Thm 5 (KabanavaRauhut)

C = (1/Ca) + ( Ca + norm(full(R) ) ) * norm(eta) / Cdj /Ca / (1 - kappa )





