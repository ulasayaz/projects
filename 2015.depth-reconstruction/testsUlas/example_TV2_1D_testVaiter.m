% Test robust recovery error in 1D signals depending on the various sampling patterns
% Author: Ulas Ayaz
% date: 2015-11-20

clear;
clc;
close all;
format long
addpath('../l1magic/Optimization')
addpath('../grammSchmidt')
addpath('../myLib')

N = 2^10 + 1; % total number of points number of points
alpha = 0.01;
optMode = 'standard'; % standard, denoising, peaking
percSamples = 0.1; % percentage of measured pixels (samples)
K = round(N * percSamples) % number of measurements
repeat_samples = 0;
repeat_img = 0;
nrCorners = 3;  % number of critical points in the piecewise linear function
epsilon = 0.01;
sampleOnlyCorners = 0;
sampleAlsoCorners = 0;
sampleAlsoNeighbors = 1;

tol = 1e-7;

%% Create piecewise linear function x (true signal)
maxValue = 2;
if repeat_img == 0
    % define a few query points (critical points)
    supp = randsample(2:N-1, nrCorners);
    xq = [1, supp, N];
    yq = (rand(1, nrCorners+2) - 0.5) * maxValue;
end
% generate linear interpolation
zGT = interp1(xq,yq,1:N,'linear')';

% let minimum be above zero
minZGT = min(zGT);
zGT = zGT - minZGT + 0.1;
yq = yq - minZGT + 0.1;

%% Compressive sensing (recovery of z with l-1 optimization)
% R is the sampling matrix that takes K elements randomly
Rfull = eye(N);
if repeat_samples == 0
    %samples = randperm(N, K);
    if sampleOnlyCorners == 0 && sampleAlsoCorners == 0 && sampleAlsoNeighbors == 0
        [samples,samplesExtreme] = distributeSamples(N,K,supp);
        samples = unique(samples);
    elseif sampleOnlyCorners == 1
        samples = xq;
        samples = unique(samples);
        samplesExtreme = samples;
    elseif sampleAlsoCorners == 1
        xq_noCorners = setdiff([1:N], xq);
        N_noCorners = length(xq_noCorners);
        samples = [xq, randsample(xq_noCorners,K)]; % not rsample(N,K)
        samples = unique(samples);
        samplesExtreme = samples;
    elseif sampleAlsoNeighbors == 1
        samples1 = xq;
        samples2 = samples1 + 1; % point after corner
        samples3 = samples1 - 1; % point after corner
        samples = [samples1 samples3(2:end) samples2(1:end-1)];
        samples = unique(samples);
        samplesExtreme = samples;
    end
    %%% choose samples at the critical points and next points
    % samples = setdiff( xq , N+1 ); % possibly exclude N+1
    % samples = union( samples, samplesExtreme );
    R = Rfull(samples, :);
    Rext = Rfull(samplesExtreme, :);
end
RGauss = randn(K,N) / sqrt(K); % Gaussian

% observations
y = R * zGT;
yext = Rext * zGT;
yG = RGauss * zGT;
K = length(y);

%% plot the signal and the sampling points
Fig  = figure('Position', [1000, 1000, 800, 400]);
subplot(121); 
plot(zGT, 'LineWidth', 5); hold on;
plot( samples, zGT(samples), '*r', 'MarkerSize', 10)
grid minor
subplot(122); 
plot(zGT, 'LineWidth', 5); hold on;
plot( samplesExtreme, zGT(samplesExtreme), '*k', 'MarkerSize', 10)
grid minor

%% create TV2 matrix
[TV2] = createFiniteDiff2(N);

zGTdiff = TV2 * zGT;

supportSet = sort(find(abs(zGTdiff) > tol));
cosupportSet = setdiff(1:N-2,supportSet);

s = numel(supportSet);

Msupp = zeros(s,N-2); % restricts a vector to supportSet
Msupp(:,supportSet) = eye(s);

%Mcosupp = zeros(N-2-s);


suppSigns = sign(zGTdiff(supportSet));

if sort(supportSet') ~= sort(supp-1)
    error('support sets do not match!');
end

%% check conditions in [KabanavaRauhut Cosparsity Survey and Vaiter et. al. "Decomposable priors" 2013]

D = TV2';
samplesComp = setdiff(1:N,samples);

Dsamp = D(samples,:);
DsampCo = D(samplesComp,:);

% if min(svd(full(DsampCo))) < tol
%     error('Null space of Dsamp is empty!');
% end

ker_D = null(full(DsampCo));
v = randn(size(ker_D,2),1);
v = v/norm(v);
t = ker_D * v;

cvx_begin
cvx_precision best
variable u(N-2,1);
%minimize( norm( Dsamp * u ) )
minimize norm(u(cosupportSet),Inf)

subject to
norm(DsampCo * u) < tol;
u(supportSet) == suppSigns;
norm(u(cosupportSet),Inf) < 1
%abs(u(cosupportSet)) < 1;
%abs(u(cosupportSet)) >= 1;
cvx_end

norm( Dsamp * u )

%% l1 minimization
% run with full samples
switch optMode
    case 'standard'
        cvx_begin
        cvx_precision best
        variable z(N,1);
        minimize( norm( TV2 * z,1) )
        subject to
        %y == R * z;
        norm(y - R * z) <= epsilon;
        cvx_end
    case 'denoising'
        cvx_begin
        cvx_precision best
        variable z(N,1);
        variable d;
        minimize( 0.5 * d^2 + alpha * norm( TV2 * z,1) )
        subject to
        norm(y - R * z) <= d;
        cvx_end 
    case 'peaking'        
        minCurv = 1e-3;
        cvx_begin
        cvx_precision best
        variable z(N,1);
        variable b_l1(N-2,1);
        variable b_curv(N-2,1);
        variable d;
        minimize( sum(b_l1) + 0.0001 * sum(b_curv) )
        subject to
        y == R * z;
        abs( TV2 * z) <= b_l1;
        b_l1 >= minCurv * ones(N-2,1) - b_curv;
        b_l1 >= 0;
        b_curv >= 0;
        cvx_end
end

%%%%%%%%%%%%%%%%%%%%




































