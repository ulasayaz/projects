% See problem 4 in the following link for more information regarding
% triangular functions as basis for piece-wise linear functions.
% 
% http://www.math.ubc.ca/~malabika/teaching/ubc/spring08/math421-510/hw4.pdf
%
% Let N = 2^m
% We would like to do compressive sensing on x = R * F_inv * z, where 
%   R is the sampling matrix,
%   F_inv is the (inverse) discrete fourier transform matrix, 
%   z is the sparse representation, 
%   x is our piece-wise linear function with only a few critical points.

clear;
clc;
close all;

addpath('../l1magic/Optimization')
addpath('../../../spgl1-1.9')

m = 10;
N = 2^m;
K = 80 % round(N / 4); % number of measurements
C = 5;  % number of critical points in the piecewise linear function


%% Create piecewise linear function x
maxValue = 5;
% x = [0, maxValue * [N-2 : -1 : 0] / (N-1)]';

% define a few quiry points (critical points)
xq = [0, randperm(N-1, C)/N, 1];
yq = (rand(1, C+2) - 0.5) * maxValue;

% generate linear interpolation
u = linspace(0,1,N);
x = interp1(xq,yq,u,'linear')';

% smooth out the function (take only the 20 largest frequency components)
fft_x = fft(x);
sorted = sort(abs(fft_x), 'descend');
fft_x = (abs(fft_x) > sorted(20)) .* fft_x;
x = ifft(fft_x);

figure(1); subplot(221); plot(linspace(0,1,N), x, 'LineWidth', 5); hold on; 

%% Take random samples
% R is the sampling matrix that takes K elements randomly 
R = eye(N);
samples = [1 randperm(N, K) N];
R = R(samples, :);

% observations
y = R * x;
subplot(221); plot( (samples-1)/(N-1), y, '*r', 'MarkerSize', 10);
title('Original Signal x and Sampled Observation y');

%% Naive Approach: connect the dots
x_connect = interp1((samples-1)/(N-1), y, u,'linear')';
subplot(222); plot(linspace(0,1,N), x_connect, 'r', 'LineWidth', 5); 
title(['Naive Approach: error=', num2str(norm(x_connect - x,2))])

%% Create Fourier basis
F = dftmtx(N);
F_inv = conj(dftmtx(N))/N;

%% Direct reconstruction 
z = fft(x); % z = F * x;
maxComponent = max(z(:))
norm0 = sum(abs(z)>0.01*max(abs(z)))
subplot(223); plot(0:1/(N-1):1, F_inv * z, 'b', 'LineWidth', 4); 
title('Basis Verification')

%% Compressive sensing (recovery of z with l-1 optimization)
%% l1-minimization
A = R * F_inv;
epsilon = 1e-5;

%% Using SPGL
opts = spgSetParms('verbosity',1);
% z_vec = spg_bp(A, y, opts);
z_vec = spg_bpdn(A, y, epsilon);
% z_vec = spg_lasso(A, y, norm0);


%%  Using l1-magic
% x0 = A'*y;
% z_vec = l1eq_pd(x0, A, [], y, epsilon);

%% Using cvx
% cvx_begin
% variable z_vec(N,1) complex;
% minimize( norm(z_vec,1) )
% subject to
%     %norm(y - A * z_vec, 2) <= epsilon;
%     y == A * z_vec
% cvx_end

% reconstruction
reconstruction = F_inv * z_vec;
subplot(224); plot(linspace(0,1,N), reconstruction, 'green', 'LineWidth', 4); 
title(['l1 Reconstruction (Measurements ', num2str(K), '/', num2str(N), '), error=', num2str(norm(reconstruction-x,2))])
% legend('Piecewise linear function x', 'Sampled observations', 'Reconstructed signal')


%% Visualize sparsity level
norm0
N
% figure(2); hist(abs(z), 20); title('Histogram of Component Values')
figure(2);
subplot(211);plot(abs(fft(x)));title('FFT of the input signal')
subplot(212);plot(abs(z_vec));title('FFT of the recovered signal')