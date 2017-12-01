% See problem 4 in the following link for more information regarding
% triangular functions as basis for piece-wise linear functions.
% 
% http://www.math.ubc.ca/~malabika/teaching/ubc/spring08/math421-510/hw4.pdf
%
% Let N = 2^m
% We would like to do compressive sensing on x = R*H'*z, where 
%   R is the sampling matrix,
%   H is the haar basis, 
%   z is the sparse representation, 
%   x is our piece-wise linear function with only a few critical points.

clear;
clc;
close all;

addpath('../l1magic/Optimization')

m = 10;
N = 2^m;
K = 200 % round(N / 4); % number of measurements
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

figure(1); 
subplot(221); plot(linspace(0,1,N), x, 'LineWidth', 5); hold on; 

%% Take random observations
% R is the sampling matrix that takes K elements randomly 
R = eye(N);
samples = randperm(N, K);
R = R(samples, :);

% observations
y = R * x;
subplot(221); plot( (samples-1)/(N-1), y, '*r', 'MarkerSize', 10)
title('Original Signal x and Sampled Observation y');

%% Naive Approach: connect the dots
x_connect = interp1((samples-1)/(N-1), y, u,'linear')';
subplot(222); plot(linspace(0,1,N), x_connect, 'r', 'LineWidth', 5); 
title('Naive Approach: connect the dots')

%% Create haar basis
H = haar_basis(N);

%% Direct reconstruction 
z = H * x;
norm0 = sum(abs(z)>0.1)
subplot(223); plot(0:1/(N-1):1, H' * z, 'b', 'LineWidth', 4); 
title('Inverse Transform H^T * z (for Basis Verification)')


%% Compressive sensing (recovery of z with l-1 optimization)

% initial guess
z0 = rand(size(z));

% solve the LP
A = R * H';
tic
zp = l1eq_pd(z0, A, [], y, 1e-3);
toc

% reconstruction
subplot(224); plot(linspace(0,1,N), H'*zp, 'green', 'LineWidth', 4); 
title(['Reconstruction with Compressive Sensing (Measurements ', num2str(K), '/', num2str(N), ')'])
%legend('Piecewise linear function x', 'Sampled observations', 'Reconstructed signal')


%% Visualize sparsity level
norm0
N
figure(2); hist(abs(z), 20); title('Histogram of Component Values')

