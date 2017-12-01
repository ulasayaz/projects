% See problem 4 in the following link for more information regarding
% triangular functions as basis for piece-wise linear functions.
% 
% http://www.math.ubc.ca/~malabika/teaching/ubc/spring08/math421-510/hw4.pdf
%
% Let N = 2^m + 1
% We would like to do compressive sensing on x = R*D*z, where 
%   R is the sampling matrix,
%   D is the triangular functions basis, 
%   z is the sparse representation, 
%   x is our piece-wise linear function with only a few critical points.

clear;
clc;
close all;

addpath('../l1magic/Optimization')
addpath('../grammSchmidt')


m = 10;
N = 2^m + 1;
K = 70 % round(N / 4); % number of measurements
C = 5;  % number of critical points in the piecewise linear function


%% Create piecewise linear function x (true signal)
maxValue = 5;
% x = [0, maxValue * [N-2 : -1 : 0] / (N-1)]';

% define a few quiry points (critical points)
xq = [0, randperm(N-1, C)/N, 1];
yq = (rand(1, C+2) - 0.5) * maxValue;

% generate linear interpolation
u = linspace(0,1,N);
x = interp1(xq,yq,u,'linear')';

figure(1); plot(linspace(0,1,N), x, 'LineWidth', 5); hold on; 

%% Create triangular function basis
% each column of D is a triangular function evaluated at different points
% therefore, each row of D is the set of triangular functions evaluated at
% the same point x, i.e., [f_0(x), f_1(x), ..., f_{N-1}(x)]
D = zeros(N);
for i = 1:N
    D(:, i) = triangle_basis(i-1, N)';
end

%% Direct reconstruction (to check sparsity w.r.t. basis)
z = inv(D) * x;
norm0 = sum(abs(z)>1e-3);
% y = D * z; % figure(2); plot(0:1/(N-1):1, y, 'r'); title('Sparsity')

%% Compressive sensing (recovery of z with l-1 optimization)
% R is the sampling matrix that takes K elements randomly 
R = eye(N);
samples = randperm(N, K);
figure(1); plot( (samples-1)/(N-1), x(samples), '*r', 'MarkerSize', 10)
R = R(samples, :);

% observations
y = R * x;

% initial guess = min energy
z0 = R'*y;

% solve the LP
A = R * D;
tic
zp = l1eq_pd(z0, A, [], y, 1e-3);
toc

% reconstruction
figure(1); plot(linspace(0,1,N), D*zp, 'r', 'LineWidth', 2); 
title(['Compressive Sensing (Measurements ', num2str(K), '/', num2str(N), ')'])
legend('Piecewise linear function x', 'Sampled observations', 'Reconstructed signal')

norm0
N

%% Try to orthonormalize D
Dortho = grammSchmidt(D);
save('Dortho.mat', 'Dortho')
figure; hold on
cmap = hsv(size(Dortho,2));  %# Creates a 6-by-3 set of colors from the HSV colormap
for i=1:10:size(Dortho,2)
    plot(Dortho(:,i),'-s','Color',cmap(i,:));
    pause
end

