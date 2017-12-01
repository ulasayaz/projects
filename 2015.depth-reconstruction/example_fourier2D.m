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

addpath('../data/traindata')
addpath('../../../spgl1-1.9')

K = 3000; % round(N / 4); % number of measurements


%% Load depth data
img = imread('3.pgm');
img = imresize(img, 0.2); % downsample image
N = prod(size(img));
figure(1); subplot(231); imshow(img); title('Orignal Depth Signal')

%% Create Fourier basis
[h,w] = size(img);
F_h = dftmtx(h);
F_w = dftmtx(w);
F_h_inv = conj(F_h)/h;
F_w_inv = conj(F_w)/w;

%% Sparsify the signal 
disp('2D FFT.')
z = F_h * double(img) * F_w;
% z = (abs(z) > 2e-3 * max(max(abs(z)))) .* z;
x = F_h_inv * z * F_w_inv;
img = uint16(x);
subplot(232); imshow(img); title('Sparsified Signal');
norm0 = sum(sum(abs(z)>0))

%% Take random measurements
% R is the sampling matrix that takes K elements randomly 
disp('Creating samples.')
samples = randperm(N, K);
y = x(samples)';

index = zeros(N, 1);
index(samples) = ones(K, 1);
mask = reshape(index, size(img));

subplot(233); imshow(mask .* x); title('Sampled Observations y');

%% Naive Approach: connect the dots
disp('Naive Approach.')
[Xq, Yq] = meshgrid(1:size(img,2), 1:size(img,1));
X_sample = Xq(samples)';
Y_sample = Yq(samples)';
Fun = scatteredInterpolant(X_sample, Y_sample, y, 'linear');
naive_result = Fun(Xq, Yq);
subplot(234); imshow(uint16(naive_result)); title(['Naive Approach, error=', num2str(norm(x - naive_result,2))])


%% Compressive sensing (recovery of z with l-1 optimization)
%% Note: computationally infeasible!

disp('Constructing matrix A...')
% %% l1-minimization
A = zeros(K, N);
for row = 1 : K
    i = Yq(samples(row));
    j = Xq(samples(row));
    A(row,:) = kron(F_w_inv(:, j).', F_h_inv(i, :));
end
disp('Constructing of A finished.')
% sampled = A * vec(z);
% norm(sampled - y, 2)

disp('l-1 Optimization...')

%% SPGL
epsilon = 0.5;
opts = spgSetParms('verbosity',1);
z_vec = spg_bp(A, y, opts);
% z_vec = spg_bpdn(A, y, epsilon, opts);
% [z_vec,R,G,INFO] = spg_bpdn(A, y, epsilon, opts);

%% cvx 
% cvx_begin
% variable z_vec(N, 1) complex;
% minimize( norm(z_vec,1) )
% subject to
%     y == A * z_vec;
% cvx_end
disp('l-1 Optimization finished.')

subplot(235)
reconstruction = F_h_inv * reshape(z_vec, [h, w]) * F_w_inv;
imshow(uint16(reconstruction));
title(['l1 optimization, error=', num2str(norm(x - reconstruction,2))])
% % reconstruction
% subplot(224); plot(linspace(0,1,N), F_inv * z_vec, 'green', 'LineWidth', 4); 
% title(['Reconstruction with Compressive Sensing (Measurements ', num2str(K), '/', num2str(N), ')'])
% % legend('Piecewise linear function x', 'Sampled observations', 'Reconstructed signal')
% 
% 

% %% Visualize sparsity level
% maxComponent = max(z(:))
% norm0 = sum(abs(z(:))>0.01*max(abs(z(:))))
% N
figure(2); 
% subplot(211); plot(sort(abs(vec(z)))); title('FFT of x')
% subplot(212); plot(sort(abs(z_vec))); title(['Reconstructed sparse signal, error=', num2str(norm(z_vec - vec(z),2))])

subplot(311); plot(abs(vec(z))); title('FFT of x')
subplot(312); plot(abs(z_vec)); title('Reconstructed sparse signal')
subplot(313); plot(abs(z_vec - vec(z))); title(['Error (2-norm=', num2str(norm(z_vec - vec(z),2)), ')'])
norm(y - A*z_vec, 2)
norm(y - A*vec(z), 2)