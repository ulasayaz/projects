% Very rough code for solving compressed sensing with a generative prior.
%
% This code will solve a problem with 2 input dimensions, where the signal 
% is [1,0];
%
% It will output the 2-dimensional latent code retured by the optimizer and
% the residual of that latent code.
%
% This code requires the adam optimizer which can be found at
% https://www.mathworks.com/matlabcentral/fileexchange/61616-adam-stochastic-gradient-descent-optimization
%
% When I ran this code with the provided problem size,
% it would give successful recovery about 75% of the time.
% 
close all; clear; clc;
addpath('DylanMuir-fmin_adam/')

ns = [2, 30, 100];  % Dimensions of input layer, first layer, second layer
m = 25;  % Dimensions of compressed measurements

d = length(ns) - 1;
k = ns(1);
x0 = zeros(k,1); x0(1) = 1;

for ii = 1:d
    W{ii} = randn(ns(ii+1), ns(ii));   % Neural Network Weights
end


xoii = x0;
for ii = 1:d
    xoii = max(0,W{ii}*xoii);        % Compute the relu'd weight matrices
    Wpxo{ii} = W{ii}(xoii >= 0, :);  % corresponding to the true solution x0
end

A = randn(m, ns(end));      % Compressed Sensing matrix
ts = zeros(d,1);

obj = @(x) empirical_risk_obj(x, x0, W, d, A, ts);

sOpt = optimset('fmin_adam');
sOpt.MaxFunEvals = 1e4;
sOpt.Display = 'iter';
sOpt.TolFun = 1e-6;

x_init = randn(k,1);     % Random initializer

tic
x_soln = fmin_adam(obj,x_init,0.01, [], [], [], [], sOpt);  % Solve with Adam optimizer
toc

% Print output
x_soln
resid_latent = norm(x_soln-x0) / norm(x0)

    