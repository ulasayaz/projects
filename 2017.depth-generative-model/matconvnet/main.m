close all; clear; clc;

run(fullfile('~', 'matconvnet-1.0-beta25', 'matlab', 'vl_setupnn.m'));

% Define network 
net.layers = {} ;
net.layers{end+1} = struct('type', 'conv', ...
                           'weights', {{0.01*randn(3,3,1,2, 'single'), zeros(1, 2, 'single')}}, ...
                           'stride', 1, ...
                           'pad', 0, ...
                           'dilate', 1) ;
net.layers{end+1} = struct('type', 'relu') ;
net = vl_simplenn_tidy(net)

% Define input data
X = single(randn(5,5,1));

% run the CNN
res = vl_simplenn(net, X) ;
res.x