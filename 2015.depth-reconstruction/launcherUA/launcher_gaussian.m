clear all
close all
clc

settings.percSamples = 0.01;
settings.subSample = 1; % only matters on real datasets
settings.lambda = 2; % parameter to weigh the diagonal TV2 norms
settings.theta = 0; % parameter to weigh the TV3 norm
settings.beta = 0; % parameter to weigh the TV norm
settings.epsilon = 0.01;
settings.addNoise = 0;
settings.mode = 'pwLinear' % load_(png,jpg,pgm), planes, pwLinear
settings.optMode = 'standard'; 
settings.optMeasure = 'uniform';
settings.prefix = [];
settings.msg = '0350';
settings.nrCorners = 1;
settings.img = '';
settings.repeat_img = 0;
settings.repeat_samples = 0;


[errorUniform,errorNaive] = example_TV2_2D(settings);  

settings.repeat_img = 1;
settings.repeat_samples = 1;
settings.optMeasure = 'Gaussian';

[errorGaussian,errorNaive] = example_TV2_2D(settings); 

