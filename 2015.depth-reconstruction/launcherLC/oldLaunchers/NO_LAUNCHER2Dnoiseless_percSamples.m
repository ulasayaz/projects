clear all
close all
clc

nrTests = 50;
conditions = [0.02:0.02:1] % percSamples
nrConditions = length(conditions);
addpath('../')
addpath('../../myLib')
addpath('./../lib')
addpath('./../testTheory')
addpath('../../libOPL')
addpath('./../testTheory')

%% settings l1
settings.epsilon = 0;
settings.optMode = 'l1'; 
settings.lambda_reg = 1;
settings.noiseMode = 'l2' 
settings.doAddNeighbors = 0;

%% common settings
settings.nrCorners = 3;
settings.maxValue = 5;
settings.addNoise = 0;
settings.mode = 'pwLinear'; % pwLinear, toy, saddle, stripe
settings.pwOption = 'diagonal'; % pyramid, diagonal
settings.d = 100; % must be odd for testing
settings.sampleMode = 'uniform';
settings.h = -1;
settings.testDiag = 1;
settings.isDebug = 0;
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

%% TESTING
countSuccess = 0;
countFailure = 0;
for test = 1:nrTests
    for cond = 1:nrConditions
        [cond test]
        settings.percSamples = conditions(cond);
        CORE_noisy_2D
    end
end
settings.N = settings.d^2;
filename = createFilename(settings,'../resultsIROS/results2D_percsamples');
save(filename)
xaxisStr = 'percent. samples';
CREATE_FIGURES_2D
