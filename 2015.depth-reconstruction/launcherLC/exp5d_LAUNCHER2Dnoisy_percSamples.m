clear all
close all
clc

nrTests = 50;
conditions = [0.02:0.02:1]; % percSamples
nrConditions = length(conditions);
ADD_LIBRARIES

%% settings l1-inf
settings.epsilon = 0.1;
settings.optMode = 'l1inf'; 
settings.lambda_reg = -1;
settings.noiseMode = 'linf'; 
settings.doAddNeighbors = 1;

%% common settings
settings.nrCorners = 3;
settings.maxValue = 5;
settings.addNoise = 1;
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
end2endTime_tmp = tic;
countSuccess = 0;
countFailure = 0;
for test = 1:nrTests
    for cond = 1:nrConditions
        settings.percSamples = conditions(cond);
        CORE_noisy_2D
    end
end
settings.N = settings.d^2;
filename = createFilename(settings,horzcat(matFolder,'results2D_newNoiseModel_percsamples'));
filename = horzcat(filename,'-DIAG05-noEff')
save(filename)
xaxisStr = 'percent. samples';
CREATE_FIGURES_2D
end2endTime = toc(end2endTime_tmp);
