clear all
close all
clc

nrTests = 50;
conditions = [1 0.5 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001] % mu
nrConditions = length(conditions);
ADD_LIBRARIES
    
%% settings l1-inf
settings.percSamples = 0.05; 
settings.optMode = 'l1inf'; 
settings.lambda_reg = -1;
settings.noiseMode = 'linf';
settings.doAddNeighbors = 1;

%% common settings
settings.epsilon = 0.1;
settings.nrCorners = 3;
settings.maxValue = 5;
settings.addNoise = 1;
settings.mode = 'pwLinear'; % pwLinear, toy, saddle, stripe
settings.pwOption = 'diagonal'; % pyramid, diagonal
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
settings.d = 100;

%% TESTING
countSuccess = 0;
countFailure = 0;
for test = 1:nrTests
    for cond = 1:nrConditions
        settings.mu = conditions(cond);
        CORE_noisy_2D
    end
end
settings.N = settings.d^2;
filename = createFilename(settings,horzcat(matFolder,'results2D_noisy_nesta_mu'));
filename = horzcat(filename,'-DIAG05')
save(filename)
xaxisStr = 'mu';
CREATE_FIGURES_2D_EXP22
