clear all
close all
clc

nrTests = 10;
conditions = [50:50:300] % d
nrConditions = length(conditions);
ADD_LIBRARIES
    
%% settings l1-inf
settings.percSamples = 0.05; 
settings.optMode = 'l1inf'; 
settings.lambda_reg = -1;
settings.noiseMode = 'linf';
settings.doAddNeighbors = 0;

%% common settings
settings.epsilon = 0;
settings.nrCorners = 3;
settings.maxValue = 5;
settings.addNoise = 0;
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

%% TESTING
countSuccess = 0;
countFailure = 0;
for test = 1:nrTests
    for cond = 1:nrConditions
        settings.d = conditions(cond);
        CORE_noisy_2D
    end
end
settings.N = settings.d^2;
filename = createFilename(settings,horzcat(matFolder,'results2D_noiseless_N'));
filename = horzcat(filename,'-DIAG05')
save(filename)
xaxisStr = 'N';
CREATE_FIGURES_2D
