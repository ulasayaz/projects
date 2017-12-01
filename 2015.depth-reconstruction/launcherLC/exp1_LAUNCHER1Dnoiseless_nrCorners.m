clear all
close all
clc

nrTests = 50;
conditions = [1:15] % percSamples
nrConditions = length(conditions);
ADD_LIBRARIES

%% settings l1inf
settings.epsilon = 0;
settings.optMode = 'l1inf'; % l1, l1inf, l1reg
settings.lambda_reg = -1;
settings.noiseMode = 'linf' % l2, linf
settings.sampleMode = 'inSegments'; % inSegments, all, onlyCorners, alsoCorners, uniform
    
%% common settings
settings.percSamples = -1;
settings.testL2 = 1;
settings.maxValue = 5;
settings.addNoise = 0;
settings.N = 2000; 
settings.doAddNeighbors = 1; % 1 or 2
settings.doAddBoundary = 1;
settings.isDebug = 0;
settings.tol = 1e-7;

%% TESTING
countSuccess = 0;
countFailure = 0;
for test = 1:nrTests
    for cond = 1:nrConditions
        settings.nrCorners = conditions(cond);
        CORE_noisy_1D
    end
end
filename = createFilename(settings,matFolder,'/results1D_noiseless');
save(filename)
xaxisStr = 'nr Corners';
CREATE_FIGURES_1D 