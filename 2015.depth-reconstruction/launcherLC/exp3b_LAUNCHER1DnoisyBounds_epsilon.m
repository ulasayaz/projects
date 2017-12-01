clear all
close all
clc

nrTests = 50;
conditions = [0.0:0.02:1] % noise level
nrConditions = length(conditions);
ADD_LIBRARIES
    
%% settings l1inf
settings.percSamples = 0.05;
settings.optMode = 'l1inf'; % l1, l1inf, l1reg
settings.lambda_reg = -1;
settings.noiseMode = 'linf' % l2, linf
settings.sampleMode = 'inSegments'; % inSegments, all, onlyCorners, alsoCorners, uniform

%% common settings
settings.testL2 = 1;
settings.nrCorners = 3; 
settings.maxValue = 5;
settings.addNoise = 1;
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
        settings.epsilon = conditions(cond);
        CORE_noisy_1D
    end
end
filename = createFilename(settings,matFolder,'results1D_s2ineq_bounds_epsilon');
save(filename)
xaxisStr = 'noise level';
CREATE_FIGURES_1D

