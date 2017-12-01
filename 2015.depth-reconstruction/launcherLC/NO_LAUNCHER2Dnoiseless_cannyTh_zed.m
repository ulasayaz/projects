clear all
close all
clc

nrTests = 1;
conditions =  [0.01:0.01:0.25] % range of canny thresholds
nrConditions = length(conditions);
ADD_LIBRARIES


%% settings l1
settings.percSamples = 1;
settings.epsilon = 0;
settings.optMode = 'l1'; % l1inf is the same for epsilon = 0
settings.lambda_reg = -1;
settings.noiseMode = 'l2'; % irrelevant, noiseless 
settings.doAddNeighbors = 1;

%% common settings
settings.nrCorners = -1; % irrelevant
settings.maxValue = -1; % irrelevant
settings.addNoise = 0; % noiseless
settings.mode = 'pwLinear'; % irrelevant
settings.pwOption = 'diagonal'; % irrelevant
settings.d = 100; % irrelevant
settings.sampleMode = 'edgesRGBrandom';
settings.h = -1;
settings.isDebug = 0;
settings.doAddBoundary = 0;
settings.doAddRandom = 0;
settings.nrNeighbors = 1;
settings.includeAllDirs = 0; % 0 or 1
settings.tol = 1e-5;    
settings.subSample = 0.2;
settings.doSaveResults = 1;
settings.resultsFolder = matFolder;
    
settings.useInputImage = 1;
settings.dataset = 'zed';
settings.imgID = 400;
settings.path = '~/Dropbox/3.FUTURE_WORK/sparseSensing/datasets/';

%% TESTING
countSuccess = 0;
countFailure = 0;
for test = 1:nrTests
    for cond = 1:nrConditions
        [cond test]
        settings.thresVect = conditions(cond);
        %% check files exist
        [depthName, rgbName] = getDepthRbgFilenames(settings.path, settings.dataset, settings.imgID);
        if ~(exist(depthName, 'file') == 2) || ~(exist(rgbName, 'file') == 2)
            fprintf('skipping image: %d (non-existent)\n', conditions(cond));
            continue
        end
        %% run actual test
        CORE_noisy_2D
    end
end
