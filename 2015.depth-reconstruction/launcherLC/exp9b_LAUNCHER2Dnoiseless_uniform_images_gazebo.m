clear all
close all
clc

nrTests = 1;
conditions =  [1:13] % range of images
nrConditions = length(conditions);
ADD_LIBRARIES

%% settings l1
settings.percSamples = 0.02;
settings.epsilon = 0;
settings.optMode = 'l1inf'; % l1inf is the same for epsilon = 0
settings.lambda_reg = -1;
settings.noiseMode = 'linf'; % irrelevant, noiseless 
settings.doAddNeighbors = 1;

%% common settings
settings.nrCorners = -1; % irrelevant
settings.addNoise = 0; % noiseless
settings.mode = 'pwLinear'; % irrelevant
settings.pwOption = 'diagonal'; % irrelevant
settings.d = 100; % irrelevant
settings.sampleMode = 'uniform';
settings.h = -1;
settings.isBounded = 1;
settings.maxValue = 10; % relevant
settings.testDiag = 1;
settings.isDebug = 0;
settings.doAddBoundary = 0;
settings.doAddRandom = 0;
settings.nrNeighbors = 1;
settings.includeAllDirs = 0; % 0 or 1
settings.tol = 1e-5; 
settings.dataset = 'gazebo';
settings.path = '~/Dropbox (MIT)/sparseSensing/datasets/';
% settings.path = '/Users/ulasayaz/Desktop/Dropbox (MIT)/sparseSensing/datasets/';
settings.imgID = -1;
settings.subSample = 0.2;
settings.useInputImage = 1;
% ===========================
settings.doSaveResults = 1;
settings.thresVect = [];
settings.resultsFolder = matFolder;
 
%% TESTING
end2endTime_tmp = tic;
countSuccess = 0;
countFailure = 0;

for test = 1:nrTests
    for cond = 1:nrConditions
        settings.imgID = conditions(cond);
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
settings.N = settings.d^2;
filename = createFilename(settings,horzcat(matFolder,'results2D_gazebo_image_edgesRBG'));
filename = horzcat(filename,'-DIAG05-noEff')
save(filename)
CREATE_FIGURES_2D
