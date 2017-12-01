clear all
close all
clc

conditions = [0.0 :0.02:1] % noise level
nrConditions = length(conditions);
ADD_LIBRARIES

%% settings l1-inf
settings.percSamples = 0.05; 
settings.optMode = 'l1inf'; 
settings.lambda_reg = -1;
settings.noiseMode = 'linf';
settings.doAddNeighbors = 1;

%% common settings
settings.nrCorners = -1; % irrelevant
settings.maxValue = -1; % irrelevant
settings.addNoise = 1; 
settings.mode = 'pwLinear'; % irrelevant
settings.pwOption = 'diagonal'; % irrelevant
settings.d = 100; % irrelevant
settings.sampleMode = 'uniform';
settings.h = -1;
settings.testDiag = 1;
settings.isDebug = 0;
settings.isBounded = 0; % cannot be used here
settings.doAddBoundary = 0;
settings.doAddRandom = 0;
settings.nrNeighbors = 1;
settings.includeAllDirs = 0; % 0 or 1
settings.tol = 1e-5; 
settings.dataset = 'ortho_proj_gazebo';
settings.path = '~/Dropbox (MIT)/sparseSensing/datasets/';
% settings.path = '/Users/ulasayaz/Desktop/Dropbox (MIT)/sparseSensing/datasets/';
settings.imgID = -1;
settings.subSample = 1;
settings.useInputImage = 1;
% ===========================
settings.doSaveResults = 0;
settings.thresVect = [];
settings.resultsFolder = matFolder;
 
%% TESTING
end2endTime_tmp = tic;
countSuccess = 0;
countFailure = 0;

nrTestsPerImage = 5;
imgIDS = kron([1:20], ones(1,nrTestsPerImage));
nrTests = length(imgIDS);
for test = 1:nrTests
    settings.imgID = imgIDS(test);
        for cond = 1:nrConditions
            settings.epsilon = conditions(cond);
            %% check files exist
            [depthName, rgbName] = getDepthRbgFilenames(settings.path, settings.dataset, settings.imgID);
            if ~(exist(depthName, 'file') == 2)
                fprintf('skipping image: %d (non-existent)\n', conditions(cond));
                continue
            end
            %% run actual test
            CORE_noisy_2D
        end
end
settings.N = settings.d^2;
filename = createFilename(settings,horzcat(matFolder,'results2D_gazeboOrtho_epsilon'));
filename = horzcat(filename,'-DIAG05-noEff')
save(filename)
xaxisStr = 'noise level';
CREATE_FIGURES_2D
end2endTime = toc(end2endTime_tmp);
