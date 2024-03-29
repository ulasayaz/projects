clear all
close all
clc

conditions = [0.02:0.02:1]  % percentage of samples
nrConditions = length(conditions);
ADD_LIBRARIES

%% settings l1
settings.epsilon = 0;
settings.optMode = 'l1inf'; % l1inf is the same for epsilon = 0
settings.lambda_reg = -1;
settings.noiseMode = 'linf'; % irrelevant, noiseless 
settings.doAddNeighbors = 0;

%% common settings
settings.nrCorners = -1; % irrelevant
settings.maxValue = -1; % irrelevant
settings.addNoise = 0; % noiseless
settings.mode = 'pwLinear'; % irrelevant
settings.pwOption = 'diagonal'; % irrelevant
settings.d = 100; % irrelevant
settings.sampleMode = 'uniform';
settings.h = -1;
settings.testDiag = 1;
settings.isDebug = 0;
settings.doAddBoundary = 0;
settings.doAddRandom = 0;
settings.nrNeighbors = 1;
settings.includeAllDirs = 0; % 0 or 1
settings.tol = 1e-5; 
settings.dataset = 'zed';
%settings.path = '~/Dropbox (MIT)/sparseSensing/datasets/';
settings.path = '/Users/ulasayaz/Desktop/Dropbox (MIT)/sparseSensing/datasets/';
settings.imgID = -1;
settings.subSample = 0.2;
settings.useInputImage = 1;
% ===========================
settings.doSaveResults = 0;
settings.thresVect = [];
settings.resultsFolder = matFolder;
 
%% TESTING
end2endTime_tmp = tic;
countSuccess = 0;
countFailure = 0;
countMissingImages = 0;

imgIDS = [400:10:1300];
nrTests = length(imgIDS);
for test = 1:nrTests
    settings.imgID = imgIDS(test);
        for cond = 1:nrConditions
            [cond test]
            settings.percSamples = conditions(cond);
            %% check files exist
            [depthName, rgbName] = getDepthRbgFilenames(settings.path, settings.dataset, settings.imgID);
            if ~(exist(depthName, 'file') == 2) || ~(exist(rgbName, 'file') == 2)
                fprintf('skipping image: %d (non-existent)\n', conditions(cond));
                countMissingImages = countMissingImages+1;
                results(test,cond).e_z = NaN; 
                continue
            end
            %% run actual test
            CORE_noisy_2D
        end
end
settings.N = settings.d^2;
filename = createFilename(settings,horzcat(matFolder,'results2D_zed_percSamples'));
filename = horzcat(filename,'-DIAG05-noEff')
save(filename)
xaxisStr = 'percent. samples';
CREATE_FIGURES_2D
end2endTime = toc(end2endTime_tmp);
