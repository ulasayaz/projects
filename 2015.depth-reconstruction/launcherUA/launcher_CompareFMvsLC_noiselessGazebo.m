clear all
close all
clc

nrTests = 11;
conditions = [0.1:0.05:1]  % percentage of samples
nrConditions = length(conditions);
addpath('../')
addpath('../../myLib')
addpath('./../lib')
addpath('./../testTheory')
addpath('../../libOPL')
addpath('./../testTheory')

%% settings l1
settings.epsilon = 0;
settings.addNoise = 0; % noiseless
settings.optMode = 'l1inf'; % l1inf is the same for epsilon = 0
settings.lambda_reg = -1;
settings.noiseMode = 'linf'; % irrelevant, noiseless 
settings.doAddNeighbors = 0;
settings.doAddBoundary = 0;
settings.doAddRandom = 0;


%% common settings
settings.nrCorners = -1; % irrelevant
settings.maxValue = -1; % irrelevant

settings.mode = 'pwLinear'; % irrelevant
settings.pwOption = 'diagonal'; % irrelevant
settings.d = 100; % irrelevant
settings.sampleMode = 'uniform';
settings.h = -1;
settings.testDiag = 1;
settings.isDebug = 0;

settings.nrNeighbors = 1;
settings.includeAllDirs = 0; % 0 or 1
settings.tol = 1e-5; 
settings.dataset = 'gazebo';
settings.path = '/Users/ulasayaz/Desktop/Dropbox (MIT)/sparseSensing/datasets/';
settings.imgID = -1;
settings.subSample = 0.2;
settings.useInputImage = 1;
% ===========================
settings.doSaveResults = 0;
settings.thresVect = [];
settings.resultsFolder = '../resultsIROS/';


%% TESTING
countSuccess = 0;
countFailure = 0;

% percSamples
for cond = 1:nrConditions
    
    % images
    for imgID = 2:11
        settings.imgID = imgID;
        
         %% check files exist
        [depthName, rgbName] = getDepthRbgFilenames(settings.path, settings.dataset, settings.imgID);
        if ~(exist(depthName, 'file') == 2) || ~(exist(rgbName, 'file') == 2)
            fprintf('skipping image: %d (non-existent)\n', conditions(cond));
            continue
        end
    
        [imgID cond]
        settings.percSamples = conditions(cond);

        %% run actual test
        fprintf('(((((((((((((percSample: %d (%g), imgID: %d, countSuccess: %d, countFailure: %d)))))))))))))\n',cond,conditions(cond),imgID,countSuccess,countFailure)  
        results(imgID,cond) = example_TV2_2D_theory_function(settings);
        
    end
end

%%

X = conditions;
Y_z = [];
Y_n = [];
Y_zdiag = [];

time_z = [];
time_n = [];
time_zdiag = [];

for cond = 1:nrConditions
    sum_z = 0;
    sum_n = 0;
    sum_zdiag = 0;
    
    sum_time_z = 0;
    sum_time_n = 0;
    sum_time_zdiag = 0;
    
    for i = 2:11
        e_z = results(i, cond).e_z;
        e_n = results(i, cond).e_n;
        e_zdiag = results(i, cond).e_z_diag;
        
        t_z = results(i, cond).cvx_time1;
        t_zdiag = results(i, cond).cvx_time_diag;
        t_n = results(i, cond).naive_time;
        
        sum_z = sum_z + e_z;
        sum_n = sum_n + e_n;
        sum_zdiag = sum_zdiag + e_zdiag;
        
        sum_time_z = sum_time_z + t_z;
        sum_time_n = sum_time_n + t_n;
        sum_time_zdiag = sum_time_zdiag + t_zdiag;
        
    end
    Y_z = [Y_z, sum_z / 10];
    Y_n = [Y_n, sum_n / 10];
    Y_zdiag = [Y_zdiag, sum_zdiag / 10];
    
    time_z = [time_z, sum_time_z / 10]
    time_n = [time_n, sum_time_n / 10];
    time_zdiag = [time_zdiag, sum_time_zdiag / 10];
end

figure;
plot(X, Y_z, '-*r');
hold
plot(X, Y_n, '-*k');
plot(X, Y_zdiag, '-*b')

figure;
plot(X, time_z, '-*r');
hold
plot(X, time_n, '-*k');
plot(X, time_zdiag, '-*b')

