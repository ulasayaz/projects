clear all
close all
clc

nrTests = 20;
conditions = [1:1:20]; % horizon
T = 50;
nrConditions = length(conditions);
ADD_LIBRARIES

%% settings l1-inf
settings.epsilon = 0.1;
settings.percSamples = 0.05;
settings.optMode = 'l1inf'; 
settings.lambda_reg = -1;
settings.noiseMode = 'linf'; 
settings.doAddNeighbors = 0;

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
countSuccess = 0;
countFailure = 0;
deltaTran = -0.5; % at each time we move 50cm away from the signal
deltaTranEps = 0.1; % we accumulate some odometry noise
for test = 1:nrTests % for each test a different signal
    
    %% generate measurements for all time stamps
    [ZGT, zGT, settings.Nx, settings.Ny, settings.N, rgb] = loadGTsignal_2D(settings);
    interFrameTrueMotion = deltaTran*ones(T-1,1);
    positionsTrue_z = cumsum([0; interFrameTrueMotion]); % integrated true motion
    interFrameOdometry = interFrameTrueMotion + deltaTranEps * (2*rand(T-1,1)-1); % plus measurement noise
    positionsNoisy_z = cumsum([0; interFrameOdometry]); % integrated noisy odometry
    for t=1:T 
       tempDataset(t).ZGT = ZGT - positionsTrue_z(t); % getting deltaTran farther from the image at each step
       tempDataset(t).zGT = zGT - positionsTrue_z(t); % getting deltaTran farther from the image at each step
       tempDataset(t).rgb = rgb;
       tempDataset(t).deltaTranEps = deltaTranEps; % odometry noise
       tempDataset(t).positionsNoisy_z = positionsNoisy_z(t); % integrated odometry from frame 1 to frame t
       tempDataset(t).positionsTrue_z = positionsTrue_z(t); % true motion from frame 1 to frame t

       [tempDataset(t).samples, tempDataset(t).K, tempDataset(t).percSamples] = ...
           createSampleSet_2D(settings, tempDataset(t).ZGT, [], rgb); 
       
       Rfull = speye( settings.N);
       tempDataset(t).R = Rfull(tempDataset(t).samples, :);
       Kt = length(tempDataset(t).samples);
       if settings.addNoise == 1          
           tempDataset(t).noise = settings.epsilon * (2*rand(Kt,1)-1); % entry-wise in [-eps,+eps]
       else
           tempDataset(t).noise = zeros(Kt,1);
       end
       tempDataset(t).y = tempDataset(t).R * tempDataset(t).zGT + tempDataset(t).noise;
    end
    
    %% run tests
    for cond = 1:nrConditions
        horizon = conditions(cond);
        
        for t=T % single test case
            fprintf('(((((((((((((t: %d, cond: %d (%g), test: %d, countSuccess: %d, countFailure: %d, horizon: %d)))))))))))))\n',t,cond,conditions(cond),test,countSuccess,countFailure,horizon)

            pass = false;
            while pass == false
                try
                    tempDataset_tH_t = tempDataset(t-horizon+1:t);
                    results(test,cond) = example_TV2_2D_theory_temporal_function(settings,tempDataset_tH_t);
                    pass = true;
                    countSuccess = countSuccess+1;
                catch ME
                    ME.identifier
                    countFailure = countFailure + 1;
                end
                close all
            end
        end
    end
end
settings.N = settings.d^2;
filename = createFilename(settings,horzcat(matFolder,'results2D_temporal_horizon'));
filename = horzcat(filename,'-DIAG05')
save(filename)
xaxisStr = 'horizon';
CREATE_FIGURES_2D
