clear all
close all
clc

nrTests = 10;
% settings.percSamples = 0.2; % only matters on real datasets
settings.subSample = 0.2; % only matters on real datasets
settings.lambda = 0; % parameter to weigh the diagonal TV2 norms
settings.theta = 0; % parameter to weigh the TV3 norm
settings.beta = 0; % parameter to weigh the TV norm
settings.epsilon = 1e-2;
settings.addNoise = 0;
settings.mode = 'load_png' 
settings.optMode = 'standard'; 
settings.msg = '';

percSamplesSet = [0.01 0.05 0.1 0.15 0.2 0.25];
for i = 1:length(percSamplesSet)
    settings.percSamples = percSamplesSet(i);
    for j=1:nrTests
        close all; clc; [i j]
        settings.msg = num2str(j);
        [errorL1(j,i),errorNaive(j,i)] = example_minTV2D(settings);
    end
end

%% plot summary
if size(errorL1,1)==1 % vector
    errorL1 = repmat(errorL1,2,1);
    errorNaive = repmat(errorNaive,2,1);
end
Fig = figure(); hold on;
plot(percSamplesSet,mean(errorL1),'-b','linewidth',3)
plot(percSamplesSet,mean(errorNaive),'-r','linewidth',3)
legend('l1','naive')

%% save results to file
filename = horzcat('./results/summary-',mode,'-sub',num2str(100*subSample),'-percSamples',num2str(100*percSamples),'-eps',num2str(log10(epsilon)))
set(Fig, 'paperunits', 'points' )
set(Fig,'position',[0, 0, 300 300]); %
set(Fig,'papersize',[300 300]);
set(Fig,'PaperPositionMode','Auto')
saveas(Fig,filename,'pdf');