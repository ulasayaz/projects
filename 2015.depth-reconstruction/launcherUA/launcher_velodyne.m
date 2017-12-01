clear all
close all
clc

settings.percSamples = 0.05;
settings.subSample = 1; % only matters on real datasets
settings.lambda = 0; % parameter to weigh the diagonal TV2 norms
settings.theta = 0; % parameter to weigh the TV3 norm
settings.beta = 0; % parameter to weigh the TV norm
settings.epsilon = 0.1;
settings.addNoise = 0;
settings.mode = 'mat' 
settings.optMode = 'standard'; 

timeRange = [250:280];
for t= timeRange
    settings.msg = sprintf('%.4d',t)
    inFilename = horzcat('./dataVelodyne/lids_corridor/lids_corridor (Frame ',settings.msg,').csv'); 
    outFilename = horzcat('./dataVelodyne/lids_corridor_frames/img',settings.msg,'.mat'); 
    readVelodyne(inFilename,outFilename);
    [errorL1(t),errorNaive(t)] = example_minTV2D(settings)
end

Fig = figure(); hold on;
plot(errorL1(timeRange),'-b','linewidth',3)
plot(errorNaive(timeRange),'-r','linewidth',3)
legend('l1','naive')

%% save results to file
filename = horzcat('./results/summary-',mode,'-sub',num2str(100*subSample),'-percSamples',num2str(100*percSamples),'-eps',num2str(log10(epsilon)))
set(Fig, 'paperunits', 'points' )
set(Fig,'position',[0, 0, 300 300]); %
set(Fig,'papersize',[300 300]);
set(Fig,'PaperPositionMode','Auto')
saveas(Fig,filename,'pdf');