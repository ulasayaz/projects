clear all
close all
clc

settings.percSamples = 0.05;
settings.subSample = 0.1; 
settings.lambda = 0; % parameter to weigh the diagonal TV2 norms
settings.theta = 0; % parameter to weigh the TV3 norm
settings.tau = 0.1; % parameter to weigh the TV norm
settings.epsilon = 0.01;
settings.addNoise = 0;
settings.mode = 'video';
settings.nrFrames = 4;
settings.optMode = 'standard'; 
settings.regenerateRandom = 1;
saveVideo = 1;

if saveVideo
    settings.writerObj = VideoWriter(horzcat('./results/','output.avi'));
    settings.writerObj.FrameRate = 2; % 10
    settings.writerObj.Quality = 50;
    open(settings.writerObj);
else
    settings.writerObj = [];
end

timeRange = [1:40]; % [300:400];
for i = 1 : length(timeRange) - settings.nrFrames + 1
    t = timeRange(i);
    settings.msg = sprintf('%.4d',t)
    [errorL1(t+settings.nrFrames-1),errorNaive(t+settings.nrFrames-1)] = example_minTV3Dtemp(settings);
    settings.regenerateRandom = 0;
end
if (saveVideo) close(settings.writerObj); end

%% plot errors
Fig = figure(); hold on;
plot(errorL1,'-b','linewidth',1)
plot(errorNaive,'-r','linewidth',1)
legend('l1','naive')

%% save results to file
filename = horzcat('./results/summary-',settings.mode,'-sub',num2str(100*settings.subSample),'-percSamples',num2str(100*settings.percSamples),'-eps',num2str(log10(settings.epsilon)))
set(Fig, 'paperunits', 'points' )
set(Fig,'position',[0, 0, 300 300]); %
set(Fig,'papersize',[300 300]);
set(Fig,'PaperPositionMode','Auto')
saveas(Fig,filename,'pdf');

%% save mat file
close all
save('./results/results.mat')