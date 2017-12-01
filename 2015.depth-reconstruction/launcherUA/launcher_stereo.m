clear all
close all
clc

settings.percSamples = 0.05;
settings.subSample = 0.1; 
settings.lambda = 0; % parameter to weigh the diagonal TV2 norms
settings.theta = 0; % parameter to weigh the TV3 norm
settings.beta = 0; % parameter to weigh the TV norm
settings.epsilon = 0.01;
settings.addNoise = 0;
settings.mode = 'load_png';
settings.optMode = 'standard'; 
saveVideo = 1;

if saveVideo
    settings.writerObj = VideoWriter(horzcat('./results/','output.avi'));
    settings.writerObj.FrameRate = 2; % 10
    settings.writerObj.Quality = 50;
    open(settings.writerObj);
    settings.writerObjGrad = VideoWriter(horzcat('./results/','gradients.avi'));
    settings.writerObjGrad.FrameRate = 2; % 10
    settings.writerObjGrad.Quality = 50;
    open(settings.writerObjGrad);
else
    settings.writerObj = [];
    settings.writerObjGrad = [];
end

timeRange = [300:1350];
for t = timeRange
    % close all
%     try
        settings.msg = sprintf('%.4d',t)
        [errorL1(t),errorNaive(t)] = example_minTV2D(settings);
%     catch ME
%         warning('error in example_minTV2D')
%         errorL1(t) = nan;
%         errorNaive(t) = nan;
%     end
end
if (saveVideo) close(settings.writerObj); end
if (saveVideo) close(settings.writerObjGrad); end

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