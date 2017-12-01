clear all
close all
clc

nrTests = 10;
subSample = 0.1; % only matters on real datasets
lambda = 0; % parameter to weigh the diagonal TV2 norms
theta = 0; % parameter to weigh the TV3 norm
beta = 0; % parameter to weigh the TV norm
epsilon = 1e-2;
mode = 'load' % load, planes, pwLinear
optMode = 'standard'; % useTVsqrt is very slow, standard, useEdges, useTVsqrt

percSamplesSet = [0.007 0.01 0.05 0.1 0.15 0.2];
for i = 1:length(percSamplesSet)
    percSamples = percSamplesSet(i);
    for j=1:nrTests
        close all; clc; [i j]
        %[errorL1(j,i),errorNaive(j,i)] = example_minTV2D(percSamples, subSample, epsilon, mode, optMode, num2str(j))
        [errorL1(j,i),errorNaive(j,i)] = example_minTV2D_ulas(percSamples, subSample, epsilon, mode, optMode, num2str(j),lambda ,theta, beta)       
    end
end

lambda = 2; % parameter to weigh the diagonal TV2 norms
theta = 0; % parameter to weigh the TV3 norm
beta = 0; % parameter to weigh the TV norm

percSamplesSet = [0.007 0.01 0.05 0.1 0.15 0.2];
for i = 1:length(percSamplesSet)
    percSamples = percSamplesSet(i);
    for j=1:nrTests
        close all; clc; [i j]
        %[errorL1(j,i),errorNaive(j,i)] = example_minTV2D(percSamples, subSample, epsilon, mode, optMode, num2str(j))
        [errorL1diag(j,i),dummy] = example_minTV2D_ulas(percSamples, subSample, epsilon, mode, optMode, num2str(j),lambda ,theta, beta)       
    end
end

save results

%% plot summary
if size(errorL1,1)==1 % vector
    errorL1 = repmat(errorL1,2,1);
    errorL1diag = repmat(errorL1diag,2,1);
    errorNaive = repmat(errorNaive,2,1);
end
ratioL1 = mean(mean(errorL1)./mean(errorNaive));
ratioL1diag = mean(mean(errorL1diag)./mean(errorNaive));
Fig = figure(); hold on;
plot(percSamplesSet,mean(errorL1),'-bo','linewidth',3)
plot(percSamplesSet,mean(errorL1diag),'-ko','linewidth',3)
plot(percSamplesSet,mean(errorNaive),'-ro','linewidth',3)
title(horzcat('lambda =', num2str(lambda), ', theta =', num2str(theta)))
legend(horzcat('l1, ', num2str(ratioL1)),horzcat('l1diag ', num2str(ratioL1diag)),'naive')
grid on
grid minor

%% save results to file
filename = horzcat('./results/summary-',mode,'-sub',num2str(100*subSample),'-percSamples',num2str(100*percSamples),'-eps',num2str(log10(epsilon)))
set(Fig, 'paperunits', 'points' )
set(Fig,'position',[0, 0, 300 300]); %
set(Fig,'papersize',[300 300]);
set(Fig,'PaperPositionMode','Auto')
saveas(Fig,filename,'pdf');