clear all
close all
clc

settings.percSamples = 0.05;
settings.subSample = 1; % only matters on real datasets
settings.lambda = 0; % parameter to weigh the diagonal TV2 norms
settings.theta = 0; % parameter to weigh the TV3 norm
settings.beta = 0; % parameter to weigh the TV norm
settings.epsilon = 0.05;
settings.addNoise = 0;
settings.mode = 'pwLinear' 
settings.optMode = 'standard'; 
settings.msg = '0350';
settings.nrCorners = 2;
settings.img = '';
settings.repeat_img = 1;
settings.repeat_samples = 0;

nrTests = 5;

nrImages = 25;

Nx = 100;
Ny = 100;

for i = 1 : nrImages
    i
    close all;
    settings.img = createPWlinear(settings.nrCorners, Nx, Ny);     
    [H,V,H1,V1,H2,V2,H3,V3] = createFiniteDiff2(Nx,Ny);% Create TV matrices
    [VZ_GT, ZH_GT, normGrad_GT, VZ_uncut_GT, ZH_uncut_GT] = computeGrad(settings.img,V,H);

    areas = segmentBinaryImageAndGetAreas(normGrad_GT);
    areas = areas(areas > 3);
    nrAreas = numel(areas);
    minArea(i) = min(areas);
    sparsity(i) = sum(zeroOneTh(normGrad_GT(:)));

    for j = 1 : nrTests
        j
        [errorL1(j,i),errorNaive(j,i)] = example_minTV2D(settings);      
    end
end

save results

%% plot summary
if size(errorL1,1)==1 % vector
    errorL1 = repmat(errorL1,2,1);
    errorNaive = repmat(errorNaive,2,1);
end
meanL1 = mean(errorL1);
meanNaive = mean(errorNaive);
Fig = figure(); hold on;

[sortMeanL1,order] = sort(meanL1);
sortMeanNaive = meanNaive(order);
sortMinArea = minArea(order);
sortSparsity = sparsity(order);

plot(sortMeanL1,'-bo','linewidth',2)
plot(sortSparsity,'-ko','linewidth',3)
plot(sortMinArea,'-ro','linewidth',3)
title(horzcat('lambda =', num2str(lambda), ', theta =', num2str(theta)))
legend(horzcat('l1, ', num2str(ratioL1)),horzcat('l1diag ', num2str(ratioL1diag)),'naive')
grid on
grid minor

%% save results to file
filename = horzcat('./results/summary-',mode,'-sub',num2str(100*subSample),'-percSamples',num2str(100*percSamples),'-eps',num2str(log10(epsilon)))
set(Fig, 'paperunits', 'points' )
%set(Fig,'position',[0, 0, 300 300]); %
set(Fig,'papersize',[300 300]);
set(Fig,'PaperPositionMode','Auto')
saveas(Fig,filename,'pdf');