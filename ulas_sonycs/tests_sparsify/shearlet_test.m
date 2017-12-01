tic;
fprintf('\n---SLExampleImageDenoising---\n');
fprintf('loading image... ');

clear;
%%settings
sigma = 30;
scales = 4;
thresholdingFactor = 3;

cs=40;

%%load data
fn='cam512.png';
%fn='cameraman.bmp';
%fn='barbara.jpg';
X = imread(fn);
X = double(X);

N=size(X,1);

X=normalize_frame(X);

%%create shearlets
shearletSystem = SLgetShearletSystem2D(0,size(X,1),size(X,2),scales);

elapsedTime = toc;
fprintf([num2str(elapsedTime), ' s\n']);
tic;
fprintf('decomposition, thresholding and reconstruction... ');

%%decomposition
coeffs = SLsheardec2D(X,shearletSystem);

%%thresholding
c_sort = sort(abs(coeffs(:)),'descend'); %sort in descening order

dim=length(c_sort);
thresh = c_sort(floor(dim/cs));

coeffs = coeffs .* (abs(coeffs) > thresh);

%%reconstruction
Xrec = SLshearrec2D(coeffs,shearletSystem);

elapsedTime = toc;
fprintf([num2str(elapsedTime), ' s\n']);

%%compute psnr
%PSNR = SLcomputePSNR(X,Xrec);
PSNR=PSNR3(X,Xrec);
SSIM=SSIM3(X,Xrec);

fprintf(['PSNR: ', num2str(PSNR) , ' db\n']);

figure;
colormap gray;
subplot(1,2,1);
imagesc(X);
title(['original image, res =', num2str(N)]);
axis image;

% subplot(1,3,2);
% imagesc(X);
% title(['noisy image, sigma = ', num2str(sigma)]);
% axis image;

subplot(1,2,2);
imagesc(Xrec);
title(['rec image, cs =',num2str(cs), ' PSNR = ',num2str(PSNR), ...
    ' db, SSIM = ',num2str(SSIM)]);
axis image;

%
%  Copyright (c) 2014. Rafael Reisenhofer
%
%  Part of ShearLab3D v1.01
%  Built Thu, 22/05/2014
%  This is Copyrighted Material