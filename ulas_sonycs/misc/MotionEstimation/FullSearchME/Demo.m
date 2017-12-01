clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo for Subpixel Motion Estimation
% version 1.6
%
% Stanley Chan
% 
% Copyright 2010
% University of California, San Diego
%
% Last modified: 
% 29 Apr, 2010
% 29 Jun, 2010
%  7 Jul, 2010
%  3 Jan, 2013 clean up demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
opts.BlockSize   = 4; % 4 is better than 8 and 16 and same as 2
opts.SearchLimit = 20; % 20 same as 10 but better than 40

% Read image
%img0 = im2double(imread('foreman001.png')); %imgs/foreman001.png
%img1 = im2double(imread('foreman002.png'));%img0(6:end,:);%imrotate(img0,2);% %imgs/foreman002.png
img0 = im2double(imread('tennis00.bmp')); 
img1 = im2double(imread('tennis31.bmp'));
%img0 = img0(1:end-5,:);%img1 = imcrop(img1,[1 1 255 255]);

% Motion estimation
tic
[MVx, MVy] = Bidirectional_ME(img0, img1, opts);
toc

% Motion Compensation
tic
imgMC = reconstruct(img0, MVx, MVy, 0.25); % 0.25 is better than 0.5 
toc

% Evaluation
[M N C] = size(imgMC);
Res  = imgMC-img1(1:M, 1:N, 1:C);
MSE  = norm(Res(:), 'fro')^2/numel(imgMC);
PSNR = 10*log10(max(imgMC(:))^2/MSE);

ResO  = img0(1:M, 1:N, 1:C)-img1(1:M, 1:N, 1:C);
MSEO  = norm(ResO(:), 'fro')^2/numel(imgMC);
PSNRO = 10*log10(max(imgMC(:))^2/MSEO);

% Show results
figure(1);
quiver(MVx(end:-1:1,:), MVy(end:-1:1,:));
title('Motion Vector Field');

figure(2);
subplot(221);
imshow(img0); title('img_0');

subplot(222);
imshow(img1); title('img_1');

subplot(223);
imshow(imgMC); title('img_M');

subplot(224); 
T = sprintf('img_M - img_1, PSNR %3g dB', PSNR);
imshow(rgb2gray(imgMC)-rgb2gray(img1(1:M, 1:N, :)), []); title(T);
