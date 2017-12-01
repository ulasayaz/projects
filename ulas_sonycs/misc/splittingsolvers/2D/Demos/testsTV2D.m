clear x z y xhat* yhat* tvoptions C
global dict Omega
% Phantom 
load phantom;
%x = x(1:2:end,1:2:end); % Subsample for a faster demo. Comment if undesired.
p = length(x);  % Image width.
x = x(:);
maxIters=200;


% Fourier samples along radial lines
L = 22;  % number of radial lines in the Fourier domain
[M,Omega] = LineMask(L,p);
n = length(Omega);  % Number of measurements n.


% Measurement (sensing) operator.
dict = 'Fourier';
tightFrame = 1;

% Observed noisy data.
sigma = 0.01;
z = FastMeasure2D(x, dict, Omega);
SNR = 20*log10(std(z)/sigma);
y = z + sigma*randn(size(z));


tvoptions.dimension = '2';	% 2 for 2D images. Other fields are set to default values.
tvoptions.numdeep = 8;		% Depth of the dyadic search for the fast TV prox.
tvoptions.lmin = min(x(:));	% The solution is truncated to lmin/lmax.
tvoptions.lmax = max(x(:));
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
%epsilon = sqrt(n)*sigma; 	
epsilon = sqrt(2*n)*sigma*sqrt(1 + 2*sqrt(2)/sqrt(2*n)); % Desired residual error. Slightly larger than sqrt(n*n)*sigma.

tic;yhatTVDR   = real(SolveTVDouglasRachford('FastCSOp2D', z, p*p, tvoptions, gamma, tightFrame, maxIters, 0, 0, 0));timeTVDR=toc
tic;yhatTVDNDR = real(SolveTVDNDouglasRachford('FastCSOp2D', y, p*p, tvoptions, epsilon, gamma, tightFrame, maxIters, 0, 0, 0));timeTVDNDR=toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display recovered images.
subplot(221); 
xft = fftshift(fft2(reshape(x,p,p)));
zft = xft.*fftshift(M);
imagesc(log(1+abs(zft)));axis image;rmaxis
title(sprintf('Projections Real Fourier SNR=%g dB',SNR));

subplot(222); 
imagesc(reshape(x,p,p));axis image;rmaxis
title(sprintf('Original phantom image'));

subplot(223); 
imagesc(real(reshape(yhatTVDR,p,p)));axis image;rmaxis
title(sprintf('TV-DR on noiseless data Iter=%d PSNR=%g dB',maxIters,psnr(x(:),yhatTVDR(:))));

subplot(224); 
imagesc(real(reshape(yhatTVDNDR,p,p)));axis image;rmaxis
title(sprintf('TVDN-DR on noisy data Iter=%d PSNR=%g dB',maxIters,psnr(x(:),yhatTVDNDR(:))));

saveas(gcf,'2D/Datasets/phantomCSTV2D.fig','fig');
