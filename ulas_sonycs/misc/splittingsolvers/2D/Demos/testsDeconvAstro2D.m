clear all
global H h dict pars1 pars2 pars3

x = fitsread('simu_sky.fits');
x = 255*x/max(x(:));
n = length(x);
maxIters=100;

% PSF
h = fftshift(fitsread('simu_psf.fits'));
% Normalize and FFT
h = h/sum(h(:));
H = fft2(h);

% Dictionary (UDWT).
qmf=MakeONFilter('Symmlet',6);
dict='UDWT2';pars1=2;pars2=qmf;pars3=0;
Cx = FastA2(x,dict,pars1,pars2,pars3);

% Coefficient vector length.
p = SizeOfDict2(n,dict,pars1,pars2,pars3);
tightFrame = p/(n*n); % The dictionary is a tight frame with constant p/(n*n).

% Observed noisy data.
BSNR   = 30;
z = real(ifft2(fft2(x).*H));
sigma = std(z(:))*10^(-BSNR/20);
y = z + sigma*randn(size(z));
Cy = FastA2(y,dict,pars1,pars2,pars3);
y = y(:);


mu = 1.9/(max(abs(H(:)))^2*p/(n*n));	% Relaxation parameter for Forward-Backward.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = n*sigma*sqrt(1 + 2*sqrt(2)/n); % Desired residual error. Slightly larger than sqrt(n)*sigma.
lambda = 0.1*sigma; % Regularization parameter.
%lambda = sigma*sqrt(2*log(p))/n; % Theoretical but should be normalized appropriately by the l2-norms of the atoms in the equivalent dictionary H*Phi.
tau = 0.01*max(abs(Cy))*p; % Radius of the l1-norm constraint in the LASS-BPDN l1-constrained form.
			  % Here we use the upper bound: ||x||_1 <= sparsity*||x||_infty = fraction*p*||x||_infty, 0 < fraction <= 1. 

tic;xhatBPDNFB = real(SolveGroupDNForwardBackward('FastConvDictOp2D', y, p, [], lambda, mu, maxIters, 0, 0, 0));timeBPDNFB=toc
tic;xhatlassoFB = real(SolveLassoProx('FastConvDictOp2D', y, p, tau, mu, 0, maxIters, 0, 0, 0, 0));timelassoFB=toc
%tic;xhatBPDNDR = real(SolveBPDouglasRachford('FastConvDictOp2D', y, p, epsilon, gamma, max(abs(H(:)))^2*p/(n*n), 0, maxIters, 0, 0, 0, 0));timeBPDNDR=toc

% Compute recovered images.
yhatBPDNFB = FastS2(xhatBPDNFB,n,dict,pars1,pars2,pars3);
yhatlassoFB = FastS2(xhatlassoFB,n,dict,pars1,pars2,pars3);

fprintf('%s\n','*'*ones(1,90));
fprintf('%40sSummary\n',' ');
fprintf('%s\n','*'*ones(1,90));
fprintf('%-10s%-15s\t%-15s\t%-15s\n',' ','Original','BPDN-FB','LASSO-FB');
fprintf('%-10s%-15g\t%-15g\t%-15g\n','||x||_0:',length(find(abs(Cx))),		...
							      length(find(abs(xhatBPDNFB))), ...
							      length(find(abs(xhatlassoFB))));   

fprintf('%-10s%-15g\t%-15g\t%-15g\n','||x||_1:',norm(Cx,1),		...
							      norm(xhatBPDNFB,1), 	...
							      norm(xhatlassoFB,1));

fprintf('%-25s\t%-15g\t%-15g\t%-15g\t%-15g\n','CPU (s):\n', timeBPDNFB, timelassoFB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display recovered images.
subplot(221); 
imagesc(reshape(x,n,n));axis image;rmaxis
title('Original image');

subplot(222); 
imagesc(reshape(y,n,n));axis image;rmaxis
title(sprintf('Degraded image: blurred and noisy BSNR=%g dB',BSNR));

subplot(223); 
imagesc(reshape(yhatBPDNFB,n,n));axis image;rmaxis
title(sprintf('Deconvolved by solving BPDN in penalized form PSNR=%g dB',psnr(x(:),yhatBPDNFB(:))));

subplot(224); 
imagesc(reshape(yhatlassoFB,n,n));axis image;rmaxis
title(sprintf('Deconvolved by solving BPDN in constrained form PSNR=%g dB',psnr(x(:),yhatlassoFB(:))));

colormap('hsv');

%saveas(gcf,'2D/Datasets/simuskyDeconv2D.fig','fig');
