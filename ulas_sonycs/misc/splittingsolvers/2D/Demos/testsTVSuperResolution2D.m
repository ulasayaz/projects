clear x z y xhat* yhat* tvoptions C
global dict Omega H
% Phantom 
load phantom;
%x = x(1:2:end,1:2:end);
p = length(x);  % Image width.
x = x(:);
maxIters=200;

% PSF
L = 4;
h = [ ones(1,L) , zeros(1, p-L)]'/L;
h = h*h';
h = h/sum(h(:));
H = fft2(h);

% Samples Omega
q = zeros(p,p);
q(1:L:end,1:L:end) = 1;
Omega = find(q(:));
n = floor(sqrt(length(Omega)));

% Measurement operator.
dict = 'Conv';
tightFrame = 0;

% Observed data.
y = FastMeasure2D(x, dict, Omega);

tvoptions.dimension = '2';	% 2 for 2D images. Other fields are set to default values.
tvoptions.numdeep = 8;		% Depth of the dyadic search for the fast TV prox.
tvoptions.lmin = min(x(:));	% The solution is truncated to lmin/lmax.
tvoptions.lmax = max(x(:));
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.

tic;yhatTVDR   = real(SolveTVDouglasRachford('FastCSOp2D', y, p*p, tvoptions, gamma, tightFrame, maxIters, 0, 0, 0));timeTVDR=toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display recovered images.
subplot(211); 
imagesc(reshape(y,n,n));axis image;rmaxis
h = gca;
set(h,'Xlim',[-(p-p/L)/2 (p+p/L)/2],'Ylim',[-(p-p/L)/2 (p+p/L)/2]);axis off
title(sprintf('Degraded image: blurred and 1:%d x 1:%d down-sampled',L,L));

subplot(223); 
imagesc(reshape(x,p,p));axis image;rmaxis
title(sprintf('Original phantom image'));

subplot(224); 
imagesc(real(reshape(yhatTVDR,p,p)));axis image;rmaxis
title(sprintf('TV-DR super-resolved image Iter=%d PSNR=%g dB',maxIters,psnr(x(:),yhatTVDR(:))));

saveas(gcf,'2D/Datasets/phantomSuperResolutionTV2D.fig','fig');

