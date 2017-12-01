clear x z y xhat* yhat* tvoptions C
global A Omega
n=100; % Number of measurements.
p=256; % Signal vector length.
maxIters=300;

% Dictionary matrix: Real Fourier.
F = FourierMat(p);
q = randperm(p/2-1)+1;
Omega = q(1:n/2)';
A = sqrt(2)*[real(F(Omega,:)); imag(F(Omega,:))];
A = [1/sqrt(p)*ones(1,p); A];
x = MakeSignal('Blocks',p)';
x = x - mean(x);
tightFrame = 1;

% Observed noisy data.
SNR = 25;
z = A*x;
sigma = std(z)*10^(-SNR/20);
y = z + sigma*randn(n+1,1);


tvoptions.dimension = '1';	% 1 for 1D vectors. Other fields are set to default values.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = sqrt(n)*sigma; 	
%epsilon = sqrt(n)*sigma*sqrt(1 + 2*sqrt(2)/sqrt(n)); % Desired residual error. Slightly larger than sqrt(n)*sigma.

tic;yhatTVDR   = real(SolveTVDouglasRachford('aprod', z, p, tvoptions, gamma, tightFrame, maxIters, 0, 0, 0));timeTVDR=toc
tic;yhatTVDNDR = real(SolveTVDNDouglasRachford('aprod', y, p, tvoptions, epsilon, gamma, tightFrame, maxIters, 0, 0, 0));timeTVDNDR=toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot recovered signals.
subplot(3,1,1); 
plot(real([z y]));
legend('Noiseless','Noisy');
title(sprintf('Projections Real Fourier SNR=%g dB',SNR));axis tight

subplot(3,1,2); 
plot([yhatTVDR x]);axis tight;
title(sprintf('TV-DR on noiseless data Iter=%d PSNR=%g dB',maxIters,psnr(x(:),yhatTVDR(:))));

subplot(3,1,3); 
plot([yhatTVDNDR x]);axis tight;
title(sprintf('TVDN-DR on noisy data Iter=%d PSNR=%g dB',maxIters,psnr(x(:),yhatTVDNDR(:))));

saveas(gcf,'1D/Datasets/testsCSTV1D.fig','fig');
