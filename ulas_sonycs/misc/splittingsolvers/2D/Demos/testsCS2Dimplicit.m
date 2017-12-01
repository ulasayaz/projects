clear x z y xhat* yhat* C
global dict Omega
n=64;   % Number of measurements n x n.
p=n*2;  % Coefficient image length.
k=50;   % Sparsity level.
maxIters=100;

% Stricly or compressible image.
Compress = 0;
if ~Compress,
	x = SparseVector(p*p, k, 'GAUSSIAN', true);
else
	s = 0.7;
	x = SparseVector(p*p, k, 'Signs', true).*(1:p*p)'.^(-1/s);
	x = x(randperm(p*p));
	x = x/max(abs(x));
end

% Random measurement (sensing) operator: Hadamard, Fourier, Real Fourier (RST), USE, etc.
dict = 'RST';
tightFrame = 1;  % For Hadamard, RST and Fourier (Parseval tight frames), (Phi Phi') = I_n. Otherwise, if unknown or frame, set to 0.
q = randperm(p*p);
Omega = q(1:n*n)';

% Observed noisy data.
PSNR = 20;
z = FastMeasure2D(x, dict, Omega);
sigma = max(abs(z))*10^(-PSNR/20);
y = z + sigma*randn(size(z));

mu = 1; 			% Relaxation parameter for ProxLasso.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = sqrt(n*n)*sigma;   	
%epsilon = sqrt(n*n)*sigma*sqrt(1 + 2*sqrt(2)/sqrt(n*n)); % Desired residual error. Slightly larger than sqrt(n*n)*sigma.
OptTol  = sigma/norm(y);
lambdaStop = 3*sigma;	    	% LARS/LASSO stops when the Lagrange multiplier <= lambdaStop. 
lssolution = 1;		    	% If the LS solution is desired, i.e. A_I^+y.


tic;xhatlassoprox = real(SolveLassoProx('FastCSOp2D', y, p*p, norm(x,1), mu, lambdaStop, maxIters, lssolution, 0, 0, 0));timelassoprox=toc
tic;xhatBPDR = real(SolveBPDouglasRachford('FastCSOp2D', z, p*p, gamma, tightFrame, 0, maxIters, lssolution, 0, 0, 0));timeBPDR=toc
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('FastCSOp2D', y, p*p, epsilon, gamma, tightFrame, lambdaStop, maxIters, lssolution, 0, 0, 0));timeBPDNDR=toc
tic;xhatlars = real(SolveLasso('FastCSOp2D', y, p*p, 'lars', maxIters, lambdaStop, 0, 0, OptTol));timelars=toc
if lssolution, xhatlars = LSsolution('FastCSOp2D', y, xhatlars, p*p, find(xhatlars)); end

fprintf('%s\n','*'*ones(1,90));
fprintf('%40sSummary\n',' ');
fprintf('%s\n','*'*ones(1,90));
fprintf('%-10s%-15s\t%-15s\t%-15s\t%-15s\t%-15s\n',' ','Original','LASSO-Prox','BP-DR','BPDN-DR','LARS');
fprintf('%-10s%-15g\t%-15g\t%-15g\t%-15g\t%-15g\n','||x||_0:',length(find(abs(x))),		...
							      length(find(abs(xhatlassoprox))), ...
							      length(find(abs(xhatBPDR))),	...
							      length(find(abs(xhatBPDNDR))),	...
							      length(find(abs(xhatlars))));   

fprintf('%-10s%-15g\t%-15g\t%-15g\t%-15g\t%-15g\n','||x||_1:',norm(x,1),		...
							      norm(xhatlassoprox,1), 	...
							      norm(xhatBPDR,1), 	...
							      norm(xhatBPDNDR,1), 	...
							      norm(xhatlars,1));

fprintf('%-25s\t%-15g\t%-15g\t%-15g\t%-15g\n','CPU (s):', timelassoprox, timeBPDR, timeBPDNDR, timelars);

subplot(231)
imagesc(reshape(real(y(1:n*n)),n,n));axis image;rmaxis
title(sprintf('Projections %s matrix PSNR=%g dB',dict,PSNR));

subplot(232)
imagesc(reshape(x,p,p));axis image;rmaxis
title(sprintf('Original image'));

subplot(233)
imagesc(reshape(xhatlassoprox,p,p));axis image;rmaxis
title(sprintf('LASSO-Prox MSE=%e',sqrt(norm(xhatlassoprox-x).^2/(p*p))));

subplot(234)
imagesc(reshape(xhatBPDR,p,p));axis image;rmaxis
title(sprintf('BP-DR (noiseless) MSE=%e',sqrt(norm(xhatBPDR-x).^2/(p*p))));

subplot(235)
imagesc(reshape(xhatBPDNDR,p,p));axis image;rmaxis
title(sprintf('BPDN-DR MSE=%e',sqrt(norm(xhatBPDNDR-x).^2/(p*p))));

subplot(236)
imagesc(reshape(xhatlars,p,p));axis image;rmaxis
title(sprintf('LARS MSE=%e',sqrt(norm(xhatlars-x).^2/(p*p))));

colormap('gray')

saveas(gcf,'2D/Datasets/spikesCS2D.fig','fig');
