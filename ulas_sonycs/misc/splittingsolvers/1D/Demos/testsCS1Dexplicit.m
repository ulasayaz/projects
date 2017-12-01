clear x z y xhat* yhat* C
global A
n=256; % Number of measurements.
p=512; % Coefficient vector length.
k=20;  % Sparsity level.
maxIters=2*k;

% Measurement matrix: Dirac + Real Fourier.
F = RSTMat(p);
q = randperm(p);
Omega = q(1:n)';
A = F(Omega,:);
x = SparseVector(p, k, 'GAUSSIAN', true);
tightFrame = 1;

% Observed noisy data.
SNR = 20;
z = A*x;
sigma = std(z)*10^(-SNR/20);
y = z + sigma*randn(n,1);

alpha = max(eig(A*A')); 	% Spectral norm of A*A' (should be 2).
mu = 1/alpha; 			% Relaxation parameter for ProxLasso.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = sqrt(n)*sigma;
%epsilon = sqrt(n)*sigma*sqrt(1 + 2*sqrt(2)/sqrt(n)); % Desired residual error. Slightly larger than sqrt(n)*sigma.
OptTol  = sigma/norm(y);
lambdaStop = 3*sigma;	    	% LARS/LASSO stopes when Lagrange multiplier <= lambdaStop. 
lssolution = 1;		    	% If the LS solution is desired, i.e. A_I^+y.

tic;xhatlassoprox = real(SolveLassoProx('aprod', y, p, norm(x,1), mu, lambdaStop, maxIters, lssolution, 0, 0, 0));timelassoprox=toc
tic;xhatBPDR = real(SolveBPDouglasRachford('aprod', z, p, gamma, tightFrame, 0, maxIters, lssolution, 0, 0, 0));timeBPDR=toc
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('aprod', y, p, epsilon, gamma, tightFrame, lambdaStop, maxIters, lssolution, 0, 0, 0));timeBPDNDR=toc
tic;xhatlars = real(SolveLasso('aprod', y, p, 'lars', maxIters, lambdaStop, 0, 0, OptTol));timelars=toc
if lssolution, xhatlars = LSsolution('aprod', y, xhatlars, p, find(xhatlars)); end

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

subplot(511)
plot([z y]);
legend('Original','Noisy');
title(sprintf('Projections Real Fourier PSNR=%g dB',SNR));axis tight

subplot(512)
stem(1:p,xhatlassoprox,'r+');hold on
plot(1:p,x,'o');axis([1 p min(x) max(x)]);hold off
title(sprintf('LASSO-Prox MSE=%e',sqrt(norm(xhatlassoprox-x).^2/p)));

subplot(513)
stem(1:p,xhatBPDR,'r+');hold on
plot(1:p,x,'o');axis([1 p min(x) max(x)]);hold off
title(sprintf('BP-DR (noiseless) MSE=%e',sqrt(norm(xhatBPDR-x).^2/p)));

subplot(514)
stem(1:p,xhatBPDNDR,'r+');hold on
plot(1:p,x,'o');axis([1 p min(x) max(x)]);hold off
title(sprintf('BPDN-DR MSE=%e',sqrt(norm(xhatBPDNDR-x).^2/p)));

subplot(515)
stem(1:p,xhatlars,'r+');hold on
plot(1:p,x,'o');axis([1 p min(x) max(x)]);hold off
title(sprintf('LARS MSE=%e',sqrt(norm(xhatlars-x).^2/p)));

saveas(gcf,'1D/Datasets/testsCS1explicit.fig','fig');
