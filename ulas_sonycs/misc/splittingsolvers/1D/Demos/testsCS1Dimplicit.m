clear x z y xhat* yhat* C
global dict Omega
n=1024; % Number of measurements.
p=n*4; % Coefficient vector length.
k=20;  % Sparsity level.
maxIters=50;

% Stricly or compressible signal.
Compress = 0;
if ~Compress,
	x = SparseVector(p, k, 'GAUSSIAN', true);
else
	%x = besselkforms_rnd(0,1E-3,1,p,1);
	%x = x/max(abs(x));
	s = 0.7;
	x = SparseVector(p, k, 'Signs', true).*(1:p)'.^(-1/s);
	x = x(randperm(p));
	x = x/max(abs(x));
end

% Random measurement (sensing) operator: Hadamard, Fourier, Real Fourier, Real sinusoid (RST), USE, etc.
dict = 'RST';
tightFrame = 1;  % e.g. for Hadamard, Fourier and RST (tight frames), (Phi Phi') = I_n. Otherwise, if unknown or frame, set to 0.
q = randperm(p);
Omega = q(1:n)';

% Observed noisy data.
PSNR = 20;
z = FastMeasure(x, dict, Omega);
sigma = max(abs(z))*10^(-PSNR/20);
y = z + sigma*randn(size(z));

mu = 1; 			% Relaxation parameter for ProxLasso (depends on the sensing matrix redundancy).
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = sqrt(n)*sigma;
%epsilon = sqrt(n)*sigma*sqrt(1 + 2*sqrt(2)/sqrt(n)); % Desired residual error. Slightly larger than sqrt(n)*sigma.
tau = max(eps,sqrt(n/(16*log(n*p)*k)))*sigma*sqrt(log(n*p));
fboptions.beta = 1.99/tightFrame;
OptTol  = sigma/norm(y);
lambdaStop = 3*sigma;	    	% LARS/LASSO stops when the Lagrange multiplier <= lambdaStop. 
lssolution = 1;		    	% If the LS solution is desired, i.e. A_I^+y.


tic;xhatlassoprox = real(SolveLassoProx('FastCSOp', y, p, norm(x,1), mu, lambdaStop, maxIters, lssolution, 0, 0, 0));timelassoprox=toc
tic;xhatBPDR = real(SolveBPDouglasRachford('FastCSOp', z, p, gamma, tightFrame, 0, maxIters, lssolution, 0, 0, 0));timeBPDR=toc
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('FastCSOp', y, p, epsilon, gamma, tightFrame, lambdaStop, maxIters, lssolution, 0, 0, 0));timeBPDNDR=toc
tic;xhatDANTDR = real(SolveDantzigDouglasRachford('FastCSOp', y, p, fboptions, tau, gamma, tightFrame, lambdaStop, maxIters, lssolution, 0, 0, 0));timeDANTDR=toc
tic;xhatlars = real(SolveLasso('FastCSOp', y, p, 'lars', maxIters, lambdaStop, 0, 0, OptTol));timelars=toc
if lssolution, xhatlars = LSsolution('FastCSOp', y, xhatlars, p, find(xhatlars)); end

fprintf('%s\n',char('*'*ones(1,90)));
fprintf('%40sSummary\n',' ');
fprintf('%s\n',char('*'*ones(1,90)));
fprintf('%-10s%-15s\t%-15s\t%-15s\t%-15s\t%-15s\t%-15s\n',' ','Original','LASSO-Prox','BP-DR','BPDN-DR','Dantzig','LARS');
fprintf('%-10s%-15g\t%-15g\t%-15g\t%-15g\t%-15g\t%-15g\n','||x||_0:',length(find(abs(x))),		...
							      length(find(abs(xhatlassoprox))), ...
							      length(find(abs(xhatBPDR))),	...
							      length(find(abs(xhatBPDNDR))),	...
							      length(find(abs(xhatDANTDR))),	...
							      length(find(abs(xhatlars))));   

fprintf('%-10s%-15g\t%-15g\t%-15g\t%-15g\t%-15g\t%-15g\n','||x||_1:',norm(x,1),		...
							      norm(xhatlassoprox,1), 	...
							      norm(xhatBPDR,1), 	...
							      norm(xhatBPDNDR,1), 	...
							      norm(xhatDANTDR,1), 	...
							      norm(xhatlars,1));

fprintf('%-25s\t%-15g\t%-15g\t%-15g\t%-15g\t%-15g\n','CPU (s):', timelassoprox, timeBPDR, timeBPDNDR, timeDANTDR, timelars);

subplot(611)
plot([z y]);
legend('Original','Noisy');
title(sprintf('Projections %s matrix PSNR=%g dB',dict,PSNR));axis tight

subplot(612)
stem(1:p,xhatlassoprox,'r+');hold on
plot(1:p,x,'o');axis([1 p min(x) max(x)]);hold off
title(sprintf('LASSO-Prox MSE=%e',sqrt(norm(xhatlassoprox-x).^2/p)));

subplot(613)
stem(1:p,xhatBPDR,'r+');hold on
plot(1:p,x,'o');axis([1 p min(x) max(x)]);hold off
title(sprintf('BP-DR (noiseless) MSE=%e',sqrt(norm(xhatBPDR-x).^2/p)));

subplot(614)
stem(1:p,xhatBPDNDR,'r+');hold on
plot(1:p,x,'o');axis([1 p min(x) max(x)]);hold off
title(sprintf('BPDN-DR MSE=%e',sqrt(norm(xhatBPDNDR-x).^2/p)));

subplot(615)
stem(1:p,xhatDANTDR,'r+');hold on
plot(1:p,x,'o');axis([1 p min(x) max(x)]);hold off
title(sprintf('Dantzig-DR MSE=%e',sqrt(norm(xhatDANTDR-x).^2/p)));

subplot(616)
stem(1:p,xhatlars,'r+');hold on
plot(1:p,x,'o');axis([1 p min(x) max(x)]);hold off
title(sprintf('LARS MSE=%e',sqrt(norm(xhatlars-x).^2/p)));

saveas(gcf,'1D/Datasets/testsCS1implicit.fig','fig');
