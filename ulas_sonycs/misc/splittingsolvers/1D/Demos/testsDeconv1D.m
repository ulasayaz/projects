clear x z y xhat* yhat* C
global H h dict pars1 pars2 pars3

n=1024;	% Signal length.
maxIters=300;

% Signal.
k = 20;
x = SparseVector(n, k, 'GAUSSIAN', true);
I = find(x);

% Coefficient vector length.
p = n;

% PSF
  % Gaussian
	%h = fftshift(normal_pdf(linspace(-2,2,n),0,0.0015))';
	%h = real(ifft(fft(h)+0.1));
	%h = fftshift(exp(-[-n/2:n/2-1].^2/(2*5^2)))';
  % Dirac
  	%h = zeros(n,1);h(1)=1;
  % Moving-Average 
	%q = 7;
	%h = shift([ ones(1,q) , zeros(1, n-q)]'/q,-floor(q/2));
  % Exponential
  	nu = 0.9;
	h = shift(exp(-abs([-n/2:n/2-1]').^nu/2),-n/2);
  % Algebraic
	%nu=1;
	%h = shift(1./(1+abs([-n/2:n/2-1]).^nu),-n/2)';
	%h = shift(1./(1+abs(([-n/2:n/2-1]/n)./0.005)).^4,-n/2)';
% Normalize and FFT
h = h/sum(h(:));
H = fft(h);
  
% Observed noisy data.
PSNR = 25;
z = real(ifft(fft(x).*H));
sigma = std(z(:))*10^(-PSNR/20);
y = z + sigma*randn(size(z));

mu = 1/(max(abs(H))^2*p/n);	% Relaxation parameter for ProxLasso.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
%epsilon = sqrt(n)*sigma;   	
epsilon = sqrt(n)*sigma*sqrt(1 + 2*sqrt(2)/sqrt(n)); % Desired residual error. Slightly larger than sqrt(n)*sigma.
OptTol  = sigma/norm(y);
lambdaStop = 3*sigma;	    	% LARS/LASSO stops when the Lagrange multiplier <= lambdaStop. 
lssolution = 1;		    	% If the LS solution is desired, i.e. A_I^+y.


% LASSO-Prox needs ~10 times more iterations to converge to the same solution. Here ||x||_1 is chosen with an oracle.
tic;xhatlassoprox = real(SolveLassoProx('FastConvOp1D', y, p, norm(x,1), mu, 0, maxIters, 0, 0, 0, 0));timelassoprox=toc
tic;xhatBPDR = real(SolveBPDouglasRachford('FastConvOp1D', y, p, gamma, 0, lambdaStop, maxIters, lssolution, 0, 0, 0));timeBPDR=toc
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('FastConvOp1D', y, p, epsilon, gamma, 0, lambdaStop, maxIters, lssolution, 0, 0, 0));timeBPDNDR=toc
tic;xhatlars = real(SolveLasso('FastConvOp1D', y, p, 'lars', maxIters, lambdaStop, 0, 0, OptTol));timelars=toc
if lssolution, xhatlars = LSsolution('FastConvOp1D', y, xhatlars, p, find(xhatlars)); end

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

t = [1:p];
subplot(511); 
plot(t,z,'k',t,y,'--b');axis tight
legend('Original','Noisy and blurred');
title(sprintf('Noisy blurred spikes PSNR_{in}=%g dB',PSNR));

subplot(512);
stem(t(xhatlassoprox~=0),xhatlassoprox(xhatlassoprox~=0),'.b');hold on;stem(t(I),x(I),'or');axis([1 p min(x) max(x)]);
title(sprintf('LASSO-Prox Iter=%d PSNR=%g dB',maxIters,psnr(x(:),xhatlassoprox(:))));

subplot(513); 
stem(t(xhatBPDR~=0),xhatBPDR(xhatBPDR~=0),'.b');hold on;stem(t(I),x(I),'or');axis([1 p min(x) max(x)]);
title(sprintf('BP-DR Iter=%d PSNR=%g dB',maxIters,psnr(x(:),xhatBPDR(:))));

subplot(514);
stem(t(xhatBPDNDR~=0),xhatBPDNDR(xhatBPDNDR~=0),'.b');hold on;stem(t(I),x(I),'or');axis([1 p min(x) max(x)]);
title(sprintf('BPDN-DR Iter=%d PSNR=%g dB',maxIters,psnr(x(:),xhatBPDNDR(:))));

subplot(515);
stem(t(xhatlars~=0),xhatlars(xhatlars~=0),'.b');hold on;stem(t(I),x(I),'or');axis([1 p min(x) max(x)]);
title(sprintf('LARS Iter=%d PSNR=%g dB',maxIters,psnr(x(:),xhatlars(:))));

saveas(gcf,'1D/Datasets/testsDeconv1D.fig','fig');
