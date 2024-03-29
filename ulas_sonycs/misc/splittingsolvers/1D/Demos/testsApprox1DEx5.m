clear x z y xhat* yhat* C
global dict pars1 pars2 pars3

n = 1024; % Signal length.
maxIters=100;

%Signal
x = InputSignal('WernerSorrows',n); % Chen et al.

% Representation dictionary: cosine packets.
D = log2(n)-4;
dict  = MakeList('CP');
pars1 = MakeList(D);
pars2 = MakeList(0);
pars3 = MakeList(0);

% Coefficient vector length.
p = SizeOfDict1(n,dict,pars1,pars2,pars3);
tightFrame = p/n; % The dictrionary is a tight frame.

% Data vector.
PSNR = 25;
sigma = max(abs(x))*10^(-PSNR/20);
y = x + sigma*randn(size(x));

mu = n/p;			% Relaxation parameter for ProxLasso.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = sqrt(n)*sigma*sqrt(1+2*sqrt(2/n));   	
OptTol  = 0;
lambdaStop = 3*sigma;	    		% LARS/LASSO stops when the Lagrange multiplier <= lambdaStop. 
lssolution = 0;		    	% If the LS solution is desired, i.e. A_I^+y.


% LASSO-Prox may need much more iterations to converge to the same solution as DR-based algorithms. 
% ||x||_1 <= R, R is chosen with an oracle (here ~300 obtained from BP-DR solution).
tic;xhatlassoprox = real(SolveLassoProx('FastDictOp1D', y, p, 300, mu, lambdaStop, maxIters, 0, 0, 0, 0));timelassoprox=toc
tic;xhatBPDR = real(SolveBPDouglasRachford('FastDictOp1D', x, p, gamma, tightFrame, lambdaStop, maxIters, 0, 0, 0, 0));timeBPDR=toc
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('FastDictOp1D', y, p, epsilon, gamma, tightFrame, lambdaStop, maxIters, 0, 0, 0, 0));timeBPDNDR=toc
tic;xhatlars = real(SolveLasso('FastDictOp1D', y, p, 'lars', n, lambdaStop, 0, 0, OptTol));timelars=toc
if lssolution, xhatBPDR   = LSsolution('FastDictOp1D', y, xhatBPDR, p, find(abs(xhatBPDR)>1E-3)); end
if lssolution, xhatBPDNDR = LSsolution('FastDictOp1D', y, xhatBPDNDR, p, find(abs(xhatBPDNDR)>1E-3)); end
if lssolution, xhatlars   = LSsolution('FastDictOp1D', y, xhatlars, p, find(abs(xhatlars))); end


fprintf('%s\n','*'*ones(1,90));
fprintf('%40sSummary\n',' ');
fprintf('%s\n','*'*ones(1,90));
fprintf('%-10s%-15s\t%-15s\t%-15s\t%-15s\n',' ','LASSO-Prox','BP-DR','BPDN-DR','LARS');
fprintf('%-10s%-15g\t%-15g\t%-15g\t%-15g\n','||x||_0:',length(find(abs(xhatlassoprox))), 	...
						       length(find(abs(xhatBPDR))),	     ...
						       length(find(abs(xhatBPDNDR))),	     ...
						       length(find(abs(xhatlars))));

fprintf('%-10s%-15g\t%-15g\t%-15g\t%-15g\n','||x||_1:',norm(xhatlassoprox,1), 	  ...
						       norm(xhatBPDR,1), 	  ...
						       norm(xhatBPDNDR,1), 	  ...
						       norm(xhatlars,1));
						       
fprintf('%-10s%-15g\t%-15g\t%-15g\t%-15g\n','CPU (s)', timelassoprox, timeBPDR, timeBPDNDR, timelars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(1));figure(1);clf					       
subplot(311);
plot(y);axis tight
title('Signal: WernerSorrows')

subplot(3,2,3);
PhasePlane(xhatlassoprox, 'CP', n, 256);
title('Phase Plane LASSO-Prox')

subplot(3,2,4);
PhasePlane(xhatBPDR, 'CP', n, 256);
title('Phase Plane BP-DR')

subplot(3,2,5);
PhasePlane(xhatBPDNDR, 'CP', n, 256);
title('Phase Plane BPDN-DR')

subplot(3,2,6);
PhasePlane(xhatlars, 'CP', n, 256);
title('Phase Plane LARS')

colormap(1-gray)
brighten(-.5)
brighten(-.5)

saveas(gcf,'1D/Datasets/testsApproxWernerSorrows.fig','fig');
