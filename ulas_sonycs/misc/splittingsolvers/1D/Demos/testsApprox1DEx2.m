clear x z y xhat* yhat* C
global dict pars1 pars2 pars3

%Signal
x = ReadSignal('Transients'); % Mallat & Zhang (1993).

n = length(x); % Signal length.
maxIters=100;

% Representation dictionary: wavelet packets.
D = log2(n);
qmf = MakeONFilter('Symmlet', 8);
dict  = MakeList('WP');
pars1 = MakeList(D);
pars2 = MakeList(qmf);
pars3 = MakeList(0);

% Coefficient vector length.
p = SizeOfDict1(n,dict,pars1,pars2,pars3);
tightFrame = p/n; % The dictrionary is a tight frame.

% Data vector.
y = x;

mu = n/p;			% Relaxation parameter for ProxLasso.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = 0;   	
OptTol  = 0;
lambdaStop = 0;	    		% LARS/LASSO stops when the Lagrange multiplier <= lambdaStop. 
lssolution = 0;		    	% If the LS solution is desired, i.e. A_I^+y.


% LASSO-Prox may need much more iterations to converge to the same solution as DR-based algorithms. 
% ||x||_1 <= R, R is chosen with an oracle (here ~120, obtained from the BP-DR solution).
tic;xhatlassoprox = real(SolveLassoProx('FastDictOp1D', y, p, 120, mu, lambdaStop, maxIters, lssolution, 0, 0, 0));timelassoprox=toc
tic;xhatBPDR = real(SolveBPDouglasRachford('FastDictOp1D', y, p, gamma, tightFrame, lambdaStop, maxIters, lssolution, 0, 0, 0));timeBPDR=toc
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('FastDictOp1D', y, p, epsilon, gamma, tightFrame, lambdaStop, maxIters, lssolution, 0, 0, 0));timeBPDNDR=toc
tic;xhatlars = real(SolveLasso('FastDictOp1D', y, p, 'lars', n, lambdaStop, 0, 0, OptTol));timelars=toc
if lssolution, xhatlars = LSsolution('FastDictOp1D', y, xhatlars, p, find(xhatlars)); end


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
figure(1);clf					       
subplot(311);
plot(y);axis tight
title('Signal: Transients')

subplot(3,2,3);
PhasePlane(xhatlassoprox, 'WP', n, qmf, 256);
title('Phase Plane LASSO-Prox')

subplot(3,2,4);
PhasePlane(xhatBPDR, 'WP', n, qmf, 256);
title('Phase Plane BP-DR')

subplot(3,2,5);
PhasePlane(xhatBPDNDR, 'WP', n, qmf, 256);
title('Phase Plane BPDN-DR')

subplot(3,2,6);
PhasePlane(xhatlars, 'WP', n, qmf, 256);
title('Phase Plane LARS')

colormap(1-gray)

saveas(gcf,'1D/Datasets/testsApproxTransients.fig','fig');

