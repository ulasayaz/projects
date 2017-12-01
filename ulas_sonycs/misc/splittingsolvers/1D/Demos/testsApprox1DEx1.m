clear x z y xhat* yhat* C
global dict pars1 pars2 pars3

n=512;	% Signal length.
maxIters=300;

% Signal.
x = InputSignal('Star');
x = x(1:n);

% Representation dictionary: overcomplete DCT/DST+DIRAC.
fineness = 2;
dict  = MakeList('DCT','DST','DIRAC');
pars1 = MakeList(fineness,fineness,0);
pars2 = MakeList(0,0,0);
pars3 = MakeList(0,0,0);

% Coefficient vector length.
p = SizeOfDict1(n,dict,pars1,pars2,pars3);
tightFrame = 0; % The dictrionary is not a tight frame because of the DST.

% Data vector.
y = x;

mu = n/p;			% Relaxation parameter for ProxLasso.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = 0;   	
OptTol  = 1E-5;
lambdaStop = 0;	    		% LARS/LASSO stops when the Lagrange multiplier <= lambdaStop. 
lssolution = 0;		    	% If the LS solution is desired, i.e. A_I^+y.


% LASSO-Prox may need much more iterations to converge to the same solution as DR-based algorithms. 
% ||x||_1 <= R, R is chosen with an oracle (here ~860, obtained from the BP-DR solution).
tic;xhatlassoprox = real(SolveLassoProx('FastDictOp1D', y, p, 860, mu, lambdaStop, maxIters, lssolution, 0, 0, 0));timelassoprox=toc
tic;xhatBPDR = real(SolveBPDouglasRachford('FastDictOp1D', y, p, gamma, tightFrame, lambdaStop, maxIters, lssolution, 0, 0, 0));timeBPDR=toc
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('FastDictOp1D', y, p, epsilon, gamma, tightFrame, lambdaStop, maxIters, lssolution, 0, 0, 0));timeBPDNDR=toc
tic;xhatlars = real(SolveLasso('FastDictOp1D', y, p, 'lars', n, lambdaStop, 0, 0, OptTol));timelars=toc
if lssolution, xhatlars = LSsolution('FastDictOp1D', y, xhatlars, p, find(xhatlars)); end


% Compute recovered signals and respective morphological components.
NumberOfDicts = LengthList(dict);
ind = 0;
for i = 1:NumberOfDicts,
	NAME = NthList(dict, i);
	PAR1 = NthList(pars1, i);
	PAR2 = NthList(pars2, i);
	PAR3 = NthList(pars3, i);
	DimDomain = SizeOfDict1(n, NAME, PAR1, PAR2, PAR3);
	C = xhatlassoprox((ind+1):(ind+DimDomain));
	yhatlassoprox(:,i) = eval(['Fast' NAME, 'Synthesis(C, n, PAR1, PAR2, PAR3)']);
	C = xhatBPDR((ind+1):(ind+DimDomain));
	yhatBPDR(:,i) = eval(['Fast' NAME, 'Synthesis(C, n, PAR1, PAR2, PAR3)']);
	C = xhatBPDNDR((ind+1):(ind+DimDomain));
	yhatBPDNDR(:,i) = eval(['Fast' NAME, 'Synthesis(C, n, PAR1, PAR2, PAR3)']);
	C = xhatlars((ind+1):(ind+DimDomain));
	yhatlars(:,i) = eval(['Fast' NAME, 'Synthesis(C, n, PAR1, PAR2, PAR3)']);
	ind = ind + DimDomain;
end

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
% Plot recovered signals.
subplot(5,1,1); 
plot(y);axis tight
title(sprintf('Original signal'));

subplot(5,3,4); 
plot([sum(yhatlassoprox,2) x]);axis tight;
title(sprintf('LASSO-Prox Iter=%d PSNR=%g dB',maxIters,psnr(x(:),sum(yhatlassoprox,2))));

subplot(5,3,7); 
plot([sum(yhatBPDR,2) x]);axis tight;
title(sprintf('BP-DR Iter=%d PSNR=%g dB',maxIters,psnr(x(:),sum(yhatBPDR,2))));

subplot(5,3,10); 
plot([sum(yhatBPDNDR,2) x]);axis tight;
title(sprintf('BPDN-DR Iter=%d PSNR=%g dB',maxIters,psnr(x(:),sum(yhatBPDNDR,2))));

subplot(5,3,13); 
plot([sum(yhatlars,2) x]);axis tight;
title(sprintf('LARS Iter=%d PSNR=%g dB',maxIters,psnr(x(:),sum(yhatlars,2))));

% Plot DCT-DST recovered component.
subplot(5,3,5); 
plot(sum(yhatlassoprox(:,1:2),2));axis tight;
title(sprintf('LASSO-Prox Iter=%d',maxIters));

subplot(5,3,8); 
plot(sum(yhatBPDR(:,1:2),2));axis tight;
title(sprintf('BP-DR Iter=%d',maxIters));

subplot(5,3,11); 
plot(sum(yhatBPDNDR(:,1:2),2));axis tight;
title(sprintf('BPDN-DR Iter=%d',maxIters));

subplot(5,3,14); 
plot(sum(yhatlars(:,1:2),2));axis tight;
title(sprintf('LARS Iter=%d',maxIters));

% Plot Dirac recovered component.
subplot(5,3,6); 
plot(1:n,yhatlassoprox(:,end));axis tight;
title(sprintf('LASSO-Prox Iter=%d',maxIters));

subplot(5,3,9); 
plot(1:n,yhatBPDR(:,end));axis tight;
title(sprintf('BP-DR Iter=%d',maxIters));

subplot(5,3,12); 
plot(1:n,yhatBPDNDR(:,end));axis tight;
title(sprintf('BPDN-DR Iter=%d',maxIters));

subplot(5,3,15); 
plot(1:n,yhatlars(:,end));axis tight;
title(sprintf('LARS Iter=%d',maxIters));

saveas(gcf,'1D/Datasets/testsApproxStar_1.fig','fig');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);clf
% Plot recovered coeffs.
subplot(4,3,1); 
plot(xhatlassoprox);axis([1 p min(xhatlassoprox) max(xhatlassoprox)]);
title(sprintf('Coeffs of LASSO-Prox'));

subplot(4,3,4); 
plot(xhatBPDR);axis([1 p min(xhatBPDR) max(xhatBPDR)]);
title(sprintf('Coeffs of BP-DR'));

subplot(4,3,7); 
plot(xhatBPDNDR);axis([1 p min(xhatBPDNDR) max(xhatBPDNDR)]);
title(sprintf('Coeffs of BPDN-DR'));

subplot(4,3,10); 
plot(xhatlars);axis([1 p min(xhatlars) max(xhatlars)]);
title(sprintf('Coeffs of LARS'));

% Plot DCT-DST recovered coeffs.
subplot(4,3,2); 
plot(xhatlassoprox(1:2*n*fineness));
axis tight;
axis1 = axis;
axis([1 2*n*fineness axis1(3) axis1(4)]);
title(sprintf('DCT-DST coeffs log scale LASSO-Prox'));

subplot(4,3,5); 
plot(xhatBPDR(1:2*n*fineness));
axis tight;
axis1 = axis;
axis([1 2*n*fineness axis1(3) axis1(4)]);
title(sprintf('DCT coeffs log scale BP-DR'));

subplot(4,3,8); 
plot(xhatBPDNDR(1:2*n*fineness));
axis tight;
axis1 = axis;
axis([1 2*n*fineness axis1(3) axis1(4)]);
title(sprintf('DCT coeffs log scale BPDN-DR'));

subplot(4,3,11); 
plot(xhatlars(1:2*n*fineness));
axis tight;
axis1 = axis;
axis([1 2*n*fineness axis1(3) axis1(4)]);
title(sprintf('DCT coeffs log scale LARS'));

% Plot Dirac recovered coeffs.
subplot(4,3,3); 
plot(xhatlassoprox(end-n+1:end));
axis tight;
axis1 = axis;
axis([1 n axis1(3) axis1(4)]);
title(sprintf('Dirac coeffs log scale LASSO-Prox'));

subplot(4,3,6); 
plot(xhatBPDR(end-n+1:end));
axis tight;
axis1 = axis;
axis([1 n axis1(3) axis1(4)]);
title(sprintf('Dirac coeffs log scale BP-DR'));

subplot(4,3,9); 
plot(xhatBPDNDR(end-n+1:end));
axis tight;
axis1 = axis;
axis([1 n axis1(3) axis1(4)]);
title(sprintf('Dirac coeffs log scale BPDN-DR'));

subplot(4,3,12); 
plot(xhatlars(end-n+1:end));
axis tight;
axis1 = axis;
axis([1 n axis1(3) axis1(4)]);
title(sprintf('Dirac coeffs log scale LARS'));

saveas(gcf,'1D/Datasets/testsApproxStar_2.fig','fig');
