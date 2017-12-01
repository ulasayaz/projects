clear x z y xhat* yhat* C
global dict pars1 pars2 pars3

n=1024;	% Signal length.
maxIters=300;

% Signal.
% Four cosines have close frequencies, not in the dictionary
fineness = 4;
t = ((1:n)' - .5) / n;freq = pi * ((1:(4*n))' - 1) / (fineness*n);
freq1 = pi * (126.55-1) / (fineness*n);
freq2 = pi * (127.55-1) / (fineness*n);
freq3 = pi * (128.55-1) / (fineness*n);
freq4 = pi * (129.55-1) / (fineness*n);
const = (2/n) ^ .5;
x1 = const * cos(pi * ((126.55  - 1) / fineness) * t);
x2 = const * cos(pi * ((127.55  - 1) / fineness) * t);
x3 = const * cos(pi * ((128.55  - 1) / fineness) * t);
x4 = const * cos(pi * ((129.55  - 1) / fineness) * t);
zerosn = zeros(n,1);
EnergyRatio = 1E1;
xdct = EnergyRatio*(x1 + x2 + x3 + x4);
xdir = SparseVector(n, 3, 'UNIFORM', true);
xdir(find(xdir)) = normsignal(xdir(find(xdir)),1,2);
pos = find(xdir);
x = xdct + xdir;

% Representation dictionary: overcomplete DCT+DIRAC.
dict  = MakeList('DCT','DIRAC');
pars1 = MakeList(fineness,0);
pars2 = MakeList(0,0);
pars3 = MakeList(0,0);

% Coefficient vector length.
p = SizeOfDict1(n,dict,pars1,pars2,pars3);
tightFrame = p/n; % The dictrionary is a tight frame with constant p/n.

% Observed noisy data.
PSNR = 15;
sigma = std(x(:))*10^(-PSNR/20);
y = x + sigma*randn(size(x));

mu = n/p;			% Relaxation parameter for ProxLasso.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = sqrt(n)*sigma;  	
%epsilon = sqrt(n)*sigma*sqrt(1 + 2*sqrt(2)/sqrt(n)); % Desired residual error. Slightly larger than sqrt(n)*sigma.
OptTol  = sigma/norm(y);
lambdaStop = 3*sigma;	    	% LARS/LASSO stops when the Lagrange multiplier <= lambdaStop. 
lssolution = 1;		    	% If the LS solution is desired, i.e. A_I^+y.


% LASSO-Prox may need much more iterations to converge to the same solution as DR-based algorithms. 
% ||x||_1 <= R, R is chosen with an oracle (here 6).
tic;xhatlassoprox = real(SolveLassoProx('FastDictOp1D', y, p, 4*EnergyRatio+norm(xdir,1), mu, lambdaStop/3, maxIters, lssolution, 0, 0, 0));timelassoprox=toc
tic;xhatBPDR = real(SolveBPDouglasRachford('FastDictOp1D', x, p, gamma, tightFrame, 0, maxIters, 0, 0, 0, 0));timeBPDR=toc
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('FastDictOp1D', y, p, epsilon, gamma, tightFrame, lambdaStop/3, maxIters, lssolution, 0, 0, 0));timeBPDNDR=toc
tic;xhatlars = real(SolveLasso('FastDictOp1D', y, p, 'lars', maxIters, lambdaStop, 0, 0, OptTol));timelars=toc
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
fprintf('%-10s%-15s\t%-15s\t%-15s\t%-15s\t%-15s\n',' ','Original','LASSO-Prox','BP-DR','BPDN-DR','LARS');
fprintf('%-10s%-15g\t%-15g\t%-15g\t%-15g\t%-15g\n','||x||_0:',4+length(pos),  ...
							      length(find(abs(xhatlassoprox))), ...
							      length(find(abs(xhatBPDR))), ...
							      length(find(abs(xhatBPDNDR))), ...
							      length(find(abs(xhatlars))));

fprintf('%-10s%-15g\t%-15g\t%-15g\t%-15g\t%-15g\n','||x||_1:',4*EnergyRatio+norm(xdir,1), ...
							      norm(xhatlassoprox,1), 	  ...
							      norm(xhatBPDR,1), 	  ...
							      norm(xhatBPDNDR,1), 	  ...
							      norm(xhatlars,1));

fprintf('%-25s\t%-15g\t%-15g\t%-15g\t%-15g\n','CPU (s):', timelassoprox, timeBPDR, timeBPDNDR, timelars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);clf
% Plot recovered signals.
subplot(5,1,1); 
plot([x y]);axis tight
legend('Original','Noisy');
title(sprintf('Noisy signal PSNR_{in}=%g dB',PSNR));

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

% Plot DCT recovered component.
subplot(5,3,5); 
plot([yhatlassoprox(:,1) xdct]);axis([1 n min(xdct) max(xdct)]);
title(sprintf('LASSO-Prox Iter=%d PSNR=%g dB',maxIters,psnr(xdct(:),yhatlassoprox(:,1))));

subplot(5,3,8); 
plot([yhatBPDR(:,1) xdct]);axis([1 n min(xdct) max(xdct)]);
title(sprintf('BP-DR Iter=%d PSNR=%g dB',maxIters,psnr(xdct(:),yhatBPDR(:,1))));

subplot(5,3,11); 
plot([yhatBPDNDR(:,1) xdct]);axis([1 n min(xdct) max(xdct)]);
title(sprintf('BPDN-DR Iter=%d PSNR=%g dB',maxIters,psnr(xdct(:),yhatBPDNDR(:,1))));

subplot(5,3,14); 
plot([yhatlars(:,1) xdct]);axis([1 n min(xdct) max(xdct)]);
title(sprintf('LARS Iter=%d PSNR=%g dB',maxIters,psnr(xdct(:),yhatlars(:,1))));

% Plot Dirac recovered component.
subplot(5,3,6); 
plot(1:n,yhatlassoprox(:,2),'-b',pos,xdir(pos),'+r');axis([1 n min(xdir) max(xdir)]);
title(sprintf('LASSO-Prox Iter=%d PSNR=%g dB',maxIters,psnr(xdir(:),yhatlassoprox(:,2))));

subplot(5,3,9); 
plot(1:n,yhatBPDR(:,2),'-b',pos,xdir(pos),'+r');axis([1 n min(xdir) max(xdir)]);
title(sprintf('BP-DR Iter=%d PSNR=%g dB',maxIters,psnr(xdir(:),yhatBPDR(:,2))));

subplot(5,3,12); 
plot(1:n,yhatBPDNDR(:,2),'-b',pos,xdir(pos),'+r');axis([1 n min(xdir) max(xdir)]);
title(sprintf('BPDN-DR Iter=%d PSNR=%g dB',maxIters,psnr(xdir(:),yhatBPDNDR(:,2))));

subplot(5,3,15); 
plot(1:n,yhatlars(:,2),'-b',pos,xdir(pos),'+r');axis([1 n min(xdir) max(xdir)]);
title(sprintf('LARS Iter=%d PSNR=%g dB',maxIters,psnr(xdir(:),yhatlars(:,2))));


saveas(gcf,'1D/Datasets/testsDenoise1D_1.fig','fig');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);clf
% Plot recovered coeffs.
subplot(4,3,1); 
stem(xhatlassoprox,'.');axis([1 p min(xhatlassoprox) max(xhatlassoprox)]);
title(sprintf('Coeffs of LASSO-Prox'));

subplot(4,3,4); 
stem(xhatBPDR,'.');axis([1 p min(xhatBPDR) max(xhatBPDR)]);
title(sprintf('Coeffs of BP-DR'));

subplot(4,3,7); 
stem(xhatBPDNDR,'.');axis([1 p min(xhatBPDNDR) max(xhatBPDNDR)]);
title(sprintf('Coeffs of BPDN-DR'));

subplot(4,3,10); 
stem(xhatlars,'.');axis([1 p min(xhatlars) max(xhatlars)]);
title(sprintf('Coeffs of LARS'));

% Plot DCT recovered coeffs.
subplot(4,3,2); 
stem(freq,log(abs(xhatlassoprox(1:n*fineness))+1),'.');
axis1 = axis;
axis([freq(120) freq(140) axis1(3) axis1(4)]);
X1 = [freq1 freq2 freq3 freq4]; X1 = [X1;X1];
Y1 = [axis1(3)*ones(1,4);axis1(4)*ones(1,4)];
hold on
plot(X1, Y1, ':b');
hold off
title(sprintf('DCT coeffs log scale LASSO-Prox'));

subplot(4,3,5); 
stem(freq,log(abs(xhatBPDR(1:n*fineness))+1),'.');
axis1 = axis;
axis([freq(120) freq(140) axis1(3) axis1(4)]);
X1 = [freq1 freq2 freq3 freq4]; X1 = [X1;X1];
Y1 = [axis1(3)*ones(1,4);axis1(4)*ones(1,4)];
hold on
plot(X1, Y1, ':b');
hold off
title(sprintf('DCT coeffs log scale BP-DR'));

subplot(4,3,8); 
stem(freq,log(abs(xhatBPDNDR(1:n*fineness))+1),'.');
axis1 = axis;
axis([freq(120) freq(140) axis1(3) axis1(4)]);
X1 = [freq1 freq2 freq3 freq4]; X1 = [X1;X1];
Y1 = [axis1(3)*ones(1,4);axis1(4)*ones(1,4)];
hold on
plot(X1, Y1, ':b');
hold off
title(sprintf('DCT coeffs log scale BPDN-DR'));

subplot(4,3,11); 
stem(freq,log(abs(xhatlars(1:n*fineness))+1),'.');
axis1 = axis;
axis([freq(120) freq(140) axis1(3) axis1(4)]);
X1 = [freq1 freq2 freq3 freq4]; X1 = [X1;X1];
Y1 = [axis1(3)*ones(1,4);axis1(4)*ones(1,4)];
hold on
plot(X1, Y1, ':b');
hold off
title(sprintf('DCT coeffs log scale LARS'));

% Plot Dirac recovered coeffs.
xmin = min(pos)-30;xmax = max(pos)+30;
subplot(4,3,3); 
stem(log(abs(xhatlassoprox(end-n+1:end))+1),'.');
hold on;
plot(pos,log(abs(xdir(pos))+1),'+r');axis tight
hold off
axis1 = axis;
axis([xmin xmax axis1(3) axis1(4)]);
title(sprintf('Dirac coeffs log scale LASSO-Prox'));

subplot(4,3,6); 
stem(log(abs(xhatBPDR(end-n+1:end))+1),'.');
hold on;
plot(pos,log(abs(xdir(pos))+1),'+r');axis tight
hold off
axis1 = axis;
axis([xmin xmax axis1(3) axis1(4)]);
title(sprintf('Dirac coeffs log scale BP-DR'));

subplot(4,3,9); 
stem(log(abs(xhatBPDNDR(end-n+1:end))+1),'.');
hold on;
plot(pos,log(abs(xdir(pos))+1),'+r');axis tight
hold off
axis1 = axis;
axis([xmin xmax axis1(3) axis1(4)]);
title(sprintf('Dirac coeffs log scale BPDN-DR'));

subplot(4,3,12); 
stem(log(abs(xhatlars(end-n+1:end))+1),'.');
hold on;
plot(pos,log(abs(xdir(pos))+1),'+r');axis tight
hold off
axis1 = axis;
axis([xmin xmax axis1(3) axis1(4)]);
title(sprintf('Dirac coeffs log scale LARS'));

saveas(gcf,'1D/Datasets/testsDenoise1D_2.fig','fig');
