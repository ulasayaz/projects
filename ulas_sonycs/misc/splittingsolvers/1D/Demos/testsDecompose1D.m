clear x z y xhat* yhat* C
global dict pars1 pars2 pars3

n=1024;	% Signal length.
maxIters=300;

% Signal.
% Cosine signal
t = ((-n/2):(n/2-1))'/n;
%xdct = cos(pi * (n/16) * t);
xdct=zeros(n,1);
xdct(1:n/4)=cos(2*pi*t(1:n/4)*32);
xdct(n-n/4+1:n)=cos(2*pi*t(1:n/4)*32);
xdct(n/2-n/8+1:n/2+n/8)=cos(2*pi*t(1:n/4)*16);
% Make the gaussian part.
xwav=exp(-1E3.*t.^2)+exp(-2E3.*(t+0.25).^2)+exp(-1.5E3.*(t-0.25).^2);
x = xdct + xwav;

% Representation dictionary (here local DCT + DWT).
qmf=MakeONFilter('Symmlet',6);
%[qmf,dqmf]=MakeBSFilter('CDF',[4 4]);
dict1='RealFourier';pars11=n/8;pars12=0;pars13=0;
dict2='UDWT';pars21=2;pars22=qmf;pars23=0;
dict=MakeList(dict1,dict2);
pars1=MakeList(pars11,pars21);
pars2=MakeList(pars12,pars22);
pars3=MakeList(pars13,pars23);
%E=computeL2norms1D(n,dict,pars1,pars2,pars3);

% Coefficient vector length.
p = SizeOfDict1(n,dict,pars1,pars2,pars3);
tightFrame = p/n; % The dictionary is a tight frame with constant p/n.

% Observed noisy data.
SNR = 15;
sigma = std(x(:))*10^(-SNR/20);
y = x + sigma*randn(size(x));

mu = n/p;			% Relaxation parameter for ProxLasso.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = sqrt(n)*sigma;  	
%epsilon = sigma*sqrt(p)*sqrt(1 + 2*sqrt(2)/sqrt(p)); % Desired residual error. Slightly larger than sqrt(n)*sigma.
OptTol  = sigma/norm(y);
lambdaStop = sigma;	    	% MCA stops when the Lagrange multiplier <= lambdaStop. 
				% Can be used in LassoProx and BPDNDR to remove coeffs with magnitude <= lambdaStop.

% LASSO-Prox may need much more iterations to converge to the same solution as DR-based algorithms. 
% ||x||_1 <= R, R is calibrated from the BPDN-DR solution (here ~ 520).
tic;xhatlassoprox = real(SolveLassoProx('FastDictOp1D', y, p, 520, mu, lambdaStop, maxIters, 0, 0, 0, 0));timelassoprox=toc
tic;xhatBPDR = real(SolveBPDouglasRachford('FastDictOp1D', x, p, gamma, tightFrame, 0, maxIters, 0, 0, 0, 0));timeBPDR=toc
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('FastDictOp1D', y, p, epsilon, gamma, tightFrame, lambdaStop, maxIters, 0, 0, 0, 0));timeBPDNDR=toc
%tic;xhatlars = real(SolveLasso('FastDictOp1D', y, p, 'lars', maxIters, lambdaStop, 0, 0, OptTol));timelars=toc
tic;yhatMCA=SolveMCA(y,dict,pars1,MakeList(1/pars11,pars22),pars3,maxIters,0,3*sigma,epsilon,0);timeMCA=toc

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
	C = eval(['Fast' NAME, 'Analysis(squeeze(yhatMCA(i,:,:)), PAR1, PAR2, PAR3)']);
	xhatMCA((ind+1):(ind+DimDomain)) = C;
	%C = xhatlars((ind+1):(ind+DimDomain));
	%yhatlars(:,i) = eval(['Fast' NAME, 'Synthesis(C, n, PAR1, PAR2, PAR3)']);
	ind = ind + DimDomain;
end


fprintf('%s\n','*'*ones(1,90));
fprintf('%40sSummary\n',' ');
fprintf('%s\n','*'*ones(1,90));
fprintf('%-10s\t%-15s\t%-15s\t%-15s\t%-15s\n',' ','LASSO-Prox','BP-DR','BPDN-DR','MCA');
fprintf('%-10s\t%-15g\t%-15g\t%-15g\t%-15g\n','||x||_0:',length(find(abs(xhatlassoprox))), ...
							      length(find(abs(xhatBPDR))), ...
							      length(find(abs(xhatBPDNDR))), ...
							      length(find(abs(xhatMCA))));

fprintf('%-10s\t%-15g\t%-15g\t%-15g\t%-15g\n','||x||_1:',norm(xhatlassoprox,1), 	  ...
							      norm(xhatBPDR,1), 	  ...
							      norm(xhatBPDNDR,1), 	  ...
							      norm(xhatMCA,1));

fprintf('%-10s\t%-15g\t%-15g\t%-15g\t%-15g\n','CPU (s):', timelassoprox, timeBPDR, timeBPDNDR, timeMCA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot recovered signals.
subplot(5,1,1); 
plot([x y]);axis tight
legend('Original','Noisy');
title(sprintf('Noisy signal SNR_{in}=%g dB',SNR));

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
plot([squeeze(sum(yhatMCA,1)) x]);axis tight;
title(sprintf('MCA Iter=%d PSNR=%g dB',maxIters,psnr(x(:),squeeze(sum(yhatMCA,1)))));

% Plot DCT recovered component.
subplot(5,3,5); 
plot([yhatlassoprox(:,1) xdct]);axis tight
title(sprintf('LASSO-Prox Iter=%d PSNR=%g dB',maxIters,psnr(xdct(:),yhatlassoprox(:,1))));

subplot(5,3,8); 
plot([yhatBPDR(:,1) xdct]);axis tight
title(sprintf('BP-DR Iter=%d PSNR=%g dB',maxIters,psnr(xdct(:),yhatBPDR(:,1))));

subplot(5,3,11); 
plot([yhatBPDNDR(:,1) xdct]);axis tight
title(sprintf('BPDN-DR Iter=%d PSNR=%g dB',maxIters,psnr(xdct(:),yhatBPDNDR(:,1))));

subplot(5,3,14); 
plot([squeeze(yhatMCA(1,:,:)) xdct]);axis tight
title(sprintf('MCA Iter=%d PSNR=%g dB',maxIters,psnr(xdct(:),squeeze(yhatMCA(1,:,:)))));

% Plot Gaussian recovered component.
subplot(5,3,6); 
plot([yhatlassoprox(:,2) xwav]);axis tight
title(sprintf('LASSO-Prox Iter=%d PSNR=%g dB',maxIters,psnr(xwav(:),yhatlassoprox(:,2))));

subplot(5,3,9); 
plot([yhatBPDR(:,2) xwav]);axis tight
title(sprintf('BP-DR Iter=%d PSNR=%g dB',maxIters,psnr(xwav(:),yhatBPDR(:,2))));

subplot(5,3,12); 
plot([yhatBPDNDR(:,2) xwav]);axis tight
title(sprintf('BPDN-DR Iter=%d PSNR=%g dB',maxIters,psnr(xwav(:),yhatBPDNDR(:,2))));

subplot(5,3,15); 
plot([squeeze(yhatMCA(2,:,:)) xwav]);axis tight
title(sprintf('MCA Iter=%d PSNR=%g dB',maxIters,psnr(xwav(:),squeeze(yhatMCA(2,:,:)))));


saveas(gcf,'1D/Datasets/testsDecompose1D.fig','fig');
