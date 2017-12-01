clear x z y xhat* yhat* C
global dict pars1 pars2 pars3

%%%%%%%%%%%%%%% 2D synthetic image: texture+gaussians %%%%%%%%%%%%%%%%%%%%%%%%
x = double(imread('barbara_512x512.png'));
x = x(end/2-end/4+1:end/2+end/4,end/2-end/4+1:end/2+end/4);
%x = x(1:2:end,1:2:end);
n = length(x);
x = x(:);
maxIters=200;

% Dictionary (here local DCT + DWT).
qmf=MakeONFilter('Symmlet',6);
dict1='UDWT2';pars11=2;pars12=qmf;pars13=0;
dict2='LDCT2iv';pars21=0;pars22=32;pars23=0;
dict=MakeList(dict1,dict2);
pars1=MakeList(pars11,pars21);
pars2=MakeList(pars12,pars22);
pars3=MakeList(pars13,pars23);
%E=computeL2norms2D(n,dict,pars1,pars2,pars3);

% Coefficient vector length.
p = SizeOfDict2(n,dict,pars1,pars2,pars3);
tightFrame = p/(n*n); % The dictionary is a tight frame with constant p/(n*n).


% Observed noisy data.
SNR = 15;
sigma = std(x(:))*10^(-SNR/20);
y = x + sigma*randn(size(x));

mu = (n*n)/p;			% Relaxation parameter for ProxLasso.
gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = sqrt(n*n)*sigma;  	
%epsilon = sqrt(n*n)*sigma*sqrt(1 + 2*sqrt(2)/sqrt(n*n)); % Desired residual error. Slightly larger than sqrt(n*n)*sigma.
OptTol  = sigma/norm(y);
lambdaStop = 3*sigma;	    	% MCA stops when the Lagrange multiplier <= lambdaStop. 
				% Can be used in LassoProx and BPDNDR to remove coeffs with magnitude <= lambdaStop.

% LASSO-Prox may need much more iterations to converge to the same solution as DR-based algorithms. 
% ||x||_1 <= R, R is chosen with an oracle (here 600).
tic;xhatlassoprox = real(SolveLassoProx('FastDictOp2D', y, p, 3.42e+07, mu, 0, maxIters, 0, 0, 0, 0));timelassoprox=toc
tic;xhatBPDR = real(SolveBPDouglasRachford('FastDictOp2D', x, p, gamma, tightFrame, 0, maxIters, 0, 0, 0, 0));timeBPDR=toc
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('FastDictOp2D', y, p, epsilon, gamma, tightFrame, 0, maxIters, 0, 0, 0, 0));timeBPDNDR=toc
tic;yhatMCA=SolveMCA(reshape(y,n,n),dict,pars1,pars2,pars3,maxIters,0,lambdaStop,0,0);timeMCA=toc

% Compute recovered signals and respective morphological components.
NumberOfDicts = LengthList(dict);
ind = 0;
for i = 1:NumberOfDicts,
	NAME = NthList(dict, i);
	PAR1 = NthList(pars1, i);
	PAR2 = NthList(pars2, i);
	PAR3 = NthList(pars3, i);
	DimDomain = SizeOfDict2(n, NAME, PAR1, PAR2, PAR3);
	C = xhatlassoprox((ind+1):(ind+DimDomain));
	yhatlassoprox(i,:,:) = eval(['Fast' NAME, 'Synthesis(C, n, PAR1, PAR2, PAR3)']);
	C = xhatBPDR((ind+1):(ind+DimDomain));
	yhatBPDR(i,:,:) = eval(['Fast' NAME, 'Synthesis(C, n, PAR1, PAR2, PAR3)']);
	C = xhatBPDNDR((ind+1):(ind+DimDomain));
	yhatBPDNDR(i,:,:) = eval(['Fast' NAME, 'Synthesis(C, n, PAR1, PAR2, PAR3)']);
	C = eval(['Fast' NAME, 'Analysis(squeeze(yhatMCA(i,:,:)), PAR1, PAR2, PAR3)']);
	xhatMCA((ind+1):(ind+DimDomain)) = C;
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
% Display recovered images.
subplot(5,1,1); 
imagesc(reshape(y,n,n));axis image;rmaxis
title(sprintf('Noisy image SNR_{in}=%g dB',SNR));

subplot(5,3,4); 
imagesc(squeeze(sum(yhatlassoprox,1)));axis image;rmaxis
title(sprintf('LASSO-Prox Iter=%d PSNR=%g dB',maxIters,psnr(x(:),sum(yhatlassoprox,1))));

subplot(5,3,7); 
imagesc(squeeze(sum(yhatBPDR,1)));axis image;rmaxis
title(sprintf('BP-DR (noiseless) Iter=%d PSNR=%g dB',maxIters,psnr(x(:),sum(yhatBPDR,1))));

subplot(5,3,10); 
imagesc(squeeze(sum(yhatBPDNDR,1)));axis image;rmaxis
title(sprintf('BPDN-DR Iter=%d PSNR=%g dB',maxIters,psnr(x(:),sum(yhatBPDNDR,1))));

subplot(5,3,13); 
imagesc(squeeze(sum(yhatBPDNDR,1)));axis image;rmaxis
title(sprintf('MCA Iter=%d PSNR=%g dB',maxIters,psnr(x(:),sum(yhatMCA,1))));

% Display recovered cartoon component.
subplot(5,3,5); 
imagesc(squeeze(yhatlassoprox(1,:,:)));axis image;rmaxis
title(sprintf('LASSO-Prox Cartoon Iter=%d',maxIters));

subplot(5,3,8); 
imagesc(squeeze(yhatBPDR(1,:,:)));axis image;rmaxis
title(sprintf('BP-DR (noiseless) Cartoon Iter=%d',maxIters));

subplot(5,3,11); 
imagesc(squeeze(yhatBPDNDR(1,:,:)));axis image;rmaxis
title(sprintf('BPDN-DR Cartoon Iter=%d',maxIters));

subplot(5,3,14); 
imagesc(squeeze(yhatMCA(1,:,:)));axis image;rmaxis
title(sprintf('MCA Cartoon Iter=%d',maxIters));

% Display texture recovered component.
subplot(5,3,6); 
imagesc(squeeze(yhatlassoprox(2,:,:)));axis image;rmaxis
title(sprintf('LASSO-Prox Texture Iter=%d',maxIters));

subplot(5,3,9); 
imagesc(squeeze(yhatBPDR(2,:,:)));axis image;rmaxis
title(sprintf('BP-DR (noiseless) Texture Iter=%d',maxIters));

subplot(5,3,12); 
imagesc(squeeze(yhatBPDNDR(2,:,:)));axis image;rmaxis
title(sprintf('BPDN-DR Texture Iter=%d',maxIters));

subplot(5,3,15); 
imagesc(squeeze(yhatMCA(2,:,:)));axis image;rmaxis
title(sprintf('MCA Texture Iter=%d',maxIters));

colormap('gray')
saveas(gcf,'2D/Datasets/BarbaraDecompose.fig','fig');
save barbaradecompose.mat
