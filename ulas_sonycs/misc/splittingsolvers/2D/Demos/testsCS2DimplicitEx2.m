% needs WaveLab, call WavePath,m to initialize path to include WaveLab 
clear x z y xhat* yhat* C
global dictCS Omega redundancy dict pars1 pars2 pars3
maxIters=100;

% Image 
%load camera;
%x = camera(1:2:end,1:2:end);
x = double(imread('Mondrian.tif'));
x = x(1:2:end,1:2:end);
n = length(x);  % Image size n x n.
x = x(:) - mean(x(:));
maxIters=1000;

% Dictionary (here DWT) (sparsifying trafo)
qmf=MakeONFilter('Haar');
dict='PO2';
pars1=0;pars2=qmf;pars3=0;
p = SizeOfDict2(n,dict,pars1,pars2,pars3);

%Sparsify image.
%cw =FastA2(reshape(x,n,n),dict,pars1,pars2,pars3);
%scw=flipud(sort(abs(cw)));thd=scw(round(n*n*5/100));
%x = FastS2(cw.*(abs(cw)>=thd),n,dict,pars1,pars2,pars3);
%x = x(:);

% Random measurement (sensing) operator.
%m = floor(n/sqrt(2));  		% Number of measurements.
%dictCS = 'RST';
%q = randperm(n*n);
%Omega = [1;q(1:m*m-1)'];


% Measurement (sensing) operator.
dictCS = 'Fourier';
L = 100;  % number of radial lines in the Fourier domain.
[M,Omega] = LineMask(L,n);
m = length(Omega);  	% Number of measurements m.
tightFrame = p/(n*n); 	% The tight frame constant is p/(n*n) x constant of the sensing matrix (e.g. 1 for Fourier or RST).
redundancy = p/(n*n);

% Observed noisy data.
SNR = 30;
z = FastMeasure2D(x, dictCS, Omega);  % apply measurement matrix
sigma = std(z)*10^(-SNR/20);
y = z + sigma*randn(size(z));  % add noise to measurement

mu = 0.99*(n*n)/p; 			% Relaxation parameter for ProxLasso.
gamma = 100;			% Relaxation parameter for Douglas-Rachford iteration.
epsilon = sqrt(m*m)*sigma;   	
%epsilon = sqrt(m*m)*sigma*sqrt(1 + 2*sqrt(2)/sqrt(m*m)); % Desired residual error. Slightly larger than sqrt(m*m)*sigma.
OptTol  = sigma/norm(y);
lambdaStop = 3*sigma;	    	% LARS/LASSO stops when the Lagrange multiplier <= lambdaStop. 
lssolution = 0;		    	% If the LS solution is desired, i.e. A_I^+y.

% A=FastCSDictOp2D
% Lasso: min || y - Ax ||_2^2 s.t. || x ||_1 <= q
tic;xhatlassoprox = real(SolveLassoProx('FastCSDictOp2D', y, p, 6.8E5, mu, lssolution*lambdaStop, maxIters, lssolution, 0, 0, 0));timelassoprox=toc
% BPDR: min || x ||_1 s.t. y = Ax
tic;xhatBPDR = real(SolveBPDouglasRachford('FastCSDictOp2D', z, p, gamma, tightFrame, 0, maxIters, lssolution, 0, 0, 0));timeBPDR=toc
% BPDNDR: min || x ||_1 s.t. || y -Ax ||_2 <= epsilon
tic;xhatBPDNDR = real(SolveBPDNDouglasRachford('FastCSDictOp2D', y, p, epsilon, gamma, tightFrame, lssolution*lambdaStop, maxIters, lssolution, 0, 0, 0));timeBPDNDR=toc

% Compute recovered images by applying dict to optimization result xhat
NAME = NthList(dict, 1);
PAR1 = NthList(pars1, 1);
PAR2 = NthList(pars2, 1);
PAR3 = NthList(pars3, 1);
C = xhatlassoprox;  
yhatlassoprox = eval(['Fast' NAME, 'Synthesis(C(:), n, PAR1, PAR2, PAR3)']);
C = xhatBPDR;
yhatBPDR = eval(['Fast' NAME, 'Synthesis(C(:), n, PAR1, PAR2, PAR3)']);
C = xhatBPDNDR;
yhatBPDNDR = eval(['Fast' NAME, 'Synthesis(C(:), n, PAR1, PAR2, PAR3)']);


fprintf('%s\n','*'*ones(1,90));
fprintf('%40sSummary\n',' ');
fprintf('%s\n','*'*ones(1,90));
fprintf('%-10s\t%-15s\t%-15s\t%-15s\n',' ','LASSO-Prox','BP-DR','BPDN-DR');
fprintf('%-10s\t%-15g\t%-15g\t%-15g\n','||x||_0:',length(find(abs(xhatlassoprox))), ...
						  length(find(abs(xhatBPDR))),	   ...
						  length(find(abs(xhatBPDNDR))));   

fprintf('%-10s\t%-15g\t%-15g\t%-15g\n','||x||_1:',norm(xhatlassoprox,1),    ...
						  norm(xhatBPDR,1),	   ...
						  norm(xhatBPDNDR,1));

fprintf('%-10s\t%-15g\t%-15g\t%-15g\n','CPU (s):',timelassoprox, timeBPDR, timeBPDNDR);

subplot(221)
xft = fftshift(fft2(reshape(x,n,n)));
zft = xft.*fftshift(M);
imagesc(log(1+abs(zft)));axis image;rmaxis
%imagesc(reshape(real(y(1:m*m)),m,m));axis image;rmaxis
title(sprintf('Projections %s m/n=%.2f SNR=%g dB',dictCS,m*m/(n*n),SNR));

subplot(222)
imagesc(reshape(x,n,n));axis image;rmaxis
title(sprintf('Original image'));

subplot(234); 
imagesc(yhatlassoprox);axis image;rmaxis
title(sprintf('LASSO-Prox Iter=%d PSNR=%g dB',maxIters,psnr(x(:),yhatlassoprox(:))));

subplot(235); 
imagesc(yhatBPDR);axis image;rmaxis
title(sprintf('BP-DR (noiseless) Iter=%d PSNR=%g dB',maxIters,psnr(x(:),yhatBPDR(:))));

subplot(236); 
imagesc(yhatBPDNDR);axis image;rmaxis
title(sprintf('BPDN-DR Iter=%d PSNR=%g dB',maxIters,psnr(x(:),yhatBPDNDR(:))));

colormap('gray')

%saveas(gcf,'2D/Datasets/testsCS2DimplicitEx2.fig','fig');
