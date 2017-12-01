% needs path to SparseLab\Utilities added to matlab path

clear x z y xhat* time*
global dict Omega
maxIters=50;
MC=10;
j=[6:14];

for ij=1:length(j)
  n=2^j(ij);  % Number of measurements.
  
  fprintf(['#measurements=', num2str(n),': ']);
  
  p=n*4;      % Coefficient vector length.
  k=ceil(0.05*n);  % Sparsity level.

for rep=1:MC
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

  % Observed data.
  z = FastMeasure(x, dict, Omega);

  gamma = 1;			% Relaxation parameter for Douglas-Rachford iteration.
  OptTol  = 1E-6;
  lssolution = 1;		% If the LS solution is desired, i.e. A_I^+y.


  tic;xhatBPDR = real(SolveBPDouglasRachford('FastCSOp', z, p, gamma, tightFrame, 0, maxIters, lssolution, OptTol, 0, 0));timeBPDR(rep,ij)=toc;
  tic;xhatlars = real(SolveLasso('FastCSOp', z, p, 'lars', maxIters, 0, 0, 0, 0, OptTol));timelars(rep,ij)=toc;
  tic;xhatLP = real(SolveBP('FastCSOp', z, p, maxIters, 0, OptTol));timeLP(rep,ij)=toc;
  tic;xhatStOMP = real(SolveStOMP('FastCSOp', z, p, 'FAR', 0.01, maxIters, 0, OptTol));timeStOMP(rep,ij)=toc;
  errBPDR(rep,ij) = norm(xhatBPDR-x);
  errlars(rep,ij) = norm(xhatlars-x);
  errLP(rep,ij)   = norm(xhatLP-x);
  errStOMP(rep,ij)= norm(xhatStOMP-x);
  
  fprintf('.');
  
end
  disp('.');
  actualloopsize=ij;
  breakopt=5;
  if ij==breakopt
      strResponse = input('continue (y/n) ?', 's');
      if strcmpi(strResponse,'n');          
          break;
      end
  end
end

j=j(1:actualloopsize);

[mean(errBPDR,1)' mean(errlars,1)' mean(errLP,1)' mean(errStOMP,1)']
semilogy(j,mean(timeBPDR,1),'-b',j,mean(timelars,1),'--b',j,mean(timeLP,1),'-.b',j,mean(timeStOMP,1),'.-b');axis tight
set(gca,'FontSize',12);
legend('BP-DR','LARS','LP-Interior Point','StOMP',0);
xlabel('log_2(m)');ylabel('Computation time (s)');
title(sprintf('\\Phi=Dirac, H=%s, m/n=%d%%, sparsity=%.f%%',dict,100*n/p,100*k/n));

%saveas(gcf,sprintf('1D/Datasets/testsCS1%stimes.fig',dict),'fig');
