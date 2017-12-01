function [xhat, numIters, activeSet] = SolveDantzigDouglasRachford(A, y, p, fboptions, tau, gamma, tightFrame, lambdaStop, maxIters, lssolution, fullPath, verbose, positivity)
% SolveBPDNDouglasRachford: Iterative Douglas-Rachford Splitting Proximal iteration to solve the BPDN problem:
%		min || x ||_1 s.t. || A'(y -Ax) ||_infty <= tau
%
%	NOTE: The iteration below is rigorously valid for tight frames.
%
%	The solution is given by the iteration :
%		x^(t+1/2) = x^(t) + c^-1 A'(Id - P_C)(y - A x^(t))
%		x^(t+1)   = x^(t) + mu_t (ST_{gamma}[2 x^(t+1/2) - x^(t)] - x^(t+1/2)),
%	where 
%		P_C(z): is the projector onto the closed convex polytope {|| A'z ||_infty <= tau}. This is given by duality and FB
%			splitting such that:
%				u^(i+1) = ST_{beta_i tau}[u^(i) + beta_i A'(z - Au^(i))],  0 < beta_i < 2/||A||^2
%				(Id - P_C)(z) = A ubar, where ubar is the w-limit of u^(i), and A ubar is a s-limit.
%
%		0 < mu_t < 2, fixed to mu_t=1 in this code (optimal convergence speed).
%		ST_{gamma}[x] is soft-thresholding with threshold gamma.
%
%	The solution to the minimization problem is given by x^(t+1/2) at convergence.
%
% Usage
%	[xhat, numIters, activeSet] = SolveDantzigDouglasRachford(A, y, p, fboptions, tau, gamma, tightFrame, lambdaStop, maxIters, lssolution, fullPath, verbose, positivity)
% Input
%	A           nxp matrix or implicit operator
%	y           n-vector
%  	p           length of x
%	fboptions   a structure containing options for FB subiterations to compute P_C:
%			beta		relaxation parameter for the forward-backward implementation case (optional)
%			iters		number of iterations (optional)
%	tau	    maximum residual correlation
%   	gamma	    a strictly positive real, default 1
%	tightFrame  constant of the tight frame if any. Set to zero for frames or tight frames with unknow constant
%   	lambdaStop  If specified (and > 0), the algorithm removes all coeffs <= lambdaStop in magnitude for the LS solution
%	maxIters
%	lssolution  return the least-squares estimate A_I^+y once the support is obtained    
%   	fullPath    1 returns entire solution path, 0 final
%   	verbose     
%	positivity  enforces positivity of the solution.
% Outputs
%	 xhat       solution of BPDN problem
%	 numIters
%	 activeSet  support of the solution 
%
global dict pars1 pars2 pars3

fastop = (ischar(A) || isa(A, 'function_handle'));

n = length(y);

if ~exist('positivity'), positivity = 0; 	end
if ~exist('verbose'), 	 verbose = 0; 		end
if ~exist('fullPath'), 	 fullPath = 0; 		end
if ~exist('lssolution'), lssolution = 0; 	end
if ~exist('maxIters'), 	 maxIters = 50; 	end
if ~exist('lambdaStop'), lambdaStop = 0;	end
if ~exist('tightFrame'), tightFrame = 0;	end
if ~exist('gamma'), 	 gamma = 1;		end
if ~exist('tau'), 	 tau = 0;		end
if ~exist('p'), 	 p = size(A,2); 	end
if ~exist('fboptions'),  
	fboptions.name = 'perform_dantzig_projection';
	if tightFrame, 	fboptions.beta = 1.99/tightFrame;
	else		fboptions.beta = 1.99/(p/n); 	
	end
	fboptions.iters = 200;
	fboptions.fbtol = 1E-3;
end
if ~isfield(fboptions,'name'), 	     fboptions.name = 'perform_dantzig_projection'; end
if ~isfield(fboptions,'beta'),	     
	if tightFrame, 	fboptions.beta = 1.99/tightFrame;
	else		fboptions.beta = 1.99/(p/n); 	  
	end
end
if ~isfield(fboptions,'iters'),      fboptions.iters = 200; end
if ~isfield(fboptions,'fbtol'),      fboptions.fbtol = 1E-3;  end

% lsqrms params
damp = 0;
atol = 1e-6;
btol = 1e-6;
conlim = 1e+10;
itnlim = 10;
show = 0;

xhat = zeros(p,1);
xhat1= zeros(p,1);
activeSet = [];
Ifull = 1:p;
iter = 1;
sols = [];

if ~fastop
 if ~tightFrame
    iA = A'*inv(A*A');
 else
    iA = A'/tightFrame;
 end
end

%while  (iter <= maxIters) & (norm(res) > OptTol*normy)

while  (iter <= maxIters)
    
    if fastop
        res = y - feval(A,1,n,p,xhat,Ifull,p);
	rproj = eval([fboptions.name '(res,A,p,fboptions,tau,xhat1,Ifull)']);
        if ~tightFrame
		corr = lsqrms(n, p, A, 1:p, p, rproj, damp, atol, btol, conlim, itnlim, show);
	else
		corr = feval(A,2,n,p,rproj,Ifull,p)/tightFrame;
	end
    else
    	res = y - A*xhat;
	rproj = eval([fboptions.name '(res,A,p,fboptions,tau,xhat1,Ifull)']);
        corr = iA*rproj;
    end
    xhat1 = xhat + corr;
    xdiff = 2*xhat1 - xhat;
    if isreal(xhat(1)) 
    	xhat  = (abs(xdiff) > gamma) .* (xdiff - gamma.*sign(xdiff)) - corr;  
    else
        xhat  = (abs(xdiff) > gamma) .* (xdiff - gamma.*angle(xdiff)) - corr; % Soft-thresholding also valid for the complex case.
    end
    
    if positivity
  	xhat1(find(xhat1 < 0)) = 0;
    end
  
    activeSet = find(abs(xhat1) > eps);
    
    if (verbose)
        fprintf('Iteration %d: |I| = %d, ||x||_1 = %g\n', iter, length(activeSet), norm(xhat1,1));
	l1norm(iter) = norm(xhat1,1);
	%rss(iter) = norm(res);
	%plot(rss);drawnow
    end

    if fullPath
        sols = [sols xhat1];
    end

    iter = iter+1;
end

if (lambdaStop)
    	xhat1(abs(xhat1) <= lambdaStop) = 0;
end

if (lssolution)
    if lambdaStop
    	xhat = LSsolution(A, y, xhat1, p, find(abs(xhat1) > lambdaStop));
    else
    	xhat = LSsolution(A, y, xhat1, p);
    end
    if fullPath
        sols = [sols xhat];
    end
else
    xhat = xhat1;
end

if (verbose)
        fprintf('Iteration %d: |I| = %d, ||x||_1 = %g\n', iter, length(activeSet), norm(xhat,1));
	l1norm(iter) = norm(xhat,1);
	plot(rss);
end

if fullPath
    xhat = sols;
end

numIters = iter;



function xproj = perform_dantzig_projection(x,A,p,fboptions,tau,uinit,Ifull)
% 	The FB splitting iteration to compute the projection onto the closed convex polytope {|| A'z ||_infty <= tau}. 
%

fastop = (ischar(A) || isa(A, 'function_handle'));
n = length(x);
iter = 1;
uhat = uinit;
uold = realmax;

while  (iter <= fboptions.iters) & (norm(uhat-uold) > fboptions.fbtol)
    
    if fastop
    	xproj = feval(A,1,n,p,uhat,Ifull,p);
        res = x - xproj;
        corr = feval(A,2,n,p,res,Ifull,p)*fboptions.beta;
    else
    	xproj = A*uhat;
    	res = x - xproj;
        corr = A'*res*fboptions.beta;
    end
    uold = uhat;
    uhat = uhat + corr;
    if isreal(uhat(1)) 
    	uhat  = (abs(uhat) > fboptions.beta*tau) .* (uhat - fboptions.beta*tau.*sign(uhat));  
    else
        uhat  = (abs(uhat) > fboptions.beta*tau) .* (uhat - fboptions.beta*tau.*angle(uhat)); % Soft-thresholding also valid for the complex case.
    end
    %if fboptions.verbose
    %   subplot(211)
    %   plot(uhat);axis tight;drawnow
    %   subplot(212);
    %   linfnorm(iter)=norm(feval(A,2,n,p,x-xproj,Ifull,p),inf)-tau;
    %   plot(linfnorm);axis([1 iter+1 -max(abs(linfnorm)) max(abs(linfnorm))]);drawnow
    %end
    iter = iter+1;
end


