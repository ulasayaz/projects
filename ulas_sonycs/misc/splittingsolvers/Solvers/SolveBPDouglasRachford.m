function [xhat, numIters, activeSet] = SolveBPDouglasRachford(A, y, p, gamma, tightFrame, lambdaStop, maxIters, lssolution, OptTol, verbose, positivity)
% SolveBPDouglasRachford: Iterative Douglas-Rachford Splitting Proximal iteration to solve the BP problem:
%		min || x ||_1 s.t. y = Ax
%	The solution is given by the iteration :
%		x^(t+1/2) = P_C(x^(t)) = x^(t) + pinv(A) (y - Ax^(t)) (P_C: projector onto the equality constraint set), 
%		x^(t+1)   = x^(t) + mu_t (ST[2 P_C(x^(t)) - x^(t)] - x^(t+1/2)),
%	where 
%		0 < mu_t < 2, fixed to mu_t=1 in this code (optimal convergence speed).
%		ST[x] is soft-thresholding with threshold gamma for gamma strictly positive.
%
%	The solution to the minimization problem is given by P_C(x^(t)) at convergence.
%
% Usage
%	[xhat, numIters, activeSet] = SolveBPDouglasRachford(A, y, p, gamma, tightFrame, lambdaStop, maxIters, lssolution, OptTol, verbose, positivity)
% Input
%	A           nxp matrix or implicit operator
%	y           n-vector
%  	p           length of x
%   	gamma	    a strictly positive real, default 1
%	tightFrame  constant of the tight frame if any. Set to zero for frames or tight frames with unknow constant
%   	lambdaStop  If specified (and > 0), the algorithm removes all coeffs <= lambdaStop in magnitude for the LS solution
%	maxIters
%	lssolution  return the least-squares estimate A_I^+y once the support is obtained    
%   	OptTol      algorithm terminates if the objective change (l1 norm) is less that OptTol
%   	verbose     
%	positivity  enforces positivity of the solution.
% Outputs
%	 xhat       solution of BP problem
%	 numIters
%	 activeSet  support of the solution 
%
global dict pars1 pars2 pars3

fastop = (ischar(A) || isa(A, 'function_handle'));

n = length(y);

if ~exist('positivity'), positivity = 0; 	end
if ~exist('verbose'), 	 verbose = 0; 		end
if ~exist('OptTol'), 	 OptTol = 0; 		end
if ~exist('lssolution'), lssolution = 0; 	end
if ~exist('maxIters'), 	 maxIters = 50; 	end
if ~exist('lambdaStop'), lambdaStop = 0;	end
if ~exist('tightFrame'), tightFrame = 0;	end
if ~exist('gamma'), 	 gamma = 1;		end
if ~exist('p'), 	 p = size(A,2); 	end

% lsqrms params
damp = 0;
atol = 1e-6;
btol = 1e-6;
conlim = 1e+10;
itnlim = 10;
show = 0;
    
xhat = zeros(p,1);
activeSet = [];
Ifull = 1:p;
iter = 1;
sols = [];
xhat1new = realmax;
xhat1 = 0;
res = y;

if ~fastop
 if ~tightFrame
    iA = A'*inv(A*A');
 else
    iA = A'/tightFrame;
 end
end

while  (iter <= maxIters) & (norm(res) > OptTol*norm(y)) %(norm(xhat1new - xhat1) > OptTol)

%while  (iter <= maxIters)
    xhat1new = xhat1;
    if fastop
        res = y - feval(A,1,n,p,xhat,Ifull,p);
	if ~tightFrame
		corr = lsqrms(n, p, A, 1:p, p, res, damp, atol, btol, conlim, itnlim, show);
	else
		corr = feval(A,2,n,p,res,Ifull,p)/tightFrame;
	end
    else
    	res = y - A*xhat;
        corr = iA*res;
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
	plot(l1norm);
    end

    %if fullPath
    %    sols = [sols xhat1];
    %end

    iter = iter+1;
end

if (lssolution)
    if lambdaStop
    	xhat = LSsolution(A, y, xhat1, p, find(abs(xhat1) > lambdaStop));
    else
    	xhat = LSsolution(A, y, xhat1, p);
    end
    %if fullPath
    %    sols = [sols xhat];
    %end
else
    xhat = xhat1;
end

if (verbose)
        fprintf('Iteration %d: |I| = %d, ||x||_1 = %g\n', iter, length(activeSet), norm(xhat,1));
	l1norm(iter) = norm(xhat,1);
	plot(l1norm);
end

%if fullPath
%    xhat = sols;
%end

numIters = iter;

