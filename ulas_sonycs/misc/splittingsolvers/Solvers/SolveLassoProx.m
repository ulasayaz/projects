function [xhat, numIters, activeSet] = SolveLassoProx(A, y, p, q, mu, lambdaStop, maxIters, lssolution, fullPath, verbose, positivity)
% SolveLassoProx: Iterative Backward-Forward Splitting Proximal iteration to solve the LASSO problem:
%		min || y - Ax ||_2^2 s.t. || x ||_1 <= q
%	The solution is given by the iteration (a projected gradient descent algorithm in this case):
%		x^(t+1) = P_q [x^(t) + mu_t (y - Ax^(t))]
%		where P_q is the projector onto the l_1-ball of radius q.
%	To implement P_q, many methods based on soft-thresholding with variable threshold: 
%			interpolation (exact), Uzawa (exact), dichotomic search (approximate).
%
% Usage
%	[xhat, numIters, activeSet] = SolveLassoProx(A, y, p, mu, q, maxIters, fullPath, verbose, positivity)
% Input
%	A           nxp matrix or implicit operator
%	y           n-vector
%  	p           length of x
%   	mu          relaxation parameter 0 < inf mu_t <= sup mu_t < 2/||A||^2
%   	lambdaStop  If specified (and > 0), the algorithm removes all coeffs <= lambdaStop in magnitude for the LS solution
%	maxIters
%	lssolution  return the least-squares estimate A_I^+y once the support is obtained    
%   	fullPath    1 returns entire solution path, 0 final
%   	verbose     
%	positivity  enforces positivity of the solution.
% Outputs
%	 sols       solution of LASSO problem
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
if ~exist('mu'), 	 mu = 1;		end
if ~exist('q'), 	 q = max(abs(y(:))); 	end
if ~exist('p'), 	 p = size(A,2); 	end

xhat = zeros(p,1);
activeSet = [];
Ifull = 1:p;
iter = 1;
mu = (mu/(max(q)^2));
sols = [];

%while  (iter <= maxIters) & (norm(res) > OptTol*normy)

while  (iter <= maxIters)
    
    if fastop
        res = y - feval(A,1,n,p,xhat.*q,Ifull,p);
        corr = q.*feval(A,2,n,p,res,Ifull,p);
    else
    	res = y - A*(xhat.*q);
        corr = q.*(A'*res);
    end
    xhat = xhat + mu*corr;
    [xhat,lambda] = perform_l1ball_projection(xhat,1);
    
    if positivity
  	xhat(find(xhat < 0)) = 0;
    end
    
    activeSet = find(abs(xhat) > eps);
    
    if (verbose)
        fprintf('Iteration %d: |I| = %d, lambda = %g, ||r|| = %g\n', iter, length(activeSet), lambda, norm(res));
	plot(xhat);drawnow
	%rss(iter) = norm(res);
	%plot(rss);drawnow
    end

    if fullPath
        sols = [sols xhat];
    end

    iter = iter+1;
end

xhat = xhat.*q;

if (lambdaStop)
    	xhat(abs(xhat) <= lambdaStop) = 0;
end

if (lssolution)
    if lambdaStop
    	xhat = LSsolution(A, y, xhat, p, find(abs(xhat) > lambdaStop));
    else
    	xhat = LSsolution(A, y, xhat, p, activeSet);
    end
    if fullPath
        sols = [sols xhat];
    end
end

if (verbose)
        fprintf('Iteration %d: |I| = %d, lambda = %g, ||r|| = %g\n', iter, length(activeSet), lambda, norm(res));
	rss(iter) = norm(res);
	plot(rss);
end

if fullPath
    xhat = sols;
end

numIters = iter;

