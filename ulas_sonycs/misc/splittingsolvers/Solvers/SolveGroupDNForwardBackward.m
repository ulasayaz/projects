function [xhat, numIters] = SolveGroupDNForwardBackward(A, y, p, blockoptions, lambda, mu, maxIters, fullPath, verbose, positivity)
% SolveGroupDNForwardBackward: Iterative Forward-Backward Splitting Proximal iteration to solve the augmented Lagrangian group-sparsity minimization problem:
%		min 0.5|| y - A x ||_2^2 + lambda sum_b || x_b ||_2 , b=1,...,B, x_b is the restriction of x to the bth block.
%
%	The solution is given by the iteration :
%		x^(t+1)   = ST_{mu_t lambda||.||_2}[x^(t) + mu_t A'(y-Ax^(t))],
%
%		0 < mu_t < 2/||A||^2.
%		p=ST_{mu_t lambda||.||_2}[x] is close to the Stein block-estimator (proximal operator of the block-l2 norm) whose solution is:
%		
%				p_b = x_b max(0,1-lambda/||x_b||_2).
%
%	The solution to the minimization problem is given by ST_(x^(t)) at convergence.
%
% Usage
%	[xhat, numIters] = SolveGroupDNForwardBackward(A, y, p, blockoptions, lambda, mu, maxIters, fullPath, verbose, positivity)
% Input
%	A           nxp matrix or implicit operator
%	y           n-vector
%  	p           length of x
%	blockoptions a structure containing options for block-shrinkage:
%			L 		the number of coeffs in the bth block
%			B 		the number of blocks
%		    			L and B are such that L*B=p
%			ind		L x B matrix of block indices
%	lambda	    regularization parameter
%   	mu	    relaxation parameter, 0 < mu_t < 2/||A||^2
%	maxIters
%   	fullPath    1 returns entire solution path, 0 final
%   	verbose     
%	positivity  enforces positivity of the solution.
% Outputs
%	 xhat       solution of the problem
%	 numIters
%
global dict

fastop = (ischar(A) || isa(A, 'function_handle'));

n = length(y);

if ~exist('positivity'), positivity = 0; 	end
if ~exist('verbose'), 	 verbose = 0; 		end
if ~exist('fullPath'), 	 fullPath = 0; 		end
if ~exist('maxIters'), 	 maxIters = 50; 	end
if ~exist('mu'), 	 mu = 1;		end
if ~exist('epsilon'), 	 epsilon = 0;		end
if ~exist('p'), 	 p = size(A,2); 	end
if ~exist('blockoptions') | isempty(blockoptions),
	blockoptions.L = 1;
	blockoptions.B = p;
	blockoptions.ind = 1:p;
end

% lsqrms params
damp = 0;
atol = 1e-6;
btol = 1e-6;
conlim = 1e+10;
itnlim = 10;
show = 0;

xhat = zeros(p,1);
Ifull = 1:p;
iter = 1;
sols = [];

%while  (iter <= maxIters) & (norm(res) > OptTol*normy)

while  (iter <= maxIters)
    
    if fastop
        res = y - feval(A,1,n,p,xhat,Ifull,p);
        corr = feval(A,2,n,p,res,Ifull,p);
    else
    	res = y - A*(xhat);
        corr = A'*res;
    end
    xhat = xhat + mu*corr;
    if blockoptions.L>1, ener = sqrt(sum(abs(xhat(blockoptions.ind)).^2,1));
    else		 ener = abs(xhat);
    end 
    mask = max(0,1-mu*lambda./ener);
    xhat = xhat(blockoptions.ind).*repmat(mask,blockoptions.L,1);
    xhat = xhat(:);
    
    if positivity
  	xhat(find(xhat < 0)) = 0;
    end
      
    if (verbose)
    	l2norm(iter) = lambda*sum(ener);
	rss(iter) = norm(res)^2/2;
	subplot(211);
	stem(xhat,'+r');drawnow
	subplot(212)
	plot([rss' l2norm' rss'/2+lambda*l2norm']);legend('||residual||_2^2','||.||_2','Objective');drawnow
    end

    if fullPath
        sols = [sols xhat];
    end

    iter = iter+1;
end


if (verbose)
    	l2norm(iter) = lambda*sum(sqrt(ener));
	rss(iter) = norm(res)^2/2;
	plot([rss' l2norm' rss'/2+lambda*l2norm']);legend('||residual||_2^2','||.||_2','Objective');drawnow
end

if fullPath
    xhat = sols;
end

numIters = iter;

