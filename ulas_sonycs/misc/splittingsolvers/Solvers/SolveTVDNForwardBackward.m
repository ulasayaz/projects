function [xhat, numIters] = SolveTVDNForwardBackward(A, y, p, tvoptions, lambda, mu, maxIters, fullPath, verbose, positivity)
% SolveTVDNForwardBackward: Iterative Forward-Backward Splitting Proximal iteration to solve the augmented Lagrangian TV minimization problem:
%		min 0.5|| y - A x ||_2^2 + lambda|| x ||_TV
%
%	The solution is given by the iteration :
%		x^(t+1)   = prox_TV[x^(t) + mu_t A'(y-Ax^(t))],
%
%		0 < mu_t < 2/||A||^2.
%		prox_TV[x] is the proximal operator of the TV norm (computed through Legendre-Fenchel duality and Moreau decomposition
%			   In 1D, it is implemented by a forward-backward iteration. In 2D, fast algorithm of CD).
%%
%	The solution to the minimization problem is given by prox_iC(x^(t)) at convergence.
%
% Usage
%	[xhat, numIters] = SolveTVDNForwardBackward(A, y, p, tvoptions, lambda, mu, maxIters, fullPath, verbose, positivity)
% Input
%	A           nxp matrix or implicit operator
%	y           n-vector
%  	p           length of x
%	tvoptions   a structure containing options for TV prox subiterations:
%			name		name of the algorithm that solves for the prox 
%					forward-backward, Chambolle, fast Chambolle-Darbon (for 2D).
%			dimension	1 for 1D vectors and 2 for 2D images (sqrt(p) x sqrt(p) images)
%			numdeep		depth of the dyadic search for the fast TV prox (only for 2D) (optional)
%			lmin,lmax	the solution is truncated to lmin/lmax (optional)
%			mu		relaxation parameter for the forward-backward implementation case (optional)
%			iterstv		number of iterations for the 1D case (optional)
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
if ~exist('tvoptions'),  
	tvoptions.dimension = '1';
	tvoptions.name = 'perform_prox_tv';
	tvoptions.mu = 0.4;
	tvoptions.iterstv = 200;
	tvoptions.tvtol = 1E-3;
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
if tvoptions.dimension == '1'
   if ~isfield(tvoptions,'name'), 	tvoptions.name = 'perform_prox_tv'; end
   if ~isfield(tvoptions,'mu'), 	tvoptions.mu = 0.4;  end
   if ~isfield(tvoptions,'iterstv'), 	tvoptions.iterstv = 200; end
   if ~isfield(tvoptions,'tvtol'), 	tvoptions.tvtol = 1E-3;  end
   if ~isfield(tvoptions,'verbose'), 	tvoptions.verbose = 0;   end
else
   if ~isfield(tvoptions,'name'), 	tvoptions.name = 'perform_chambollefast_tv'; end
   if ~isfield(tvoptions,'numdeep'), 	tvoptions.numdeep = 8; 	    end
   if ~isfield(tvoptions,'lmin'), 	tvoptions.lmin = min(y(:)); end
   if ~isfield(tvoptions,'lmax'), 	tvoptions.lmax = max(y(:)); end

% Options for old forward-backward code.
% if ~isfield(tvoptions,'mu'),	       tvoptions.mu = 0.25;      end
% if ~isfield(tvoptions,'iterstv'),    tvoptions.iterstv = 400; end
% if ~isfield(tvoptions,'tvtol'),      tvoptions.tvtol = 1E-3;  end
end

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
    if isreal(xhat(1)) 
    	xhat = eval([tvoptions.name tvoptions.dimension 'D(xhat,mu*lambda,tvoptions)']);
    else
        disp('The complex case is not implemented yet for TV minimization.');
    	return
    end
    
    if positivity
  	xhat1(find(xhat1 < 0)) = 0;
    end
      
    if (verbose)
    	tvnorm(iter) = norm(perform_gradvector(real(xhat)),1);
        fprintf('Iteration %d: ||x||_TV = %g\n', iter, tvnorm(iter));
	rss(iter) = norm(res);
	plot(rss);drawnow
    end

    if fullPath
        sols = [sols xhat];
    end

    iter = iter+1;
end


if (verbose)
	fprintf('Iteration %d: |I| = %d, ||c||_1 = %g\n', iter, norm(C,1));
	l1norm(iter) = norm(C,1);
	plot(l1norm);
end

if fullPath
    xhat = sols;
end

numIters = iter;

