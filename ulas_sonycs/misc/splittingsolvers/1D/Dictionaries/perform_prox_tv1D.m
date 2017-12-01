function x = perform_prox_tv1D(y,lambda,tvoptions)
% perform_prox_tv: Projected-gradient algorithm for TV-regularized minimization problem (ROF):
%		min 0.5|| y - x ||_2 + lambda || grad(x) ||_1
% The proximity operator of TV is computed through that of its Legendre-Fenchel conjugate which amounts 
% to a projected gradient iteration.
%
%
%	The solution is given by the iteration :
%		prox_TV(y) = (Id - prox_TV^*)(y) = y - P_C(y),
%	where C={div u: |u| <= 1}.
%	The projected (onto C) gradient iteration is given by:
%		u^(t+1)   = P_W[u^(t) + mu_t grad(div u^(t) - y/lambda)],
%		where W={u: |u| < 1}, and 0 < inf_t mu_t <= sup_t mu < 1/2.
%	and P_W(u) = u if u \in W and u/|u| otherwise.
%	The iteration converges to u such that div u = P_C(y).
%	The solution to the minimization problem is then given by x = y - lambda div(u) at convergence.
%
% Usage
%	x = perform_prox_tv1D(y,lambda,tvoptions)
% Input
%	y           n-vector
%	lambda	    regularization parameter
%	tvoptions   a structure containing options for TV prox subiterations:
%   			mu	    relaxation parameter
%			iterstv
%			tol	    stops when the ||x^(t+1) - x^(t)|| <= tvtol
%   			verbose     
% Outputs
%	 x	    solution of the poroblem
%


n = length(y);
y = y(:);

%if ~exist('lambda'),	 lambda = 1;		end
%if ~exist('tvoptions'),  
%	tvoptions.mu = 0.5;
%	tvoptions.iterstv = 50;
%	tvoptions.tvtol = 1E-3;
%	tvoptions.verbose = 0;
%else
%	if ~isfield(tvoptions,'mu'),	     tvoptions.mu = 0.5;     end
%	if ~isfield(tvoptions,'iterstv'),    tvoptions.iterstv = 50;  end
%	if ~isfield(tvoptions,'tvtol'),      tvoptions.tvtol = 1E-3;  end
%	if ~isfield(tvoptions,'verbose'),    tvoptions.verbose = 0;   end
%end

u = y/norm(y);
xold = 0;
%lambda = sqrt(n)*sigma;

for iter=1:tvoptions.iterstv
	yproj = perform_divvector(u,1);
	v     = perform_gradvector(yproj - y/lambda,1);
	u     = (u + tvoptions.mu*v);
	u     = u - (abs(u) > 1) .* (u - sign(u));
	
	%xnew = y - lambda*yproj;
	%if norm(xnew(:) - xold(:)) <= tvoptions.tvtol, break; end
	%xold = xnew;
	
	%if tvoptions.verbose,
	%	x = y - lambda*yproj;
	%	tvnorm(iter) = norm(perform_gradvector(x,1),1);
	%	subplot(211)
	%	plot(norm(y-x)^2+lambda*tvnorm);drawnow
	%	subplot(212)
	%	plot(x);drawnow
	%end
	% The following lines may be useful for the constrained problem.
	  %yproj = perform_divvector(u);
	  %lambda = sqrt(n)*sigma / norm(yproj);
end

x = y - lambda*yproj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = perform_gradvector(u,h)
%  Computes the forward finite differences, ie.
%           ud(i) = (u(i+1)-u(i))/(h), with special
%           care at boundaries (ud(n) = 0).
%

ud = [diff(u);0]/h;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = perform_divvector(u,h)
%  Computes the backward finite differences, ie.
%           ud(i) = (u(i)-u(i-1))/(h), with special
%           care at boundaries (ud(1) = u(1) and ud(n) = -u(n-1)).
%

ud = [u(1);diff(u(1:end-1));-u(end-1)]/h;
