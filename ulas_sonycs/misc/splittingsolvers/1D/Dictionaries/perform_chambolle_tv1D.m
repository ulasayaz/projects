function x = perform_chambollefast_tv1D(y,lambda,tvoptions)
% perform_chambolle_tv: Chambolle projection algorithm for TV-regularized minimization problem (ROF):
%		min 0.5|| y - x ||_2 + lambda || grad(x) ||_1 
%
%
%	The solution is given by the iteration :
%		u^(t+1/2) = u^(t) + mu grad[(div u^(t) - y/lambda)]
%		u^(t+1)   = u^(t+1/2) / ( 1 + mu |grad[(div u^(t) - y/lambda)]|),
%
%		with mu <= 1/4 (but too pessimistic in practice, convergence speed is better for 1/2).
%	The solution to the minimization problem is given by x = lambda div(u) at convergence.
%
% Usage
%	x = perform_chambolle_tv1D(y,lambda,tvoptions)
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
	u     = (u + tvoptions.mu*v) ./ (1 + tvoptions.mu*abs(v));
	
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
