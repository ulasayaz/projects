function x = perform_chambolle_tv2D(y,lambda,tvoptions)
% perform_chambolle_tv: Chambolle projection algorithm for TV-regularized minimization problem (ROF):
%		min 0.5|| y - x ||_2 + lambda || grad(x) ||_1 
%
%
%	The solution is given by the iteration :
%		u^(t+1/2) = u^(t) + mu grad[(div u^(t) - y/lambda)]
%		u^(t+1)   = u^(t+1/2) / ( 1 + mu |grad[(div u^(t) - y/lambda)]|),
%
%		with mu <= 1/8 (but too pessimistic in practice, convergence speed is better for 1/4).
%	The solution to the minimization problem is given by x = lambda div(u) at convergence.
%
% Usage
%	x = perform_chambolle_tv2D(y,lambda,tvoptions)
% Input
%	y           nxn image
%	lambda	    regularization parameter
%	tvoptions   a structure containing options for TV prox subiterations:
%   			mu	    relaxation parameter
%			iterstv
%			tvtol	    stops when the ||x^(t+1) - x^(t)|| <= tvtol
%   			verbose     
% Outputs
%	 x	    solution of the problem
%


n = floor(sqrt(prod(size(y))));
y = reshape(y, n, n);

%if ~exist('lambda'),	 lambda = 1;		end
%if ~exist('tvoptions'),  
%	tvoptions.mu = 0.25;
%	tvoptions.iterstv = 50;
%	tvoptions.tvtol = 1E-3;
%	tvoptions.verbose = 0;
%else
%	if ~isfield(tvoptions,'mu'),	     tvoptions.mu = 0.25;     end
%	if ~isfield(tvoptions,'iterstv'),    tvoptions.iterstv = 50;  end
%	if ~isfield(tvoptions,'tvtol'),      tvoptions.tvtol = 1E-3;  end
%	if ~isfield(tvoptions,'verbose'),    tvoptions.verbose = 0;   end
%end

u = zeros(2,n,n);
xold = 0;
%lambda = n*sigma;

for iter=1:tvoptions.iterstv
	yproj = perform_div(u,1);
	v     = perform_grad(yproj - y/lambda,1);
	den   = 1 + tvoptions.mu*sqrt(sum(v.^2,1));
	u(1,:,:) = (u(1,:,:) + tvoptions.mu*v(1,:,:)) ./ den;
	u(2,:,:) = (u(2,:,:) + tvoptions.mu*v(2,:,:)) ./ den;
	
	%xnew = y - lambda*yproj;
	%if norm(xnew(:) - xold(:)) <= tvoptions.tvtol, break; end
	%xold = xnew;
	
	% The following lines may be useful for the constrained problem.
	  %yproj = perform_div(u,1);
	  %lambda = n*sigma / norm(yproj(:));
end

x = y(:) - lambda*yproj(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = perform_grad(u,h)
%  Computes the forward finite differences, ie.
%	    ud(1,i,j) = (u(i+1,j)-u(i,j))/(h),
%	    ud(2,i,j) = (u(i,j+1)-u(i,j))/(h), with special
%           care at boundaries (ud(n) = 0).
%	    returns |ud|.
%

[n,m] = size(u);
ud = zeros(2,n,m);
ud(1,:,:) = [diff(u,1,1);zeros(1,m)]/h;
ud(2,:,:) = [diff(u,1,2),zeros(n,1)]/h;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = perform_div(u,h)
%  Computes the backward finite differences, ie.
%           ud(i,j) = (u(i,j)-u(i-1,j))/(h) + (u(i,j)-u(i,j-1))/(h), with special
%           care at boundaries (ud(1,j) = u(1,j) and ud(n,j) = -u(n-1,j)),
%			       (ud(i,1) = u(i,1) and ud(i,n) = -u(i,n-1))
%

u1 = squeeze(u(1,:,:));
u2 = squeeze(u(2,:,:));
ud = [u2(:,1),diff(u2(:,1:end-1),1,2),-u2(:,end-1)]/h ...
    +[u1(1,:);diff(u1(1:end-1,:),1,1);-u1(end-1,:)]/h;
