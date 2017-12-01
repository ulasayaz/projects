function x = perform_nesterov_tv2D(y,lambda,tvoptions)
% perform_nesterov_tv: Nesterov algorithm for TV-regularized minimization problem (ROF):
%		min 0.5|| y - x ||_2 + lambda || grad(x) ||_1
% The proximity operator of TV is computed through the Fenchel dual problem.
%
%
%	The solution is given by the iteration :
%		prox_TV(y) = (Id - prox_TV^*)(y) = y - P_C(y),
%	where C={div u: |u| <= 1}.
%	The iteration converges to u such that div u = P_C(y).
%	The solution to the minimization problem is then given by x = y - lambda div(u) at convergence.
%
% Usage
%	x = perform_prox_tv2D(y,lambda,tvoptions)
% Input
%	y           nxn image
%	lambda	    regularization parameter
%	tvoptions   a structure containing options for TV prox subiterations:
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

u    = 0;
u0   = 0;
At   = 0;
zeta = zeros(2,n,n);
mu   = 0.25;
yn   = y/lambda;
%unew = 0;
%xnew = 0;
%lambda = n*sigma;

for iter=1:tvoptions.iterstv
	v     = u0-zeta;
	den   = sqrt(sum(v.^2,1));
	v(1,den > 1) = v(1,den > 1) ./ den(den > 1)';
	v(2,den > 1) = v(2,den > 1) ./ den(den > 1)';
	at    = (mu+sqrt(mu*mu+4*mu*At))/2;
	w     = (At*u + at*v)/(at+At);
	u     = w - mu/2*perform_grad(yn - perform_div(w,1),1);
	den   = sqrt(sum(u.^2,1));
	u(1,den > 1) = u(1,den > 1) ./ den(den > 1)';
	u(2,den > 1) = u(2,den > 1) ./ den(den > 1)';
	At    = At+at;
	x     = y - lambda*perform_div(u,1);
	zeta  = zeta + at*perform_grad(x/lambda,1);
	
	%xnew = y - lambda*yproj;
	%if norm(xnew(:) - xold(:)) <= tvoptions.tvtol, break; end
	%xold = xnew;
	%xd = perform_grad(x,1);
	%energ(iter,:) = [0.5*norm(x - y,'fro')^2+lambda*sum(sum(sqrt(sum(xd.^2,1)))) 0.5*norm(x,'fro')^2 - 0.5*norm(y(:))^2 sqrt(sum(sum(sum((u - unew).^2,1)))) norm(x - tvoptions.xinf,'fro')];
	%unew = u;
	%xnew = x;
end

x = x(:);

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
