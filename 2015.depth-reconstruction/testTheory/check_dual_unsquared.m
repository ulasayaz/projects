% %% l1 minimization
cvx_quiet true
cvx_begin
cvx_precision best
variable z(N,1);
minimize( norm( TV2 * z,1) )
subject to
norm(y - R * z) <= epsilon;
cvx_end
pstar = norm( TV2 * z,1)

%% solve dual
cvx_begin
cvx_precision best
variable eta_a(N-2,1);
variable eta_b(K,1);
variable eta_e;
maximize( -eta_b' * y - eta_e * epsilon)
subject to
norm(eta_a, inf) <= 1;
norm(eta_b, 2) <= eta_e;
R' * eta_b == - TV2' * eta_a;
cvx_end
dstar = -eta_b' * y - eta_e * epsilon

%% solve regularized least unsquared norm
cvx_begin
cvx_precision best
variable a(N-2,1);
variable b(K,1);
variable x(N,1);
minimize( norm(a,1) + eta_a' * (TV2 * x - a) + eta_e' * (norm(b)-epsilon) + eta_b' * (R*x-y-b) )
cvx_end
%% ISSUE: cannot be applied if this lacks uniqueness
%x = 1e4 * rand(N,1);
Lvalue = norm(a,1) + eta_a' * (TV2 * x - a) + eta_e' * (norm(b)-epsilon) + eta_b' * (R * x - y -b) 