if N > 500
    disp('skipped duality check since N is large!!!!!!')
else
    
    disp('=============================================')
    disp('check duality')
    
    %% l1 minimization
    cvx_begin
    cvx_quiet true
    cvx_precision best
    variable z(N,1);
    minimize( norm( TV2 * z,1) )
    subject to
    (y - R * z)' * (y - R * z) <= epsilon^2;
    cvx_end
    pstar = norm( TV2 * z,1)
    
    if epsilon > 1e-5
        %% solve dual
        % M = [eye(K)  y; y' epsilon^2];
        % cvx_begin
        % cvx_precision best
        % variable eta_a(N-2,1);
        % variable eta_b(K,1);
        % variable eta_e;
        % minimize( [eta_b; 2*eta_e]' * M * [eta_b; 2*eta_e])
        % subject to
        % norm(eta_a, inf) <= 1;
        % R' * eta_b == - TV2' * eta_a;
        % eta_e > 0;
        % %% ISSUE: M not positive semidefinite
        % cvx_end
        % dstar = (-1/eta_e) * [eta_b; 2*eta_e]' * M * [eta_b; 2*eta_e]
        
        cvx_begin
        cvx_quiet true
        cvx_precision best
        variable eta_a(N-2,1);
        variable eta_b(K,1);
        maximize( - epsilon * norm(eta_b) - eta_b' * y)
        subject to
        norm(eta_a, inf) <= 1;
        R' * eta_b == - TV2' * eta_a;
        cvx_end
        dstar = - epsilon * norm(eta_b) - eta_b' * y
        % the following was eliminated in closed form
        eta_e = norm(eta_b) / (2*epsilon);
        
        %% solve lagrangian
        cvx_begin
        cvx_quiet true
        cvx_precision best
        variable a(N-2,1);
        variable b(K,1);
        variable x0(N,1);
        minimize( norm(a,1) + eta_a' * (TV2 * x0 - a) + eta_e' * (norm(b)-epsilon^2) + eta_b' * (R*x0-y-b) )
        cvx_end
        Lvalue0 = norm(a,1) + eta_a' * (TV2 * x0 - a) + eta_e' * (norm(b)-epsilon^2) + eta_b' * (R*x0-y-b)
        
        %% solve lagrangian, 2
        cvx_begin
        cvx_quiet true
        cvx_precision best
        variable x(N,1);
        minimize( norm(TV2 * x,1) + eta_e * ((R*x-y)' * (R*x-y) - epsilon^2))
        cvx_end
        Lvalue = norm(TV2 * x,1) + eta_e * ((R*x-y)' * (R*x-y) - epsilon^2)
        
        %% solve regularized least squares
        lambda = 1/(2*eta_e);
        cvx_begin
        cvx_quiet true
        cvx_precision best
        variable xls(N,1);
        minimize( 1/2 * (R*xls-y)' * (R*xls-y) + lambda * norm(TV2 * xls,1))
        cvx_end
    else
        x = zeros(N,1);
        xls = zeros(N,1);
    end
    
end