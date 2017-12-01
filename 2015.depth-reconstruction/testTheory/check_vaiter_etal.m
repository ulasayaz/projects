%% NOTE:
% D* = TV2
% Phi = R



zlinftest = zlinf;
%zlinftest = zGT;
%zlinftest = zNaiveTight;
%c = 0.8;
%zlinftest = c * zlinf + (1-c) * zlinf_s2;

%% find a solution which curvy only around the critical points

[Itest,Jtest,DItest,DJtest,sItest] = getVaiterMat(TV2,zlinftest,tol);

cardTest = length(Itest);
curvy = [];
for i =1:cardTest
    xright = find(samples >= Itest(i),1, 'first');
    
    xleft = find(samples < Itest(i) ,1,'last');
    
    %curvy = union(curvy,samples(xleft):samples(xright));
    xleft = Itest(i) - 2 ;
    xright = Itest(i) + 2;
    curvy = union(curvy,xleft : xright);
end

flat = setdiff(1:N-2,curvy);

curvy'

cvx_begin
%cvx_quiet true
cvx_precision best
variable z(N,1);
minimize( norm( TV2 * z,1) )
subject to
norm( TV2 * z,1) <= tvMin;
norm(y - R * z, Inf) <= epsilon;
%TV2(curvy,:) * z  >= tol
%TV2(curvy,:) * z  >= tol
TV2(flat,:) * z == 0;
cvx_end

norm(TV2 * z,1)

[Iz,Jz,DIz,DJz,sIz] = getVaiterMat(TV2,z,tol);

figure; clf
plot(zGT,'-g', 'LineWidth', 3); hold on;
plot(z,'-.k', 'LineWidth', 3); 
plot(samples,y,'o','markersize',10)
legend('zGT','z')
grid minor

%%

%zlinftest = z;

if N > 1000
    disp('skipped vaiter check since N is large!!!!!!')
else
    disp('starting check Vaiter')
    
    %% check condition (H0 and Hj)
    % zGT -> I,J, DI, DJ ; zlinf -> Iopt, Jopt, DIopt, DJopt 
    [I,J,DI,DJ,sI] = getVaiterMat(TV2,zGT,tol);
    Phi = R;
    
    ker_Phi = null(full(Phi));
    ker_Dstar = null(full(TV2)); % size (N-2) x N
    U = null(full(DJ')); % size = N x (N-|J|)  % cospace
    % size of DJ = N x |J|
    
    % Kernel(Phi) \intersect Kernel(D*) = {0}
    intKer = null([ker_Phi  ker_Dstar]);
    if isempty(intKer)==0
        warning('Condition H0 is not satisfied')
    end
    % Kernel(Phi) \intersect Kernel(DJ*) = {0}
    intKerJ = null([ker_Phi  U]);
    if isempty(intKerJ)==0
        warning('Condition HJ is not satisfied')
    end
    
    %% check splitting argument: xhat = x0 + O(w)
    [Iopt,Jopt,DIopt,DJopt,sIopt] = getVaiterMat(TV2,zlinftest,tol);
    Uopt = null(full(DJopt')); % size = N x (N-|J|)  % cospace
    
    size(Phi * Uopt)
    
    % little test
    C = length(Iopt)+3; C = 199;
    Tsupp = randperm(N,C); %Tsupp = samples;
    G = Uopt(Tsupp,:);
    min(svd(full(G' * G)))
    
    
    K = length(samples);
    % find active set
    v = (Phi * zlinftest - y);
    v_max = norm(v , inf); % this should be equal to epsilon, if some constraint it tight
    if ~isequalWithTol(v_max,epsilon,0.01,1e-5) warning('constraint is violated'); end 
    
    activeSet = find( abs(v) > v_max - 1e-4 & abs(v) < v_max + 1e-4 ); % close to y+eps
    inactiveSet = setdiff([1:length(v)],activeSet);
   
    
    % check Jopt \subset J
    if ~all(ismember(Jopt, J))
        warning('Jopt is not a subset of J')
        extraEntries = setdiff(Jopt,J);
        JoptReduced = setdiff(Jopt, extraEntries); % add the missing entries
        if ~all(ismember(JoptReduced, J)) 
           error('Jopt is not a subset of J after adding entries') 
        end
    end
    
    % Kernel(Phi) \intersect Kernel(DJopt*) = {0}
    intKerJopt = null([ker_Phi  Uopt]);
    if isempty(intKerJopt)==0
        warning('Condition HJopt is not satisfied')
    end
    
    %% find sigma and alpha (a) that witness optimality of zinf
    cvx_begin
    cvx_quiet true
    cvx_precision best
    variable u(length(Jopt),1)
    variable a(K,1)
    minimize( norm( DIopt * sIopt + DJopt * u + 2 * ve * Phi' * diag(a) * v) )
    subject to
    norm( u, inf ) <= 1;
    sum(a) == 1;
    a >= 0;
    a(inactiveSet) == 0;
    cvx_end
    
    if norm( DIopt * sIopt + DJopt * u + 2 * ve * Phi' * diag(a) * v) > 1e-4
        error('not able to find alpha and sigma to witness optimality')
    end
    
    %% Check Theorem 1:
    if (norm(DJ' * U) > 1e-6) error('wrong null space for DJ'); end
    if (norm(DJopt' * Uopt) > 1e-6) error('wrong null space for DJopt'); end
    if (norm(DJopt' * U) > 1e-6) warning('wrong claim on Jopt \subset J => xgt=U*b'); end
    if (norm(DJopt' * zGT) > 1e-6) warning('xgt not in null space of DJopt'); end
    if (norm(DJ' * zGT) > 1e-6) error('xgt not in null space of DJ'); end
    
    %% try to get to xopt = x0 + .. f(alpha)
    %  K = length(samples);
    %  r = K-1; % floor(K * 1); % the size of active set
    %  pos = randsample(1:K,r);
    %  a = rand(r,1); a = a/norm(a,1); a = abs(a); % create the convex hull coefficients
    %  alpha = zeros(K,1); alpha(pos) = a;
    D_a = diag(a); % diagonal matrix with coefficients alpha
    
    M_Uopt_noalpha = full(Uopt' * Phi' * Phi * Uopt); % = invertible if HJopt satisfied
    M_Uopt_alpha = full(Uopt' * Phi' * D_a * Phi * Uopt);
   
    M_U_noalpha = full(U' * Phi' * Phi * U); % = invertible if HJ satisfied
    M_U_alpha = full(U' * Phi' * D_a * Phi * U); 
    
    if min(svd(M_U_noalpha)) < 1e-4
        warning('M_U_noalpha not invertible') 
    end
    if min(svd(M_Uopt_alpha)) < 1e-4 
        warning('M_Uopt_alpha not invertible')
    end
    if min(svd(M_Uopt_noalpha)) < 1e-4 
        warning('M_Uopt_noalpha not invertible')
    end
    if min(svd(M_U_alpha)) < 1e-4 
        warning('M_U_alpha not invertible')
    end
    
    A_J_alpha = Uopt * inv(M_Uopt_alpha) * Uopt';
    error_term = A_J_alpha * (Phi' * D_a * noise - DIopt * sIopt);
    mismatch = norm( zlinftest - zGT - error_term ) ;
    if isnan(mismatch) || mismatch > 1e-4
        warning('splitting of the signal and error does not work');
    end
    
    %% find the optimal perturbation to make M_Uopt_alpha invertible
    beta = Uopt' * zlinftest;
    if norm(zlinftest - Uopt * beta) > 1e-4
        warning('coefficients beta calculated wrongly!')
    end
    
    dim = length(beta);
    Id = eye(dim);

     sig = 3;
     norm( Uopt * inv(M_Uopt_alpha + sig*Id) * (sig*Id) * beta )
    
     
%     cvx_begin
%     cvx_quiet true
%     cvx_precision best
%     variable sig
%     minimize( norm( Uopt * inv(M_Uopt_alpha + sig*Id) * (sig*Id) * beta ) )
%     cvx_end
   
    
    %% OLD checks from Vaiter

    A_J_ = U * inv(M_U_noalpha) * U' ;
    Omega_J_ = pinv(full(DJ)) * (Phi' * Phi * A_J_ - eye(N) ) * DI;
    
    
    cvx_begin
    cvx_quiet true
    cvx_precision best
    variable u(length(J),1);
    minimize( norm( Omega_J_ * sI - u, inf ) )
    subject to
    norm( DJ * u ) <= 1e-6;
    % DJ * u == zeros(N,1); % infeasible, for numerical reasons
    cvx_end
    
    ICs = norm( Omega_J_ * sI - u, inf )
    % RC is hard to compute: requires solving nonconvex problem
    wRCI = norm(Omega_J_,inf)
    % recall that: IC < RC < wRC
    
    % cvx_quiet false
    % cvx_begin
    % cvx_precision best
    % variable u2(size(U,2),1);
    % minimize( norm( Omega_J_ * sI - U * u2, inf ) )
    % cvx_end
    
    %% Check Lemma 1:
    cvx_begin
    cvx_quiet true
    cvx_precision best
    variable sigm(length(J),1);
    minimize( norm( DJ * sigm + DI * sI ) )
    subject to
    norm(sigm,inf) <= 1;
    cvx_end
    isZeroLemma1 = norm( DJ * sigm + DI * sI )
end
