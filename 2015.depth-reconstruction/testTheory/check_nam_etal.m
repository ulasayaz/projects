% condition (25), only evaluated at zGT
nonSamples = setdiff([1:N],samples);
Rnull = Rfull(nonSamples, :);
if norm(Rnull * R',inf) > 1e-6 || norm(R * Rnull',inf) > 1e-6
    error('wrong null space basis 1')
end
if size(R,2)<100 && rank(full([R; Rnull]))~= N  % only checked if matrices are sufficiently small
    error('wrong null space basis 2')
end

TV2_zGT = abs(TV2 * zGT);
Lambda = find(TV2_zGT < tol);
Omega_Lambda = TV2(Lambda,:);
Lambda_c = setdiff([1:N-2],Lambda);
Omega_Lambda_c = TV2(Lambda_c,:);

if N < 200
    c1 = norm(Omega_Lambda * zGT);
    if c1 > 1e-5
        error('wrong cosupport set')
    end
    eq25 = norm( pinv(full(Rnull * Omega_Lambda')) * (Rnull * Omega_Lambda_c') * sign(Omega_Lambda_c * zGT), inf) % should be less than 1
    eq26 = norm( pinv(full(Rnull * Omega_Lambda')) *  (Rnull * Omega_Lambda_c'), inf)
else
    disp('skipping computation of eq 25 and 26 to save time!!!!!!')
    eq25 = -1;
    eq26 = -1;
end

M1 = Rnull * Omega_Lambda_c';
M2 = pinv(full(Rnull * Omega_Lambda'));

% the following is commented out as it does not work - assumes wrong
% property for pseudoinverse
% if( norm(M2 - pinv(full(Omega_Lambda')) * pinv(full(Rnull)) ) > 1e-6) warning('pseudoinverse check failed 0'); end
% if( norm(M2 - inv(Omega_Lambda * Omega_Lambda') * Omega_Lambda * Rnull' * inv(Rnull*Rnull')) > 1e-6) warning('pseudoinverse check failed 1'); end
% if( norm(M2 - inv(Omega_Lambda * Omega_Lambda') * Omega_Lambda * Rnull') > 1e-6) warning('pseudoinverse check failed 2'); end
% eqMine = norm( inv(Omega_Lambda * Omega_Lambda') *  Omega_Lambda * Rnull' * Rnull * Omega_Lambda_c', inf)

figure
subplot(1,3,1)
spy(Rnull * Rnull'); title('Rn * Rnt')
if (norm(Rnull * Rnull' - eye(N-K)) > 1e-6) error('identity 1 failed'); end
subplot(1,3,2)
spy(Rnull' * Rnull); title('Rnt * Rn')
% This must be false: if (norm(Rnull' * Rnull - eye(N)) > 1e-6) error('identity 2 failed'); end
subplot(1,3,3)
spy(Omega_Lambda * Omega_Lambda'); title('OL * OLt')
if rank(full(Rnull))~=size(Rnull,1)
    error('not full rank Rnull')
end
if rank(full(Omega_Lambda))~=size(Omega_Lambda,1)
    warning('not full rank Omega_Lambda')
end
if abs(norm( M2 * M1, inf) - eq26) > 1e-6
    error('norm check failed')
end
% norm(M1, inf) * norm(M2, inf) % very loose bound

figure
disp('==========================================')
subplot(1,4,1)
spy(M1)
max1 = max(sum(abs(M1),2)); % linf operator norm
title('M1')
nnz1 = nnz(M1)
norm1 = norm(M1, inf)
if (abs(max1 - norm1) > 1e-6) error('norm check failed 1'); end
disp('==========================================')
subplot(1,4,2)
M2(abs(M2)<1e-6) = 0;
spy(M2)
title('M2')
max2 = max(sum(abs(M2),2)); % linf operator norm
nnz2 = nnz(M2)
norm2 = norm(M2, inf)
if (abs(max2 - norm2) > 1e-6) error('norm check failed 2'); end
disp('==========================================')
subplot(1,4,3)
M3= Rnull * Omega_Lambda';
spy(M3)
title('M3 (no pinv)')
max3 = max(sum(abs(M3),2)); % linf operator norm
nnz3 = nnz(M3)
norm3 = norm(M3, inf)
if (abs(max3 - norm3) > 1e-6) error('norm check failed 3'); end
% if (abs(nnz3 - (3 * size(M3,1)-2)) > 1e-6) error('M3 not tridiagonal'); end
disp('==========================================')
subplot(1,4,4)
M4 = M2 * M1;
spy(M4)
title('M2 * M1')
max4 = max(sum(abs(M4),2)); % linf operator norm
nnz4 = nnz(M4)
norm4 = norm(M4, inf)
if (abs(max4 - norm4) > 1e-6) error('norm check failed 4'); end
disp('==========================================')

figure
subplot(1,3,1); hold on
title('Rnull')
spy(Rnull)
subplot(1,3,2); hold on
M5 = Rnull * TV2';
title('Rnull * TV2t')
spy(M5)
subplot(1,3,3); hold on
M5pinv = pinv(full(M5));
title('pinv - Rnull * TV2t')
M5pinv(abs(M5pinv) < 1e-8) = 0;
spy(M5pinv)

TV2t = TV2';
if norm(M5 - TV2t(nonSamples,:),inf) > tol
    error('structure of M5 was not understood correctly')
end
