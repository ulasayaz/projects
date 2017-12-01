clear all
close all
clc

N = 9;
[TV2] = createFiniteDiff2(N);

TV2' * ones(N-2,1);

cut = [1 N] % [1 2] % randperm(N,2)
kept = setdiff([1:N],cut);

TV2t = TV2';
TV2tcut = TV2t(kept,:);

v = zeros(N,1);
v(1) = -1;  v(2) = -v(1);  v(end-1) = 1; v(end) = -v(end-1);
v = v(kept)
sol = TV2tcut \ v
if norm(TV2tcut * sol - v) > 1e-7
    warning('system does not admit solution')
end

%full(TV2t)
norm(inv(TV2tcut),inf);

full(inv(TV2tcut))
% pinv(full(TV2'))
cond(TV2tcut)

%% Siegel's lemma: max bound on sol is:
% N = nr of rows of our linear system
% M = nr of rows -1 assuming that we take 3 sample per segment
M = N-1;
B = 2; % max entry in the matrix
(N * B)^(M/(N-M))
D = 1 % I don't know how to set this
(1/D)*(sqrt(det(TV2 * TV2')))^(1/(N-M))


eig(full(TV2 * TV2'))

T = TV2t(2:end-1,:)
v = zeros(size(T,1),1);
v(1)=1;
v(end)=-1
inv(T) * v
invT = inv(T);
full([invT(:,1) invT(:,end)])

% cvx_begin
% cvx_quiet true
% cvx_precision best
% variable sigm(N-2,1);
% minimize( norm( TV2tcut * sigm - v ) )
% subject to
% norm(sigm,inf) <= 1;
% cvx_end
% norm( TV2tcut * sigm - v )

invT
n = size(invT,1)
for i=1:size(invT,1)
    for j=1:size(invT,2)
        if i<=j
          invTactual(i,j) = (-1)^(2*i-1) * (i) * (n-j+1) / (n+1);
        else
          invTactual(i,j) = (-1)^(2*j-1) * (j) * (n-i+1) / (n+1);
        end
    end
end
invTactual

norm(invT -invTactual )

TV2t = full(TV2')
TVblock= full(TV2t(end-3:end,end-3:end))
inv(TVblock)
v = zeros(size(TVblock,1),1);
v(1)=1;
v(2)=-1
TVblock * TVblock