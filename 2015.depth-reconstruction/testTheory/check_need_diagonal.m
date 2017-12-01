clear all
close all
clc
addpath('../')
addpath('../../myLib')
addpath('./../lib')
addpath('./../testTheory')
addpath('../../libOPL')
addpath('./../testTheory')

N = 3;
[TV2] = createFiniteDiff2(N);
I3 = eye(3);
maxTests = 10000;

M = full([kron(I3,TV2) ; kron(TV2,I3)])
nrCol = 0;
for test = 1:maxTests
    toFix = randperm(N^2,3);
    Mcut = M(:,setdiff([1:N^2], toFix));
    %% skip cases with collinear points
    if isCollinear3x3(toFix, N)
       nrCol = nrCol+1;
       continue
    end
    %% when non collinear, check rank
    if rank(Mcut) == min(size(Mcut))
       error('found invertible instance') 
    end
end
fprintf('nrCollinear instances: %d (over: %d)\n',nrCol, maxTests)

Mdiag = [M ; [-1/4 0 1/4, 0 0 0, 1/4 0 -1/4]]
nrColDiag = 0;
for test = 1:maxTests
    toFix = randperm(N^2,3);
    Mcut = Mdiag(:,setdiff([1:N^2], toFix));
    %% skip cases with collinear points
    if isCollinear3x3(toFix, N)
        nrColDiag = nrColDiag+1;
       continue
    end
    %% when non collinear, check rank
    if rank(Mcut) < min(size(Mcut))
       warning('found non invertible instance with diag (i)')
       test, I, J
       pause
    end
end
fprintf('nrCollinear instances: %d (over: %d)\n',nrColDiag, maxTests)

% Mdiag = [M ; [0 -1/2 1/2, -1/2 1 -1/2, 1/2 -1/2 0]]
% for test = 1:5000
%     toFix = randperm(N^2,3);
%     Mcut = Mdiag(:,setdiff([1:N^2], toFix));
%     if rank(Mcut) < min(size(Mcut))
%        warning('found non invertible instance with diag (ii)') 
%        break
%     end
% end