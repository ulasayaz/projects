% d is the number or rows = number of columns in the image
dd = d^2;
if size(TV2,1) ~= 2*(d^2 - 2*d) || size(TV2,2) ~= dd 
    error('misunderstood size of TV2')
end
boundCols = blockToMatIndices([1 2 d-1 d], d);
boundRows = [1:d:dd, 2:d:dd, d-1:d:dd, d:d:dd];
rowsToCut = union(boundCols,boundRows);
TV2t = TV2';
rowsTokeep = setdiff([1:dd],rowsToCut);

Delta = full(TV2t(rowsTokeep,:));

if size(Delta,1) ~= dd - 8*d+16 || size(Delta,2) ~= 2*(d^2 - 2*d)
    error('misunderstood size of Delta')
end

nnzCols = find(sum(abs(Delta)) > 1e-5);

Delta_nnz = Delta(:,nnzCols);

if size(Delta_nnz,1) ~= dd - 8*d+16 || size(Delta_nnz,2) ~= 2*dd-12*d+16% = 2*(d^2 - 2*d)-(8*d-16)
    error('misunderstood size of Delta_nnz')
end

if size(Delta_nnz,2) - size(Delta_nnz,1) ~= d*(d-4)
    error('misunderstood diff size of Delta_nnz')
end

%% Delta_nnz is full row rank and misses d*(d-4) rows, hence we need d(d-4) independent vectors
% to form the basis of the null space of Delta_nnz
if size(Delta_nnz,1) ~= rank(Delta_nnz);
    error('misunderstood rank of Delta_nnz')
end

tv = createFiniteDiff2(d-2,d-2);
Delta_nnz_expected = full([kron(eye(d-4),tv)    kron(tv,eye(d-4))]);
if norm(Delta_nnz_expected - Delta_nnz) > 1e-4
    error('misunderstood Delta_nnz_expected')
end

rDelta = size(Delta_nnz,1);
cDelta_v = size(Delta_nnz,2)/2;
cDelta_h = size(Delta_nnz,2)/2;
if size(kron(eye(d-4),tv),2) ~= cDelta_v || size(kron(tv,eye(d-4)),2) ~= cDelta_h
    error('misunderstood cDelta_v or cDelta_h')
end

%% Create null space
v1_v = ones(1,d-2);
v2_v = [1:d-2]; 

% 4(d-4) vectors
V1_v = [kron(eye(d-4),v1_v) zeros(d-4, cDelta_h)];
V2_v = [kron(eye(d-4),v2_v) zeros(d-4, cDelta_h)];
V1_h = [zeros(d-4, cDelta_v) kron(v1_v,eye(d-4)) ];
V2_h = [zeros(d-4, cDelta_v) kron(v2_v,eye(d-4)) ];

% number of vectors we need in total: d*(d-4), i.e, we need other (d-4)^2 vectors (9 in a 7x7 example) 

% basis of the null space 
Nmat = [V1_v; V2_v; V1_h; V2_h; ];
% [1:cDelta_v+cDelta_h] does not add to the null space

if norm(Delta_nnz * Nmat') > 1e-5 || rank(Nmat) ~= size(Nmat,1)
    error('misunderstood null space')
end

format short
MN = [Delta_nnz; Nmat];
N2 = null(MN,'r')';
if norm(MN * N2') > 1e-5 || rank(N2) ~= size(N2,1)
    error('misunderstood null space2')
end
