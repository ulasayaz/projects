function C = fdct3d_forward_mod(X)
[m,n,p] = size(X);

wmaxlev=max(0,floor(log2(max(size(X))))); % wmaxlev method of MultiResTrafo2
J=floor(log2(size(X))); %size_dyadic
lev=2;
wcoarsestlev=min(J)-lev;

% nbscales=wmaxlev-wcoarsestlev+1; 
nbscales = floor(log2(min([m,n,p])))-2; %default
%nbscales=6;
nbdstz_coarse = 8;
allcurvelets = 1;
C = fdct3d_forward_mex(m,n,p,nbscales,nbdstz_coarse,allcurvelets,X);



