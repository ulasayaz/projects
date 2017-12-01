function xls = LSsolution(A, y, x, p, activeSet)
% Return the LS solution from the support of x, i.e.
%		min || y -A_Ix_I ||_2 .
%
% Input
%	A           nxp matrix or implicit operator
%	y           n-vector
%	x	    input vector whose support is I
%  	p           length of x
%	activeSet   array containing the support of the solution, default: estimated from x
% Outputs
%	 xls        LS solution

global dict pars1 pars2 pars3

fastop = (ischar(A) || isa(A, 'function_handle'));

n = length(y);

if ~exist('p'), 	 p = size(A,2); 	end
if ~exist('activeSet') | length(activeSet)==p,  activeSet = find(abs(x)  > 3*max(MAD(x),min(abs(x(find(x))))-eps)); end
% lsqrms params
damp = 0;
atol = 1e-6;
btol = 1e-6;
conlim = 1e+10;
itnlim = 10;
show = 0;

xls = zeros(p,1);
if fastop
    N = length(activeSet);
    [xls(activeSet), istop, itn] = lsqrms(n, N, A, activeSet, p, y, damp, atol, btol, conlim, itnlim, show);
else
    AI = A(:,activeSet);
    xls(activeSet) = AI \ y;
end

