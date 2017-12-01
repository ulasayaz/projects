function img = FastUDWT2Synthesis(c, n, D, qmf, par3)
% FastUDWT2Synthesis -- Synthesis operator for UDWT
%  Usage:
%	img = FastUDWT2Synthesis(c, D, qmf)
%  Inputs:
%	c	the coefs, a column vector
%	D	the depth of analysis
%	qmf	the quadrature mirror filter
%  Outputs:
%	img	the synthesized image, a nxn image
%  See Also:
%	mrdwt, mirdwt, FastUDWT2TAnalysis

if isstr(D),
   D = str2num(D);
end

J = nextpow2(n);
if isstr(D) D=str2num(D); end
scale = J-D;
ll = reshape(c(end-n*n+1:end),n,n);
wc = reshape(c(1:end-n*n),n,3*scale*n);
img = mirdwt(ll, wc, qmf, scale)*(3*scale+1);
