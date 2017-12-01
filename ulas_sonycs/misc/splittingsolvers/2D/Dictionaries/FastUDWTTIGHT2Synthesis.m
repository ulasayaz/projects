function img = FastUDWTTIGHTSynthesis(c, n, D, qmf, scign)
% FastUDWTTIGHTSynthesis -- Synthesis operator for Undecimated DWT (TI)
%			The transform corresponds to a tight frame with constant 1, but
%			the atoms are not normalized to a unit l_2 norm but to 2^(-j/2).
%			Tge synthesis is then the adjoint of the analysis operator. 
%  Usage:
%	x = FastUDWTTIGHTSynthesis(wcb, n, scale, qmf, scign)
%  Inputs:
%	c	the coefs, a column vector
%	D	the depth of analysis
%	qmf	the quadrature mirror filter
%  Outputs:
%	img	the synthesized image, a nxn image
%  See Also:
%	FastUDWTTIGHT2Analysis, mrdwttight, mirdwttight
%

if isstr(D),
   D = str2num(D);
end

J = nextpow2(n);
if isstr(D) D=str2num(D); end
scale = J-D;
ll = reshape(c(end-n*n+1:end),n,n);
wc = reshape(c(1:end-n*n),n,3*scale*n);
img = mirdwttight(ll, wc, qmf, scale);
