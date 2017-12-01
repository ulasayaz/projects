function wcb = FastUDWTTIGHTAnalysis(x, scale, qmf, par3)
% FastUDWTTIGHTAnalysis -- Analysis operator for Undecimated DWT (TI)
%			   The transform corresponds to a tight frame with constant 1, but
%			   the atoms are not normalized to a unit l_2 norm but to 2^(-j/2). 
%  Usage:
%	wcb = FastUDWTTIGHTAnalysis(x, scale, qmf)
%  Inputs:
%	x	the signal, a column vector
%	scale	the coarsest decomposition scale
%	qmf	the quadrature mirror filter
%  Outputs:
%	wc	the coefs, a structure array
%  See Also:
%	FastUDWTTIGHTSynthesis, mrdwttight, mirdwttight
%

[n,J] = dyadlength(x);
if isstr(scale) scale=str2num(scale); end
D = J-scale;
[ll,wc,L] = mrdwttight(x,qmf,D);
wcb = [wc(:);ll(:)];


