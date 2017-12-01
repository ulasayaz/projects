function x = FastUDWTTIGHTSynthesis(wcb, n, scale, qmf, scign)
% FastUDWTTIGHTSynthesis -- Synthesis operator for Undecimated DWT (TI)
%			The transform corresponds to a tight frame with constant 1, but
%			the atoms are not normalized to a unit l_2 norm but to 2^(-j/2).
%			Tge synthesis is then the adjoint of the analysis operator. 
%  Usage:
%	x = FastUDWTTIGHTSynthesis(wcb, n, scale, qmf, scign)
%  Inputs:
%	wcb	the coefs, a structure array
%	scale	the coarsest decomposition scale
%	qmf	the quadrature mirror filter
%	scign   the number of detail scales to be ignored in the reconstruction
%  Outputs:
%	x	the synthesized signal, a column vector
%  See Also:
%	FastUDWTTIGHTAnalysis, mrdwttight, mirdwttight
%

J = nextpow2(n);
if isstr(scale) scale=str2num(scale); end
D = J-scale;
wcb = reshape(wcb,n,D+1);
x = mirdwttight(wcb(:,end),wcb(:,1:D),qmf,D);
x = x(:);
