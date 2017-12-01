function x = FastUDWTSynthesis(wcb, n, scale, qmf, scign)
% FastUDWTSynthesis -- Synthesis operator for Undecimated DWT (TI)
%  Usage:
%	x = FastUDWTSynthesis(wcb, n, scale, qmf, scign)
%  Inputs:
%	wcb	the coefs, a structure array
%	scale	the coarsest decomposition scale
%	qmf	the quadrature mirror filter
%	scign   the number of detail scales to be ignored in the reconstruction
%  Outputs:
%	x	the synthesized signal, a column vector
%  See Also:
%	FastUDWTAnalysis, mrdwt, mirdwt
%

J = nextpow2(n);
if isstr(scale) scale=str2num(scale); end
scale = J-scale;
wcb = reshape(wcb,n,scale+1);
x = mirdwt(wcb(:,end),wcb(:,1:scale),qmf,scale);
x = x(:)*(scale+1);
