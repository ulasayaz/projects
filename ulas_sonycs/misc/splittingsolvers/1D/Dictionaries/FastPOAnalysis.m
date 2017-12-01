function c = FastPOAnalysis(x, par1, par2, par3)
% FastPOSynthesis -- Fast Analysis Operator for Periodized-Orthognal
%			Wavelets Dictionary
%  Usage:
%	c = FastPOAnalysis(x, L, qmf)
%  Inputs
%    x    1-d signal; length(x) = 2^J
%    L    Coarsest Level of V_0;  L << J
%    qmf  quadrature mirror filter (orthonormal)
%  Outputs
%    c    1-d wavelet transform of x.
%
%  Description
%    1. qmf filter may be obtained from MakeONFilter   
%    2. usually, length(qmf) < 2^(L+1)
% 
%  See Also
%    FastPOSynthesis, FWT_PO, IWT_PO, MakeONFilter
%

n = length(x);
c = FWT_PO(x, par1, par2);
c = c(:);
