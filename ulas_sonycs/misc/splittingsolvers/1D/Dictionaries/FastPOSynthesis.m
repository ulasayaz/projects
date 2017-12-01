function x = FastPOSynthesis(c, n, par1, par2, par3)
% FastPOSynthesis -- Fast Synthesis Operator for Periodized-Orthognal
%			Wavelets Dictionary
%  Usage:
%	x = FastPOSynthesis(c, L, qmf)
%  Inputs
%    c    1-d wavelet transform of x.
%    L    Coarsest Level of V_0;  L << J
%    qmf  quadrature mirror filter (orthonormal)
%  Outputs
%    x    1-d signal; length(x) = 2^J
%
%  Description
%    1. qmf filter may be obtained from MakeONFilter   
%    2. usually, length(qmf) < 2^(L+1)
% 
%  See Also
%    FastPOAnalysis, FWT_PO, IWT_PO, MakeONFilter
%

n = length(c);
x = IWT_PO(c, par1, par2);
x = x(:);
