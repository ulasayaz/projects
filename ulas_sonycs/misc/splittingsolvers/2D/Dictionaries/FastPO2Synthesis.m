function img = FastPO2Synthesis(c, n, scale, qmf, par3)
% FastPO2Synthesis -- Fast Synthesis Operator for Periodized-Orthognal
%			Wavelets Dictionary
%  Usage:
%	img = FastPOSynthesis(c, n, L, qmf)
%  Inputs
%    c    1-d wavelet transform of x.
%    L    Coarsest Level of V_0;  L << J
%    qmf  quadrature mirror filter (orthonormal)
%  Outputs
%    img    2-d image; length(x) = 2^J
%
%  Description
%    1. qmf filter may be obtained from MakeONFilter   
%    2. usually, length(qmf) < 2^(L+1)
% 
%  See Also
%    FastPO2Analysis, FWT2_PO, IWT2_PO, MakeONFilter
%

if isstr(scale),
   scale = str2num(scale);
end

c = reshape(c, n, n);
img = IWT2_PO(c, scale, qmf);
