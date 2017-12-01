function x = FastPOSpinSynthesis(c, n, par1, par2, par3)
% FastPOSynthesis -- Fast Synthesis Operator for Periodized-Orthognal with cycle spining
%			Wavelets Dictionary
%  Usage:
%	x = FastPOSSpinynthesis(c, L, qmf)
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

x = IWT_PO(c(1:n), par1, par2);
c = c(n+1:end);
for i=1:par3
  x = x + circshift(IWT_PO(c(1:n), par1, par2),-i);
  c = c(n+1:end);
  x = x + circshift(IWT_PO(c(1:n), par1, par2),i);
  c = c(n+1:end);
end

