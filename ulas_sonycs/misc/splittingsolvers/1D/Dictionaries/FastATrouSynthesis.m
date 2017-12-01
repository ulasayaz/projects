function x = FastATrouSynthesis(c, n, D, par2, par3)
% FastATrouSynthesis -- Synthesis operator of a trous algotithm
%  Usage:
%	x = FastATrouSynthesis(c, D)
%  Inputs:
%	c	the coefs, a column vector
%	D	the depth of wavelet packet
%  Outputs:
%	x	the synthesized signal, a column vector
%  See Also:
%	FastATrouAnalysis

J = nextpow2(n);
L = (J-D) + 1;
pkt = reshape(c, n, L);
x = IWT_ATrou(pkt, J-D);
x = x(:)*L;
