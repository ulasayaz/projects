function x = FastWPSynthesis(c, n, D, qmf, par3)
% FastWPSynthesis -- Synthesis operator for Wavelet Packet dictionary
%  Usage:
%	x = FastWPSynthesis(c, D, qmf)
%  Inputs:
%	c	the coefs, a column vector
%	D	the depth of wavelet packet
%	qmf	the quadrature mirror filter
%  Outputs:
%	x	the synthesized signal, a column vector
%  See Also:
%	FastWPAnalysis, WPAnalysis, FWPSynthesis
%

m = length(c);
L = D + 1;
pkt = reshape(c, m/L, L);
x = FWPSynthesis(pkt, qmf);
x = x(:);
