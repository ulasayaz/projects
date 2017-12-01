function c = FastWPAnalysis(x, D, qmf, par3)
% FastWPAnalysis -- Analysis operator for Wavelet Packet dictionary
%  Usage:
%	c = FastWPAnalysis(x, D, qmf)
%  Inputs:
%	x	the signal, a column vector
%	D	the depth of wavelet packet
%	qmf	the quadrature mirror filter
%  Outputs:
%	c	the coefs, a column vector
%  See Also:
%	FastWPSynthesis, WPAnalysis, FWPSynthesis
%

c = WPAnalysis(x, D, qmf);
c = c(:);
