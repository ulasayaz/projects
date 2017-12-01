function c = FastCPAnalysis(x, D, par2, par3)
% FastCPAnalysis -- Analysis operator for Wavelet Packet dictionary
%  Usage:
%	c = FastCPAnalysis(x, D)
%  Inputs:
%	x	the signal, a column vector
%	D	the depth of wavelet packet
%  Outputs:
%	c	the coefs, a column vector
%  See Also:
%	FastCPSynthesis, CPAnalysis, FCPSynthesis

c = CPAnalysis(x, D, 'Sine');
c = c(:);
