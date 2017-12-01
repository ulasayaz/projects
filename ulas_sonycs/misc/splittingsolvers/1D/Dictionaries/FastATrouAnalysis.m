function c = FastATrouAnalysis(x, D, par2, par3)
% FastATrouAnalysis -- Analysis operator of the WT by a trous algorithm
%  Usage:
%	c = FastATrouAnalysis(x, D, qmf)
%  Inputs:
%	x	the signal, a column vector
%	D	the depth of wavelet packet
%  Outputs:
%	c	the coefs, a column vector
%  See Also:
%	FastATrouSynthesis
%

[n,J] = dyadlength(x);
pkt = FWT_ATrou(x, D);
c = pkt(:);

