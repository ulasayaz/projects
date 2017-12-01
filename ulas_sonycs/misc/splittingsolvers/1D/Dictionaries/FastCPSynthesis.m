function x = FastCPSynthesis(c, n, D, par2, par3)
% FastCPSynthesis -- Synthesis operator for Wavelet Packet dictionary
%  Usage:
%	x = FastCPSynthesis(c, D)
%  Inputs:
%	c	the coefs, a column vector
%	D	the depth of wavelet packet
%  Outputs:
%	x	the synthesized signal, a column vector
%  See Also:
%	FastCPAnalysis, CPAnalysis, FCPSynthesis
%

m = length(c);
L = D + 1;
temp = reshape(c, m/L, L);
x = FCPSynthesis(temp,'Sine');
x = x(:);
