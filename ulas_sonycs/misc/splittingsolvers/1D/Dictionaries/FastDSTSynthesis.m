function x = FastDSTSynthesis(c, n, fineness, par2, par3)
% FastDstAnalysis -- Synthesis operator for (overcomplete) Dst dictionary
%  Usage 
%    	x = FastDSTSynthesis(c, n, finess)
%  Inputs:
%	c	the coef
%   fineness	how many times finer the dictionary is compared to the 
%		standard one, eg: 1, 2, 4, 8 ...
%  Output:
%	x	the synthesized signal
%  See Also:
%	FastDSTAnalysis
%


m = length(c); n = m/fineness;
f = (1:(m-1))';
const = [n; .5 * (n +  (1 - cos(2*pi*f/fineness)) ./ (2*sin(pi*f/m)))];
const = const .^ .5;

n2 = 2 * n;
c = c ./ const;
z = zeros(4*m, 1);
z(1:m) = c;
y = fft(z);
x = imag(y(2:2:n2));
x = x(:);

