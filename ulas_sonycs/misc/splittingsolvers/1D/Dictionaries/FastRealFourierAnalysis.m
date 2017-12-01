function c = FastRealFourierAnalysis(x,w,pars2,pars3)
% FastRealFourierAnalysis -- Real Fourier transform for 1D signals
%  Usage
%    c = FastRealFourierAnalysis(x,w) 
%  Inputs
%    x        	1-d signal:  length(x)=2^J
%    w		window width
%  Outputs
%    c    	Real Fourier coefficients
%  Description
%    c contains coefficients of the Real Fourier Decomposition.
% See Also
%   FastRealFourierSynthesis, fft, ifft
%		

    
	[n,J] = dyadlength(x);
	
	d = floor(n/w);
	
	c = zeros(n,2);


	for p=0:d-1	
	  xft    = fft(x(p*w+1:(p+1)*w))*sqrt(2/w);
	  c(p*w+1:(p+1)*w,1) = real(xft(:));
	  c(p*w+1:(p+1)*w,2) = imag(xft(:));
	end
	
	c = c(:);
