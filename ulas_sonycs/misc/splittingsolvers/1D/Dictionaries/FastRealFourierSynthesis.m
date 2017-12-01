function x = FastRealFourierSynthesis(c,n,w,lfign,pars3)
% FastRealFourierSynthesis -- 1D signal synthesis from Real Fourier coefficients
%  Usage
%    x = FastRealFourierSynthesis(c,n,w)
%  Inputs
%    c    	Real Fourier coefficients
%    w		window width
%    lfign	low-pass band to ignore
%  Outputs
%    x        	1-d signal:  length(x)=2^J
%  Description
%    c contains coefficients of the Real Fourier Decomposition.
% See Also
%   FastRealFourierAnalysis, fft, ifft
%	

    d = floor(n/w);
    
    c = reshape(c,n,2);
    x = zeros(n,1); 
    
    if lfign,
      for p=0:d-1
    	c(p*w+1:p*w+1+floor(w*lfign),:) = 0;
      end
    end
	
    for p=0:d-1
     xft = sqrt(2)*(c(p*w+1:(p+1)*w,1) + i*c(p*w+1:(p+1)*w,2));
     x(p*w+1:(p+1)*w) = sqrt(w)*real(ifft(xft(:)));
    end
    
