function c = FastRealFourier2Analysis(img,w,pars2,pars3)
% FastRealFourier2Analysis -- Real Fourier transform for 2D images
%  Usage
%    c = FastRealFourierAnalysis(img,w) 
%  Inputs
%    img	2-d image:  length(img)=2^J
%    w		window width
%  Outputs
%    c    	Real Fourier coefficients
%  Description
%    c contains coefficients of the Real Fourier Decomposition.
% See Also
%   FastRealFourier2Synthesis, fft2, ifft2
%		

    
	[n,J] = quadlength(img);
	
	d = floor(n/w);
	
	c = zeros(n,2*n);


	for p1=0:d-1
	 for p2=0:d-1	
	  xft    = fft2(img(p1*w+1:(p1+1)*w,p2*w+1:(p2+1)*w))*sqrt(2)/w;
	  c(p1*w+1:(p1+1)*w,p2*w+1:(p2+1)*w)     = real(xft);
	  c(p1*w+1:(p1+1)*w,p2*w+1+n:(p2+1)*w+n) = imag(xft);
	 end
	end
	
	c = c(:);
	
