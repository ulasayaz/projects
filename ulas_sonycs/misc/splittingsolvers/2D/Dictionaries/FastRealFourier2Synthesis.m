function x = FastRealFourier2Synthesis(c,n,w,pars2,pars3)
% FastRealFourier2Synthesis -- 2D synthesis from Real Fourier coefficients
%  Usage
%    img = FastRealFourierSynthesis(c,n,w)
%  Inputs
%    c    	Real Fourier coefficients
%    w		window width
%  Outputs
%    img	2-d image:  length(x)=2^J
%  Description
%    c contains coefficients of the Real Fourier Decomposition.
% See Also
%   FastRealFourier2Analysis, fft2, ifft2
%	

    d = floor(n/w);
    
    c = reshape(c,n,2*n);
    x = zeros(n,n); 
    
    for p1=0:d-1
     for p2=0:d-1
     xft = sqrt(2)*(c(p1*w+1:(p1+1)*w,p2*w+1:(p2+1)*w) + i*c(p1*w+1:(p1+1)*w,p2*w+1+n:(p2+1)*w+n));
     x(p1*w+1:(p1+1)*w,p2*w+1:(p2+1)*w) = w*real(ifft2(xft));
     end
    end
    
