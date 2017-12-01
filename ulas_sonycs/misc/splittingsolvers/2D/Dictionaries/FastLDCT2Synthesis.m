function img = FastLDCT2Synthesis(ldct,n,w,lfign,pars3)
% FastLDCT2Synthesis -- Local inverse DCT transform
%  Usage
%    img = FastLDCT2Synthesis(ldct,w) 
%  Inputs
%    ldct   	2D Local DCT
%    w        	width of window
%    lfign	restore the local zero frequency component in each square
%  Outputs
%    img	2D reconstructed n by n image
%  Description
%    The matrix img contains image reconstructed from the Local DCT Decomposition.
% See Also
%   FastLDCT2Analysis, idct2
%

	
	ldct = reshape(ldct,n,n);
	
	d = floor(n/w);

	img = zeros(n,n);
%
	
	for p1=0:d-1
	 for p2=0:d-1
	     ldctp = ldct(p1*w+1:(p1+1)*w,p2*w+1:(p2+1)*w);
	     if lfign
	     	ldctp(1,1) = ldct(p1+1,p2+1);
	     end
	     c = idct2(ldctp);	   				% DCT analysis
	     img(p1*w+1:(p1+1)*w,p2*w+1:(p2+1)*w) = c;		% store
	   end
	end


    
