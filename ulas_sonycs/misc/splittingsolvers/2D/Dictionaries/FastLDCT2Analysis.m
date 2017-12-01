function ldct = FastLDCT2Analysis(img,w,lfign,pars3)
% FastLDCTAnalysis -- 2D Local transform using a DCT dictionary
%  Usage
%    ldct = FastLDCT2Analysis(img,w) 
%  Inputs
%    img	n by n image, n = 2^J 
%    w        	width of window
%    lfign	ignore the local zero frequency component in each square
%  Outputs
%    ldct    	Local DCT coefficients
%  Description
%    The ldct contains coefficients of the Local DCT Decomposition.
% See Also
%   FastLDCT2Synthesis
%

	[n,J] = quadlength(img);
	
	d = floor(n/w);

	ldct = zeros(n,n);
%
	
	for p1=0:d-1
	 for p2=0:d-1
	     imgp = img(p1*w+1:(p1+1)*w,p2*w+1:(p2+1)*w);
	     c = dct2(imgp);	   				% DCT analysis
	     ldct(p1*w+1:(p1+1)*w,p2*w+1:(p2+1)*w) = c;% store
	     if lfign
	     	ldct(p1*w+1,p2*w+1) = 0;  % The zero frequency component is set to zero
	        ldct(p1+1,p2+1)     = c(1,1); % store the zero frequency component
	     end
	   end
	end

