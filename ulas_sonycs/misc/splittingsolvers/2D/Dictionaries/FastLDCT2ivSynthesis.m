function img = FastLDCT2ivSynthesis(coef,n,bellname,w,pars3)
% FastLDCT2Synthesis -- 2D Local inverse DCT iv transform
%  Usage
%    img = FastLDCT2Synthesis(ldct,w) 
%  Inputs
%    coef   	2D Local DCT structure array
%    w        	width of window
%    bellname name of bell to use, defaults to 'Sine'
%  Outputs
%    img	    2D reconstructed n by n image
%  Description
%    The matrix img contains image reconstructed from the Local DCT Decomposition.
% See Also
%   FastLDCT2Analysis, idct2
%

	if nargin < 4 | bellname==0,
	  bellname = 'Sine';
	end
	
	J = nextpow2(n);
	coef = reshape(coef,n,n);
	

	img = zeros(n,n);
%
	nbox = floor(n/w);
    lfign=0.25;
    %lfign=-1;
    for boxcnt1=0:nbox-1
        for boxcnt2=0:nbox-1
            coef(boxcnt1*w+1:boxcnt1*w+1+floor(w*lfign),boxcnt2*w+1:boxcnt2*w+1+floor(w*lfign)) = 0;
        end
    end

    for ncol=1:n
           ldct = coef(:,ncol);
           x = FastLDCTivSynthesis(ldct,n,bellname,w);
           img(:,ncol) = x;
       end
       
    for nrow=1:n
           ldct = img(nrow,:);
           x = FastLDCTivSynthesis(ldct,n,bellname,w);
           img(nrow,:) = x';
    end


