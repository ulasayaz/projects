function coef = FastLDCT2ivAnalysis(img,bellname,w,pars3)
% FastLDCT2ivAnalysis -- Analyze image into 2-d cosine packet coefficients at a given depth (window width)
%  Usage
%    coef = FastLDCT2ivAnalysis(img,w)
%  Inputs
%    img      2-d image to be transformed into basis
%    w        width of window
%    bellname name of bell to use, defaults to 'Sine'
%  Outputs
%    coef     2-d Local DCT iv coeffs
%
%  Description
%    Once a cosine packet basis has been selected (at a given depth), 
%    this function may be used to expand a given
%    image in that basis.
%

   	if nargin < 3 | bellname==0,
	  bellname = 'Sine';
	end
    
	[n,J] = quadlength(img);
	
	d = floor(log2(n/w));

%
% CP image at depth d
%
%

       coef = zeros(n,n);
       
       for nrow=1:n
           ldct = FastLDCTivAnalysis(img(nrow,:),bellname,w);
           coef(nrow,:) = ldct(:)';
       end
       
       for ncol=1:n
           ldct = FastLDCTivAnalysis(coef(:,ncol),bellname,w);
           coef(:,ncol) = ldct(:);
       end
       
       coef = coef(:);
       
