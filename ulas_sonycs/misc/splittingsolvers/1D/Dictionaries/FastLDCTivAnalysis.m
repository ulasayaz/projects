function coef = FastLDCTivAnalysis(x,bellname,w,par3)
% FastLDCTivAnalysis -- Local DCT iv transform (orthogonal fixed folding)
%  Usage
%    ldct = FastLDCTivAnalysis(x,bell,w) 
%  Inputs
%    x        1-d signal:  length(x)=2^J
%    w        width of window
%    bell     name of bell to use, defaults to 'Sine'
%  Outputs
%    coef     1-d Local DCT iv coefficients
%  Description
%    The vector coef contains coefficients of the Local DCT Decomposition.
% See Also
%   FastLDCTivSynthesis, CPAnalysis, FCPSynthesis, fold, unfold, dct_iv, packet
%

   	if nargin < 3 | bellname==0,
	  bellname = 'Sine';
	end
    
	[n,J] = dyadlength(x);
	
	d = floor(log2(n/w));

%
% CP image at depth d
%
%

%
% taper window
%
	m = n / 2^d /2;
	[bp,bm] = MakeONBell(bellname,m);
%
% packet table
%
	n  = length(x);
	x  = ShapeAsRow(x);
	coef = zeros(n,1);
%
	   nbox = 2^d;
	   for b=0:(nbox-1)
		   if(b == 0) ,                             % gather packet and
			   xc = x(packet(d,b,n));           % left, right neighbors
			   xl = edgefold('left',xc,bp,bm);  % taking care of edge effects
		   else
			   xl = xc;
			   xc = xr;          
		   end
		   if (b+1 < nbox)
			   xr = x(packet(d,b+1,n));
		   else
			   xr = edgefold('right',xc,bp,bm);
		   end
		   y = fold(xc,xl,xr,bp,bm);    % folding projection
		   c = dct_iv(y);               % DCT-IV
		   coef(packet(d,b,n)) = c';  % store
	   end
       
