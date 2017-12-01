function x = FastLDCTivSynthesis(coef,n,bellname,w,pars3)
% FastLDCTivSynthesis -- Synthesize signal from local DCT iv coefficients (orthogonal fixed folding)
%  Usage
%    sig = FastLDCTivSynthesis(ldct,bellname,w)
%  Inputs
%    coef       local DCT iv coefficients
%    w		width of window
%    bell       name of bell to use, defaults to 'Sine'
%  Outputs
%    x          signal whose orthonormal local DCT iv coeff's are ldct
%
%  See Also
%   FastLDCTivAnalysis, CPAnalysis, FCPSynthesis, fold, unfold, dct_iv, packet
%

	[n,J] = dyadlength(coef);
	d = floor(log2(n/w));
	
%
% Create Bell
%
	if nargin < 4 | bellname==0,
	  bellname = 'Sine';
	end
	m = n / 2^d /2;
	[bp,bm] = MakeONBell(bellname,m);
	
	nbox = floor(n/w);
	%lfign=0.25;
    	lfign=-1;
    	for boxcnt=0:nbox-1
            coef(boxcnt*w+1:boxcnt*w+1+floor(w*lfign)) = 0;
    	end
%
%
%
		x = zeros(1,n);
		for b=0:(2^d-1),
			   c = coef(packet(d,b,n));
			   y = dct_iv(c);
			   [xc,xl,xr] = unfold(y,bp,bm);
			   x(packet(d,b,n)) = x(packet(d,b,n)) + xc;
			   if b>0,
				   x(packet(d,b-1,n)) = x(packet(d,b-1,n)) + xl;
			   else
			       x(packet(d,0,n))   = x(packet(d,0,n)) + edgeunfold('left',xc,bp,bm);
			   end
			   if b < 2^d-1,
				   x(packet(d,b+1,n)) = x(packet(d,b+1,n)) + xr;
			   else         
			       x(packet(d,b,n))   = x(packet(d,b,n)) + edgeunfold('right',xc,bp,bm);
			   end
		 end
		
		x = x(:);

