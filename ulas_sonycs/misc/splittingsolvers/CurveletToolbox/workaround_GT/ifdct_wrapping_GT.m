function X = ifdct_wrapping_GT(C, isreal, allcurvelets, nbscales, nbangles_coarse)

% ifdct_wrapping - Inverse curvelet transform
%
% Input
%     C         Curvelet coefficients
%     isreal    Type of transform
%                   0: complex
%                   1: real
%    allcurvelets (bool) 
% Output
%     X         A double precision matrix
%
% See also ifdct_wrapping in the fdct_wrapping_matlab/ directory.

  m = C{1}{2}(1);
  n = C{1}{2}(2);
  J = nextpow2(n);
  n = 2^J;
  
  if ~exist('nbscales','var')
  	nbscales = floor(log2(min(m,n)))-3;
  end
  if ~exist('nbangles_coarse','var')
  	nbangles_coarse = 16;
  end
  if ~exist('allcurvelets','var')
    allcurvelets = true;
  end
  
  if allcurvelets
      finest=1;
  else
      finest=2;
  end

% the m.file used here (ifdct_wrapping_workaround)
% has its own internal routine for the real case.
% The following lines must be commented out (else quality reduced)
%   if(isreal)
%     C = fdct_wrapping_r2c(C);
%   end
  
  % call mex function
  %X = ifdct_wrapping_mex(m,n,nbscales, nbangles_coarse, allcurvelets, C);
  X = ifdct_wrapping_workaround(C, isreal, finest, nbscales, nbangles_coarse,m,n);
  

  
