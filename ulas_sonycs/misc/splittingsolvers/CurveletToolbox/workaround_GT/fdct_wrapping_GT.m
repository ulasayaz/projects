function C = fdct_wrapping_GT(X, isreal, allcurvelets, nbscales, nbangles_coarse)

% fdct_wrapping - Forward curvelet transform
%
% Inputs
%     X         a double precision matrix
%     isreal    Type of transform
%                   0: complex
%                   1: real
%     allcurvelets (bool) 
%     nbscales  Coarsest decomposition scale (default: floor(log2(min(m,n)))-3)
%
% Output
%     C         Curvelet coefficients
%
% See also fdct_wrapping.m in the fdct_wrapping_matlab/ directory.
  
  [m,n] = size(X);
  [n,J] = quadlength(X);
  
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
  
  %call mex function
  %C = fdct_wrapping_mex(m,n,nbscales, nbangles_coarse, allcurvelets, double(X));
 
  C = fdct_wrapping_workaround(double(X), isreal, finest, nbscales, nbangles_coarse);

% the m.file used here (fdct_wrapping_workaround)
% has its own internal routine for the real case.
% The following lines must be commented out (else quality reduced)  
%   if(isreal)
%     C = fdct_wrapping_c2r(C);
%   end

  C{1}{2} = [m n];
   
