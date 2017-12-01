function X = FastCURVWRAP2Synthesis(C,n,L,pars2,pars3)
% ifdct_wrapping: Returns the inverse curvelet transform
%		 implemented using the wrapped IFCT (see Curvelab at curvelets.org)
%  Usage:
%    X =  ifdct_wrapping(C,isreal,L);
%  Inputs:
%    C     Vector which contains the curvelet coefficients at
%    dyadic scales. 
%    L     Coarsest scale; L <<J
%  Outputs:
%    X	   n by n image, n = 2^J 
%
%  Description
%    This essentially the wrapping version of CurveLab, using a decimated
%    rectangular grid. The transform is a numerical
%    isometry and can be inverted by its adjoint (up to a constant).
%    In this verison, we have wavelet coefficients at the coarsest
%    and finest scales.
%    See curvelets.org for more details on CurveLab.
%  See Also
%    ifdct_wrapping, FastCURVWRAPAnalysis


       J = nextpow2(n);
       if isstr(L) L=str2double(L); end
       IsImageReal = 1;
       C = C*sqrt(length(C)/(n*n));
                
%
%      Create cell array for compatibility with CurveLab inversion
%
       c = fdct_wrapping_range(n, J-L+1);
       CW{1}{2} = [n n]; 
       
       ind = 0;
       m = c{1}.sw(1);
       n = c{1}.sw(2);
       CW{1}{1} = reshape(C(ind+1:ind+m*n),m,n);
       ind = m*n;
       
       for j = 2:length(c),
       for w=1:c{j}.nbangles,
	     % Copy the coefficients back
	     m = c{j}.sw(w,1);
	     n = c{j}.sw(w,2);
	     CW{j}{w} = reshape(C(ind+1:ind+m*n),m,n);
	     ind = ind + m*n;
       end
       end
            
       % Apply the inverse transform
       X = real(ifdct_wrapping(CW, IsImageReal, J-L+1));
       
      

      
	      
	      
