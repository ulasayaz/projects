function C = FastCURVWRAP2Analysis(Img,L,pars2,pars3)
% fdct_wrapping: Returns the curvelets tight frame coefficients of an image
%		 implemented using the wrapped FCT (see Curvelab at curvelets.org)
%  Usage:
%    C =  fdct_wrapping(Img,isreal,L);
%  Inputs:
%    X     n by n image, n = 2^J 
%    L     Coarsest scale; L <<J
%  Outputs:
%    C coefficients vector which contains the curvelet coefficients at
%    dyadic scales. 
%
%  See Also
%    fdct_wrapping, FastCURVWRAP2Synthesis


       [n,J] = quadlength(Img);
       if isstr(L) L=str2double(L); end
       IsImageReal = isreal(Img);
          
       CW = fdct_wrapping(Img, IsImageReal, J-L+1);
      
%
%      Create data structure
%
       % Coarsest scale approx coeffs
       [c,nall] = fdct_wrapping_range(n, J-L+1);
       C = zeros(nall,1);
       C(1:prod(c{1}.sw)) = CW{1}{1}(:);
       ind = prod(c{1}.sw);
       
       for j = 2:length(CW),
       for w=1:length(CW{j})
	     C(ind+1:ind+prod(c{j}.sw(w,:))) = CW{j}{w}(:);
	     ind = ind + prod(c{j}.sw(w,:));
       end
       end
       
       C = C*sqrt(length(C)/(n*n));
       
       
      

      
	      
	      
