function imr = FastWAVEATOM2Synthesis(wa,n,pat,pars2,pars3)
% FastWAVEATOM2Synthesis -- (pseudo-)inverse 2D Wave Atom transform (frame with redundancy 2, symmetric version)
%  Usage
%    wa = FastWAVEATOM2Analysis(img,pat)
%  Inputs
%    wa    	Wave Atom coefficients, a structure array (atoms are normalized to unit l_2 norm)
%    pat:  	type of frequency partition which satisfies parabolic scaling relationship equal to 'p' or 'q'
%  Outputs
%    imr	n by n image, n = 2^J 
%  Description
%    imr contains reconstruction from coefficients of the Wave Atom Decomposition.
% See Also
%   FastWAVEATOM2Analysis
%

	if length(wa)~=2*n*n
		disp('Improper Wave Atom array structure size.');
		return;
	end
	
	wa = reshape(wa,2*n,n);
	
	wa(1:2,1:2) = 0;
	wa([1:2]+n,[1:2]) = 0;
%
	imr = sqrt(2)*real(iatfm2sym(wa,pat));

