function wa = FastWAVEATOM2Analysis(img,pat,pars2,pars3)
% FastWAVEATOM2Analysis -- 2D Wave Atom transform (frame with redundancy 2, symmetric version)
%  Usage
%    wa = FastWAVEATOM2Analysis(img,pat)
%  Inputs
%    img	n by n image, n = 2^J 
%    pat:  	type of frequency partition which satsifies parabolic scaling relationship equal to 'p' or 'q'
%  Outputs
%    wa    	Wave Atom coefficients,
%  Description
%    The wa contains coefficients of the Wave Atom Decomposition.
% See Also
%   FastWAVEATOM2Synthesis
%
	
	
	wa = sqrt(2)*real(fatfm2sym(img,pat));
	
	wa = wa(:);


