function c = FastUDWTTIGHT2Analysis(img, D, qmf, par3)
% FastUDWTTIGHT2Analysis -- Analysis operator for Undecimated DWT (TI)
%			   The transform corresponds to a tight frame with constant 1, but
%			   the atoms are not normalized to a unit l_2 norm but to 2^(-j/2). 
%  Usage:
%	c = FastUDWTTIGHT2Analysis(x, D, qmf)
%  Inputs:
%	img	the image, a nxn matrix
%	D	the depth of analysis
%	qmf	the quadrature mirror filter
%  Outputs:
%	c	the coefs, a column vector
%  See Also:
%	mrdwt, mirdwt, FastUDWTTIGHT2Synthesis
%

if isstr(D),
   D = str2num(D);
end

[n,J] = quadlength(img);
if isstr(D) D=str2num(D); end
scale = J-D;
[ll,wc] = mrdwttight(img, qmf, scale);
c = [wc(:);ll(:)];


