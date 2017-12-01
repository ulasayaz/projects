function c = FastUDWT2Analysis(img, D, qmf, par3)
% FastUDWT2Analysis -- Analysis operator for UDWT
%  Usage:
%	c = FastUDWT2Analysis(x, D, qmf)
%  Inputs:
%	img	the image, a nxn matrix
%	D	the depth of analysis
%	qmf	the quadrature mirror filter
%  Outputs:
%	c	the coefs, a column vector
%  See Also:
%	mrdwt, mirdwt, FastUDWT2Synthesis
%

if isstr(D),
   D = str2num(D);
end

[n,J] = quadlength(img);
if isstr(D) D=str2num(D); end
scale = J-D;
[ll,wc] = mrdwt(img, qmf, scale);
c = [wc(:);ll(:)];

