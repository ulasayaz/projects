function c = FastA2(x, NameOfDict, par1, par2, par3)
% FastA2 -- the ANALYSIS operator for a dictionary:
%			c = \Phi^T * x
%  Usage:
%	c = FastA2(x, NameOfDict, par1, par2, par3)
%  Inputs:
% 	x		the image, a n by n matrix
%	NameOfDict	name of the dictionary
%	par1,par2,par3	the parameters of the dictionary
%
%	Use 'help dictionary' for dictionary objects: NameOfDict,par1,par2,par3
%  Outputs:
%	c		the coefs, a structure array
%  See Also:
%	FastA, FastS, FastS2, MakeList
%

NumberOfDicts = LengthList(NameOfDict);
if NumberOfDicts == 1,
	C = eval(['Fast' NameOfDict 'Analysis(x, par1, par2, par3)']);
	c = C(:);
else
	c = [];
	for i = 1:NumberOfDicts,
		NAME = NthList(NameOfDict, i);
		PAR1 = NthList(par1, i);
		PAR2 = NthList(par2, i);
		PAR3 = NthList(par3, i);
		C = eval(['Fast' NAME, 'Analysis(x, PAR1, PAR2, PAR3)']);
		c = [c; C(:)];
	end
end
