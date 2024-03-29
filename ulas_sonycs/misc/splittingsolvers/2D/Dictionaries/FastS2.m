function x = FastS2(c, n, NameOfDict, par1, par2, par3)
% FastS2 -- the SYNTHESIS operator for a dictionary:
%			x = \Phi * c
%  Usage:
%	x = FastLS2(c, n, NameOfDict, par1, par2, par3)
%  Inputs:
% 	c		the coefs, a structure array
%	NameOfDict	name of the dictionary
%	par1,par2,par3	the parameters of the dictionary
%
%	Use 'help dictionary' for dictionary objects: NameOfDict,par1,par2,par3
%  Outputs:
%	x		the synthesized image, a n by n matrix
%  See Also:
%	FastA2, FastA, FastS, MakeList
%

NumberOfDicts = LengthList(NameOfDict);
if NumberOfDicts == 1,
	cmmdstr = ['Fast' NameOfDict 'Synthesis(c, n, par1, par2, par3)'];
	x = eval(cmmdstr);
else
	x = 0;
	ind = 0;
	for i = 1:NumberOfDicts,
                NAME = NthList(NameOfDict, i);
                PAR1 = NthList(par1, i);
                PAR2 = NthList(par2, i);
                PAR3 = NthList(par3, i);
		DimDomain = SizeOfDict2(n, NAME, PAR1, PAR2, PAR3);
		C = c((ind+1):(ind+DimDomain));
                x = x +  eval(['Fast' NAME, 'Synthesis(C, n, PAR1, PAR2, PAR3)']);
		ind = ind + DimDomain;
	end
end
