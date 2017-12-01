function [m, L] = SizeOfDict2(n, NameOfDict, par1, par2, par3)
% SizeOfDict -- the size of a (merged) dictionary
%  Usage:
%	[m, L] = SizeOfDict2(n, NameOfDict, par1, par2, par3)
%  Inputs:
%	n		the image size (n x n)
%	NameOfDict	name of the dictionary
%	par1,par2,par3	the parameters of the dictionary
%
%	Use 'help dictionary' for dictionary objects: NameOfDict,par1,par2,par3
%  Outputs:
%	m		# of basis functions
%	L		# redundancy = m/n
% 

NumberOfDicts = LengthList(NameOfDict);
m = 0;
L = 0;
for i = 1:NumberOfDicts,
	NAME = NthList(NameOfDict, i);
        PAR1 = NthList(par1, i);
        PAR2 = NthList(par2, i);
        PAR3 = NthList(par3, i);

	if strcmp(NAME, 'DIRAC2'),
		L0 = 1;
		m0 = n*n;
	elseif strcmp(NAME, 'LDCT2'),
		L0 = 1;
		m0 = n*n;
	elseif strcmp(NAME, 'LDST2'),
		L0 = 1;
		m0 = n*n;
	elseif strcmp(NAME, 'LDCT2iv'),
		L0 = 1;
		m0 = n*n;
	elseif strcmp(NAME, 'RealFourier2'),
		L0 = 2;
		m0 = L0*n*n;
	elseif strcmp(NAME, 'WAVEATOM2'),
		L0 = 2;
		m0 = L0*n*n;
	elseif strcmp(NAME, 'PO2'),
		L0 = 1;
		m0 = n*n;
	elseif strcmp(NAME, 'UDWT2'),
		if isstr(PAR1),
		   PAR1 = str2num(PAR1);
		end
		J  = nextpow2(n);
		L0 = 3*(J-PAR1) + 1;
		m0 = L0*n*n;
	elseif strcmp(NAME, 'UDWTTIGHT2'),
		if isstr(PAR1),
		   PAR1 = str2num(PAR1);
		end
		J  = nextpow2(n);
		L0 = 3*(J-PAR1) + 1;
		m0 = L0*n*n;
	elseif strcmp(NAME, 'CURVWRAP2'),
		if isstr(PAR1),
		   PAR1 = str2num(PAR1);
		end
		[c,m0] = fdct_wrapping_range(n,nextpow2(n)-PAR1+1);
		L0 = m0/(n*n);
	end

	m = m + m0; L = L + L0;
end
