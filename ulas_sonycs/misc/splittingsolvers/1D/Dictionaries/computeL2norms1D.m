function E=computeL2norms1D(n,dict,pars1,pars2,pars3)
% computeL2norms -- compute (approximate) L2 norms of the atoms:
%
%  Usage:
%	E=computeL2norms(n,dict,pars1,pars2,pars3,method)
%  Inputs:
% 	n		the image size (n x n matrix)
%	dict		name of the dictionary
%	par1,par2,par3	the parameters of the dictionary
%
%	Use 'help dictionary' for dictionary objects: dict,par1,par2,par3
%  Outputs:
%	E		the L2 norms, a vector of the same size as the coeff vector provided by FastA
%  See Also:
%	FastA, FastS, MakeList
%

NumberOfDicts = LengthList(dict);
J  = nextpow2(n);
E = [];

for i = 1:NumberOfDicts,
	NAME = NthList(dict, i);
        PAR1 = NthList(pars1, i);
        PAR2 = NthList(pars2, i);
        PAR3 = NthList(pars3, i);

	if strcmp(NAME, 'DIRAC'),
		m0 = n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'LDCT'),
		m0 = n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'LDST'),
		m0 = n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'LDCTiv'),
		m0 = n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'RealFourier'),
		m0 = 2*n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'PO'),
		m0 = n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'UDWT'),
		L0 = (J-PAR1) + 1;
		m0 = L0*n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'UDWTTIGHT'),
		L0 = (J-PAR1) + 1;
		E0 = zeros(n,L0);
  		for j=1:J-PAR1 
		  E0(:,j) = 2^(-j/2);
  		end
 		E0(:,end) = 2^(-(J-PAR1)/2);
	end

	E = [E;E0(:)];
	clear E0 tmp
end


