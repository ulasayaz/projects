function E=computeL2norms2D(n,dict,pars1,pars2,pars3)
% computeL2norms -- compute (approximate) L2 norms of the atoms:
%
%  Usage:
%	E=computeL2norms(n,h,dict,pars1,pars2,pars3,method)
%  Inputs:
% 	n		the image size (n x n matrix)
%	dict		name of the dictionary
%	par1,par2,par3	the parameters of the dictionary
%	method		estimation method (Dirac or AWGN)
%
%	Use 'help dictionary' for dictionary objects: dict,par1,par2,par3
%  Outputs:
%	E		the L2 norms, a vector of the same size as the coeff vector provided by FastA2
%  See Also:
%	FastA, FastS, FastS2, MakeList
%

NumberOfDicts = LengthList(dict);
J  = nextpow2(n);
E = [];

for i = 1:NumberOfDicts,
	NAME = NthList(dict, i);
        PAR1 = NthList(pars1, i);
        PAR2 = NthList(pars2, i);
        PAR3 = NthList(pars3, i);

	if strcmp(NAME, 'DIRAC2'),
		m0 = n*n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'LDCT2'),
		m0 = n*n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'LDST2'),
		m0 = n*n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'LDCT2iv'),
		m0 = n*n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'RealFourier2'),
		L0 = 2;
		m0 = L0*n*n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'WAVEATOM2'),
		L0 = 2;
		m0 = L0*n*n;
		E0 = ones(m0,1);
	elseif strcmp(NAME, 'PO2'),
		m0 = n*n;
  		E0 = ones(m0,1);
	elseif strcmp(NAME, 'UDWT2'),
		if isstr(PAR1),
		   PAR1 = str2num(PAR1);
		end
		L0 = 3*(J-PAR1) + 1;
		m0 = L0*n*n;
  		E0 = ones(m0,1);
	elseif strcmp(NAME, 'UDWTTIGHT2'),
		if isstr(PAR1),
		   PAR1 = str2num(PAR1);
		end
		L0 = 3*(J-PAR1) + 1;
		m0 = L0*n*n;
		E0 = zeros(n,m0);
  		for j=1:J-PAR1 
		  E0(:,[(j-1)*3*n+1:j*3*n]) = 2^(-j);
  		end
 		E0(:,end-n+1:end) = 2^(-(J-PAR1));
	elseif strcmp(NAME, 'CURVWRAP2'),
		if isstr(PAR1),
		   PAR1 = str2num(PAR1);
		end
		[c,m0] = fdct_wrapping_range(n,J-PAR1+1);
		L0 = m0/(n*n);
		E0 = ones(m0,1)/sqrt(L0);
	end

	E = [E;E0(:)];
	clear E0 tmp
end
