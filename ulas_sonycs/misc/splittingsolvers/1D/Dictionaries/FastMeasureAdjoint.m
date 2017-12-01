function x = FastMeasureAdjoint(y,p,dict,Omega)
%
% Apply the adjoint measurement matrix implicit operator.
%
% Input
%	y           input data
%	dict        sensing matrix name
%	Omega	    nx1 vector denoting which measurements are kept
% Output
%	y	    output vector

if strcmp(upper(dict),'FOURIER')
	X = zeros(p,1);
	X(Omega) = y;
	x = sqrt(p)*real(ifft(X(:))); % Taking only the real part may accumulate numerical errors.
elseif strcmp(upper(dict),'RST')
	X = zeros(p,1);
	X(Omega) = y;
	x = Inv_RST(X(:));
elseif strcmp(upper(dict),'REALFOURIER') % Redundancy = 2, [real part; imag part] instead of complex Fourier measurements.
	X = zeros(p,1);
	X(Omega) = sqrt(2)*(y(1:end/2) + i*y(end/2+1:end));
	x = sqrt(p)*real(ifft(X(:))); % Taking only the real part may accumulate numerical errors.
elseif strcmp(upper(dict),'HADAMARD')
	X = zeros(p,1);
	X(Omega) = y;
	x = Inv_FHT(X(:));
elseif strcmp(upper(dict),'DIRAC')
	X = zeros(p,1);
	X(Omega) = y;
	x = X(:);
else
	disp('Uknown measurement matrix');
end







    
    
