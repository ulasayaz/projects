function y = FastMeasure(x,dict,Omega)
%
% Apply the measurement matrix implicit operator.
%
% Input
%	x           input vector
%	dict        sensing matrix name
%	Omega	    nx1 vector denoting which measurements are kept
% Output
%	y	    observed data

if strcmp(upper(dict),'FOURIER')
	p = length(x);
	X = 1/sqrt(p)*fft(x);
	y = X(Omega);
elseif strcmp(upper(dict),'RST')
	X = RST(x);
	y = X(Omega);
elseif strcmp(upper(dict),'REALFOURIER') % Redundancy = 2, [real part; imag part] instead of complex Fourier measurements.
	p = length(x);
	X = sqrt(2/p)*fft(x);
	y = [real(X(Omega)); imag(X(Omega))];
elseif strcmp(upper(dict),'HADAMARD')
	X = FHT(x);
	y = X(Omega);
elseif strcmp(upper(dict),'DIRAC')
	y = x(Omega);
else
	disp('Uknown measurement matrix');
end

