function y = FastMeasure2D(x,dict,Omega)
%
% Apply the measurement matrix implicit operator.
%
% Input
%	x           input image
%	dict        sensing matrix name
%	Omega	    nx1 vector denoting which measurements are kept
% Output
%	y	    observed data

global H

if strcmp(upper(dict),'FOURIER')
	p = floor(sqrt(length(x)));
	x = reshape(x,p,p);
	X = 1/p*fft2(x);
	X = X(:);
	y = X(Omega);
elseif strcmp(upper(dict),'RST') %Real Sinusoid Transform cf. sparselab
	p = floor(sqrt(length(x)));
	x = reshape(x,p,p);
	X = RST2(x);                %cf. sparselab
	X = X(:);
	y = X(Omega);
elseif strcmp(upper(dict),'REALFOURIER') % Redundancy = 2, [real part; imag part] instead of complex Fourier measurements.
	p = floor(sqrt(length(x)));
	x = reshape(x,p,p);
	X = 1/p*fft2(x);
	X = X(:);
	y = [real(X(Omega)); imag(X(Omega))];
elseif strcmp(upper(dict),'HADAMARD') 
	p = floor(sqrt(length(x)));
	x = reshape(x,p,p);
	X = FHT2(x);                % cf. sparselab
	X = X(:);
	y = X(Omega);
elseif strcmp(upper(dict),'CONV')
	p = floor(sqrt(length(x)));
	x = reshape(x,p,p);
	X = real(ifft2(fft2(x).*H));
    	X = X(:);
	y = X(Omega);
elseif strcmp(upper(dict),'DIRAC')
	y = x(Omega);
else
	disp('Uknown measurement matrix');
end

