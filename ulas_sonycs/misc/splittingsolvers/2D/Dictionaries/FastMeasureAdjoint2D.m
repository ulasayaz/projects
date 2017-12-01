function x = FastMeasureAdjoint2D(y,p,dict,Omega)
%
% Apply the adjoint measurement matrix implicit operator.
%
% Input
%	y           input data
%	dict        sensing matrix name
%	Omega	    nx1 vector denoting which measurements are kept
% Output
%	y	    output image

global H

if strcmp(upper(dict),'FOURIER')
	p = floor(sqrt(p));
	X = zeros(p*p,1);
	X(Omega) = y;
	X = reshape(X,p,p);
	x = p*real(ifft2(X)); % Taking only the real part may accumulate numerical errors.
	x = x(:);
elseif strcmp(upper(dict),'RST')
	p = floor(sqrt(p));
	X = zeros(p*p,1);
	X(Omega) = y;
	X = reshape(X,p,p);
	x = Inv_RST2(X);
	x = x(:);
elseif strcmp(upper(dict),'REALFOURIER') % Redundancy = 2, [real part; imag part] instead of complex Fourier measurements.
	p = floor(sqrt(p));
	X = zeros(p*p,1);
	X(Omega) = (y(1:end/2) + i*y(end/2+1:end));
	X = reshape(X,p,p);
	x = p*real(ifft2(X)); % Taking only the real part may accumulate numerical errors.
	x = x(:);
elseif strcmp(upper(dict),'HADAMARD')
	p = floor(sqrt(p));
	X = zeros(p*p,1);
	X(Omega) = y;
	X = reshape(X,p,p);
	x = Inv_FHT2(X);
	x = x(:);
elseif strcmp(upper(dict),'CONV')
	p = floor(sqrt(p));
	X = zeros(p*p,1);
	X(Omega) = y;
	X = reshape(X,p,p);
	x = real(ifft2(fft2(X).*conj(H)));
	x = x(:);
elseif strcmp(upper(dict),'DIRAC')
	p = floor(sqrt(p));
	X = zeros(p*p,1);
	X(Omega) = y;
	x = X(:);
else
	disp('Uknown measurement matrix');
end







    
    
