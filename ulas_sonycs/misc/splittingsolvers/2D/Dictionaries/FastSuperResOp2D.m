function y = FastSuperResOp2D(mode, n, p, x, I, dim)

global H h Omega


if mode == 1
    
    z = zeros(dim,1);
    z(I) = x;
    X = fft2(reshape(z,floor(sqrt(dim)),floor(sqrt(dim))));
    y = real(ifft2(X.*H));
    y = y(:);
    y = y(Omega);
    
elseif mode == 2
    
    p = floor(sqrt(p));
    z = zeros(p*p,1);
    z(Omega) = x;
    X = fft2(reshape(z,p,p));
    y = real(ifft2(X.*conj(H)));
    y = y(:);
    y = y(I);
    
end
