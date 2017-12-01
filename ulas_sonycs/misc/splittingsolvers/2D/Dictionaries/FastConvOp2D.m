function y = FastConvOp2D(mode, n, p, x, I, dim)

global H h


if mode == 1
    
    z = zeros(dim,1);
    z(I) = x;
    X = fft2(reshape(z,floor(sqrt(dim)),floor(sqrt(dim))));
    y = real(ifft2(X.*H));
    y = y(:);
    
elseif mode == 2
    
    X = fft2(reshape(x,floor(sqrt(n)),floor(sqrt(n))));
    y = real(ifft2(X.*conj(H)));
    y = y(:);
    y = y(I);
    
end
