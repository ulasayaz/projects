function y = FastOp1D(mode, n, p, x, I, dim)

global H h


if mode == 1
    
    z = zeros(dim,1);
    z(I) = x;
    X = fft(z);
    y = real(ifft(X.*H));
    y = y(:);
    
elseif mode == 2
    
    X = fft(x);
    y = real(ifft(X.*conj(H)));
    y = y(:);
    y = y(I);
    
end%
