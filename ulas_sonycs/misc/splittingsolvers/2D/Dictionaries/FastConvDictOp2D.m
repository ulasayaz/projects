function y = FastConvDictOp2D(mode, n, p, x, It, dim)

global h H dict pars1 pars2 pars3

if mode == 1
    
    z = zeros(dim,1);
    z(It) = x;
    xr = FastS2(z,floor(sqrt(n)),dict,pars1,pars2,pars3);
    X = fft2(xr);
    y = real(ifft2(X.*H));
    y = y(:);
    
elseif mode == 2
    
    X = fft2(reshape(x,floor(sqrt(n)),floor(sqrt(n))));
    xf = real(ifft2(X.*conj(H)));
    y = FastA2(xf,dict,pars1,pars2,pars3);
    y = y(:);
    y = y(It);
    
end
