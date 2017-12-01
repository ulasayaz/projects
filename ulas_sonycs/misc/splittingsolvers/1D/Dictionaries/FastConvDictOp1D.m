function y = FastConvOp1D(mode, n, p, x, It, dim)

global h H dict pars1 pars2 pars3

if mode == 1
    
    z = zeros(dim,1);
    z(It) = x;
    xr = FastS1(z,n,dict,pars1,pars2,pars3);
    X = fft(xr);
    y = real(ifft(X.*H));
    y = y(:);
    
elseif mode == 2
    
    X = fft(x);
    xf = real(ifft(X.*conj(H)));
    y = FastA1(xf,dict,pars1,pars2,pars3);
    y = y(:);
    y = y(It);
    
end
