function y = FastDictOp1D(mode, n, p, x, It, dim)

global dict pars1 pars2 pars3

if mode == 1
    
    z = zeros(dim,1);
    z(It) = x;
    xr = FastS2(z,floor(sqrt(n)),dict,pars1,pars2,pars3);
    y = xr(:);
    
elseif mode == 2
    
    y = FastA2(reshape(x,floor(sqrt(n)),floor(sqrt(n))),dict,pars1,pars2,pars3);
    y = y(:);
    y = y(It);
    
end
