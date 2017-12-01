function y = FastCSDictOp2D(mode, n, p, x, It, dim)

global dictCS Omega redundancy dict pars1 pars2 pars3

if mode == 1
    
    z = zeros(dim,1);
    z(It) = x;
    xr = FastS2(z,floor(sqrt(p/redundancy)),dict,pars1,pars2,pars3);
    y = FastMeasure2D(xr(:),dictCS,Omega);
    y = y(:);
    
elseif mode == 2
    
    xm = FastMeasureAdjoint2D(x,p/redundancy,dictCS,Omega);
    y = FastA2(reshape(xm,floor(sqrt(p/redundancy)),floor(sqrt(p/redundancy))),dict,pars1,pars2,pars3);
    y = y(:);
    y = y(It);
    
end

