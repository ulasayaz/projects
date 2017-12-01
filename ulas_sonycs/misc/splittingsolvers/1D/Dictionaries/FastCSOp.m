function y = FastCSOp(mode, n, cardI, x, I, dim)

global dict Omega


if mode == 1
    
    z = zeros(dim,1);
    z(I) = x;
    y = FastMeasure(z,dict,Omega);
    y = y(:);
    
elseif mode == 2
    
    y = FastMeasureAdjoint(x,dim,dict,Omega);
    y = y(I);
    
end

