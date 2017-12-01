function y = aprod(mode, n, p, x, I, dim)

global A

if mode == 1
    
    z = zeros(dim,1);
    z(I) = x;
    y = A*z;
    y = y(:);
    
elseif mode == 2
    
    y = A'*x;
    y = y(:);
    y = y(I);
    
end

