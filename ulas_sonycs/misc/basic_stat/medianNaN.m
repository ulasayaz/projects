function y = medianNaN(x)
%median value neglecting NaNs
% Autor: G. Troll
% SMT 2011


x=x(isfinite(x));
y=median(x);

end
