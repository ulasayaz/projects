function [y,t] = perform_l1ball_projection(x,lambda)


% perform_l1ball_projection - compute the projection on the L1 ball
%
%   y = perform_l1ball_projection(x,lambda);
%
%   x is the projection of y on the set {a \ sum_i |a_i| = lambda }
%
%   Copyright (c) 2007 Gabriel Peyre

if lambda<0
    error('lambda should be > 0');
end
if lambda==0
    y = x*0;
    t = 0;
    return;
end

n = length(x(:));
% compute the thresholded L1 norm at each sampled value
s0 = sort( abs(x(:)) );
s = cumsum( s0(end:-1:1) ); s = s(end:-1:1);
s = s - s0 .* (n:-1:1)';
% compute the optimal threshold by interpolation
[i,tmp] = max( find(s>lambda) );
if isempty(i)
    y = x; 
    t = 0;
    return;
end
i = i(end);
t = ( s(i+1)-lambda )/( s(i+1)-s(i) ) * (s0(i)-s0(i+1)) + s0(i+1);
% do the actual thresholding
y = x;
y = (abs(x) > t) .* (x - t.*sign(x));
%y(abs(x)<t) = 0;
%y(abs(x)>=t) = y(abs(x)>=t) - sign(x(abs(x)>=t))*t;




