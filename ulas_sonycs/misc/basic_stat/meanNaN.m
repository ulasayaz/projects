function [y,N] = meanNaN(x,dim)
%MEAN   Average or mean value neglecting NaNs
% Autor: G. Troll
% SMT 2011


if nargin==1,
    % Determine which dimension SUM will use
    dim = min(find(size(x)~=1));
    if isempty(dim), dim = 1; end
    
end

h=isfinite(x);
N=sum(h,dim);
N(N==0)=NaN;

x(~h)=0;

if isvector(x)
    try
        y= double(kahan_sum(x))/N;  % reduces the numerical error  
    catch
        y= sum(x,dim,'double')./N;
    end
else    
    y= sum(x,dim,'double')./N;  % compute sum in double precision
end

if isa(x,'single')
    y=single(y);
end


end
