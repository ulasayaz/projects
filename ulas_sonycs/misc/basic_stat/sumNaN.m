function s = sumNaN(x,dim)
%sumNaN  sum neglecting NaNs
% Autor: G. Troll
% SMT 2012


if nargin==1,
    % Determine which dimension SUM will use
    dim = find(size(x)~=1,1,'first');
    if isempty(dim), dim = 1; end
    
end

x(isnan(x))=0;
if isvector(x)
    try
        s= double(kahan_sum(x));  % reduces the numerical error   
    catch
        s= sum(x,dim,'double');
    end
else    
    s= sum(x,dim,'double');  % compute sum in double precision
end

if isa(x,'single')
    s=single(s);
end


end
