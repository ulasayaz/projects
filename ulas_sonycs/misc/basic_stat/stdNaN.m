function y = stdNaN(x, mu, dim)
% std of s with possible NaNs
% x and the mean mu must have the same dimensions (size)
% NaNs are acceptable
% Autor: G. Troll
% SMT 2011

xdim=length(size(x));
if nargin<=1 && isvector(x)
    mu=meanNaN(x)*ones(size(x));
end
if nargin<3,
    % Determine which dimension SUM will use
    dim = min(find(size(x)~=1));
    if isempty(dim), dim = 1; end
end
ok=all(size(x)==size(mu));
if  ~ok && xdim<=2 && isvector(mu)
    % bringe mu in die richtige Form
    L=length(mu);
    if dim==1 && L==size(x,2)
        mu=repmat(mu,size(x,1),1);
    elseif dim==2 && L==size(x,1)
        mu=repmat(mu,1,size(x,2));
    end
end
if  ~ok && ~all(size(x)==size(mu))
    % Vorbedingung immer noch verletzt
    error('stdNaN','arguments must have the same dimensions');
end

if ~isa(x,'double')
    x=double(x);
    mu=double(mu);
end

h=isfinite(x-mu);
N=sum(h,dim);
N(N==0)=NaN;


x(~h)=0;
mu(~h)=0;

if isvector(x)
    try
        y=( kahan_sum((x-mu).^2) ./(N-1) ).^0.5;  % reduces the numerical error
    catch
        y=( sum((x-mu).^2,dim,'double') ./(N-1) ).^0.5;
    end
else
    y=( sum((x-mu).^2,dim,'double') ./(N-1) ).^0.5;% compute sum in double precision
end


end

