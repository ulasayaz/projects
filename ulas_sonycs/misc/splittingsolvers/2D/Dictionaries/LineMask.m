function [M,Omega] = LineMask(L,n,method)
% LineMask.m
%
% Returns the indicator of the domain in 2D fourier space for the 
% specified line geometry.
% Usage :  [M,Omega] = LineMask(L,N)
%
% Written by : Jalal Fadili
% 

if ~exist('method','var'),
      method = 1;
end


if method==1
thc = linspace(0, pi-pi/L, L);
%thc = linspace(pi/(2*L), pi-pi/(2*L), L);

M = zeros(n);

% full mask
for ll = 1:L

      if ((thc(ll) <= pi/4) | (thc(ll) > 3*pi/4))
	yr = round(tan(thc(ll))*(-n/2:n/2-1))+n/2+1;
	for nn = 1:n
		      M(yr(nn),nn) = 1;
		end
      else 
	xc = round(cot(thc(ll))*(-n/2:n/2-1))+n/2+1;
	for nn = 1:n
	      M(nn,xc(nn)) = 1;
	end
      end

end

M = M(1:n,1:n);

else


Theta = linspace(0,pi,L+1); Theta(end) = [];
M = zeros(n,n);
t = linspace(-1,1,3*n)*n;
for theta = Theta
    x = round(t.*cos(theta)) + n/2+1; y = round(t.*sin(theta)) + n/2+1;
    I = find(x>0 & x<=n & y>0 & y<=n); x = x(I); y = y(I);
    M(y+(x-1)*n) = 1;
end

end

M = fftshift(M);
Omega = find(M);
