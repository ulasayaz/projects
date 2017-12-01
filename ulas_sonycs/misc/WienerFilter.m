function [xfiltered] = WienerFilter(x,N)
% x = 3 dim data
% WienerFilter filters each frame separately with 
% Matlab function  wiener2 with window size (N by N)


xfiltered = zeros([size(x)]);

L = size(x,3);

for k=1:L
    xfiltered(:,:,k) = wiener2(x(:,:,k),[N N]);
end

end

