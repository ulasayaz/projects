function [m,r] = SSIM3(orig,test)
% ssim values of two data cubes
% returns the mean and the full array 

L=size(orig,3);

r=zeros(L,1);

K = [0.01 0.03];
window = fspecial('gaussian', 11, 1.5);	
LL=1;

for j=1:L
    r(j)=ssim(orig(:,:,j),test(:,:,j),K,window,LL);
end

m=mean(r);

end

