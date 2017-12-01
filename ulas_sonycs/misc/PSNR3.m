function [m,r] = PSNR3(orig,test)
% psnr values of two data cubes
% returns the mean and the full array 

L=size(orig,3);

r=zeros(L,1);

for j=1:L
    r(j)=psnr(orig(:,:,j),test(:,:,j));
end

m=mean(r);

end

