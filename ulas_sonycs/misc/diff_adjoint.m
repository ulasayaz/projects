function [xadj] = diff_adjoint(x)
% adjoint of diff() operator of Matlab

[M,N,L]=size(x);

xadj=zeros(M,N,L+1);

xadj(:,:,1)=-x(:,:,1);

for j=2:L
    xadj(:,:,j)=x(:,:,j-1)-x(:,:,j);    
end

xadj(:,:,L+1)=x(:,:,L);

end

