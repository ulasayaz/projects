function  [wfinal]=sbhe(k,hd,mt)
%
% Function implementing measurement matrix as given in  
% "Fast compressive imaging using scrambled block Hadamard ensemble."
%  by Lu Gan, Thong Do and Trac Tran
%
% Written by Igor Carron.
%
% Input:
%     k : number of hadamard matrix blocks
%     hd: size of Hadamard matrix (needs to be divisible by 2)
%     mt: number of columns of wfinal
% Output:
%     wfinal: matrix with k*nd columns and mt rows
%
% k=10;
% hd = 32;
% mt = 30;
ntotal = hd * k;
wh = hadamard(hd);
st1='';
for i=1:k-1
	st1 = strcat(st1,'wh,');
end
st1= strcat(st1,'wh');
eval(['w = blkdiag(' st1 ');']);
u = randperm(ntotal);
for ii=1:ntotal
	wnew(:,ii)=w(:,u(ii));
end
%figure(1)
%spy(w)
%figure(2)
%spy(wnew)
vv = [];
while size(vv,1)~=mt
vv = [];
v2 = cat(2,vv,round(1+rand(mt,1)*(ntotal-1)));
vv = unique(v2);
end
m = size(vv,1);
for j=1:m
	wfinal(j,:) = wnew(vv(j),:);
end
%spy(wfinal)
