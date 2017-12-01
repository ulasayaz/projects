function sd = factor_subplots( L )
% factorize L plots into subplots

sd=[];
if L>=1
    fac=factor(L);
    if length(fac)>2
        f=cumprod(fac);
    else
        f=fac;
    end
    r=sqrt(L);
    idx=find(f>=r,1);
    sd=zeros(1,2);
    sd(2)=f(idx);
    sd(1)=L/sd(2);
    if sd(1)==1
        sd(2)=ceil(r);
        sd(1)=ceil(L/sd(2));
    end
end

end

