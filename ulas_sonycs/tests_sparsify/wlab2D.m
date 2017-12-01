function [inv_xw] = wlab2D(x,cs,c,L)
% performs wavelet compression with Wavelab functions fwt2_po & iwt2_po
% ratio=compression rate in decimals
% c=wave type 
% L=depth

% L=3 by default
if nargin < 4 
    L=3;
end

%D4, C3, S8 wavelets

if c==1
    wave='Daubechies',num=4;
elseif c==2
    wave='Coiflet',num=3;
elseif c==3
    wave='Symmlet',num=8;
end
%title(wave);

qmf = MakeONFilter(wave,num);
xw = FWT2_PO(x,L,qmf);

xw_sort = sort(abs(xw(:)),'descend'); %sort in descening order

dim=length(x)^2;
wthresh = xw_sort(floor(dim/cs));

s_xw = xw .* (abs(xw) > wthresh);

inv_xw = IWT2_PO(s_xw,L,qmf);

end

