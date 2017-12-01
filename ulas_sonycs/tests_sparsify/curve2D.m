function [inv_xc] = curve2D(x,cs,mode,r)
% uses default curvelet transform for 2D with allcurvelets=0
%x = 2D data 
%r=1 for real, r=0 for complex. 1 for default
%mode=1 for usfft, mode=2 for wrapping

if mode==1
    forward=@fdct_usfft;
    inverse=@ifdct_usfft;
elseif mode==2
    forward=@fdct_wrapping;
    inverse=@ifdct_wrapping;
end

if nargin < 4
    r=1;
end

xc=forward(x,r); %for fdct_wrapping

% LN=max(0,floor(log2(max(size(x))))); % wmaxlev method of MultiResTrafo2
% lev=LN-obj.wcoarsestlev; % default: floor(log2(min(m,n)))-2)
% xc=forward(x,r,true,lev); % for fdct_wrapping_GT

% Get threshold value
cfs =[];
for s=1:length(xc)
    for w=1:length(xc{s})
        cfs = [cfs; abs(xc{s}{w}(:))];
    end
end
cfs = sort(cfs); cfs = cfs(end:-1:1);
nb = round(length(cfs)/cs);
cutoff = cfs(nb);

% set small coefficients to zero
for s=1:length(xc)
    for w=1:length(xc{s})
        xc{s}{w} = xc{s}{w} .* (abs(xc{s}{w})>cutoff);
    end
end

inv_xc = inverse(xc,r);

end

