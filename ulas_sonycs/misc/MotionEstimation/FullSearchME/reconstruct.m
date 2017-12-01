function g = reconstruct(img0, MVx, MVy, pel)
% g=~(img0,MVx,MVy,[pel]) Integer pixel motion compensation;
% constructs a motion compensated frame of img0 according to the motion
% vectors specified by MVx and MVy
% if <<pel>> is missing no resize is done and blocksize=1 assumed.
% 
% Stanley Chan
% 29 Apr, 2010
% 10 Feb, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin >=4
    img = imresize(img0, 1/pel, 'bilinear');
    BlockSize  = floor(size(img,1)/size(MVx,1));
    [m, n, C]    = size(img);
    M          = floor(m/BlockSize)*BlockSize;
    N          = floor(n/BlockSize)*BlockSize;
    f          = img(1:M, 1:N, 1:C);
    
    MVxmap = imresize(MVx, BlockSize);
    MVymap = imresize(MVy, BlockSize);
    Dx = round(MVxmap*(1/pel));
    Dy = round(MVymap*(1/pel));
    g          = zeros(M, N, C);
    
    [xgrid, ygrid] = meshgrid(1:N, 1:M);
    
    X = min(max(xgrid+Dx, 1), N);
    Y = min(max(ygrid+Dy, 1), M); % bug fixed
    
    idx = (X(:)-1)*M + Y(:);
    for coloridx = 1:C
        fc = f(:,:,coloridx);
        g(:,:,coloridx) = reshape(fc(idx), M, N);
    end
    g = imresize(g, pel);
    
else  % assume img0, MVx, MVy all have the same size
    [M, N, C]=size(img0);
    g= zeros(M, N, C);
    [xgrid, ygrid] = meshgrid(1:N, 1:M);
    
    Dx = round(MVx);
    Dy = round(MVy);
    X = min(max(xgrid+Dx, 1), N);
    Y = min(max(ygrid+Dy, 1), M);
    
    idx = (X(:)-1)*M + Y(:);
    for coloridx = 1:C
        fc = img0(:,:,coloridx);
        g(:,:,coloridx) = reshape(fc(idx), M, N);
    end
    
end
