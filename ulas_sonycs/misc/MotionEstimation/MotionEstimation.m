function VF = MotionEstimation (videodata, method,opts)
% This function outputs motion vector field for video data in the forward direction
% using one of the following methods
% method = 'OF' - Optical Flow
% method = 'FSBM' - Full Search Block Matching
%
% videodata contains input video of size H x B x L where L is the no. of
% video frames
%
% VF is the output data container (4-D array) of dimension H x B x L-1 x 3
% containing Horizontal Vector Field in 1st dim.,
% Vertical vector field in 2nd dim. and
% Motion Compensated frame in the 3rd dim.

if nargin <3 || isempty(opts)
    opts=struct;
end
sz = size(videodata);
L = sz(3);
VF = zeros (sz(1),sz(2), L-1, 3); % 1 for Vx, 2 for Vy, 3 for MC frame (4-D array of motion vectors)

label=['motion estimation ',method,' ...'];
multiWaitbar(label, 0); 

for j = 1:L-1,        
    
    img1 = videodata(:,:,j);
    img2 = videodata(:,:,j+1);
    if strcmp(method,'OF'),
        % transform image to grey values in [0,1]
        img1 = (img1-min(img1(:)))/(max(img1(:))-min(img1(:)));
        img2 = (img2-min(img2(:)))/(max(img2(:))-min(img2(:)));
        
        % set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
        if ~isfield(opts,'alpha') || isempty(opts.alpha)
            opts.alpha = 0.012;
        end
        if ~isfield(opts,'ratio') || isempty(opts.ratio)
            opts.ratio = 0.75;
        end
        if ~isfield(opts,'minWidth') || isempty(opts.minWidth)
            opts.minWidth = 20;
        end
        if ~isfield(opts,'nOuterFPIterations') || isempty(opts.nOuterFPIterations)
            opts.nOuterFPIterations = 7;
        end
        if ~isfield(opts,'nInnerFPIterations') || isempty(opts.nInnerFPIterations)
            opts.nInnerFPIterations = 1;
        end
        if ~isfield(opts,'nSORIterations') || isempty(opts.nSORIterations)
            opts.nSORIterations = 30;
        end
        
        para = [opts.alpha,opts.ratio,opts.minWidth,opts.nOuterFPIterations,...
            opts.nInnerFPIterations,opts.nSORIterations];
        
        % this is the core part of calling the mexed dll file for computing optical flow
        % it also returns the time that is needed for two-frame estimation
        %tic;
        [vx,vy,warpI2] = Coarse2FineTwoFrames(img1,img2,para);
        VF(:,:,j,1) = vx;
        VF(:,:,j,2) = vy;
        VF(:,:,j,3) = warpI2;
        % reconstruction with the other method FSBM:
        % VF(:,:,j,3) = reconstruct(img2, VF(:,:,j,1), VF(:,:,j,2), 1);
        
        %toc
        %figure;imagesc(im1)
        %figure;imagesc(im2)
        %figure;imagesc(warpI2)
        %clear flow;
        %flow(:,:,1) = vx;
        %flow(:,:,2) = vy;
        %imflow = flowToColor(flow);
        %figure;imshow(imflow);
    elseif strcmp(method,'FSBM'),
        % Parameters
        if ~isfield(opts,'BlockSize') || isempty(opts.BlockSize)
            opts.BlockSize   = 6; % 4 is better than 8 and 16 and same as 2
        end
        if ~isfield(opts,'SearchLimit')  || isempty(opts.SearchLimit)
            opts.SearchLimit = 20; % 20 same as 10 but better than 40
        end
        if ~isfield(opts,'subpixel')  || isempty(opts.subpixel)
            opts.subpixel = 0.25;% 0.25 is better than 0.5
        end
        % Full Search Block Matching Motion Estimation
        %tic
        [MVx, MVy] = Bidirectional_ME(img2, img1, opts); %to keep same signing convention to OF
        MVx = imresize(MVx, opts.BlockSize);
        MVy = imresize(MVy, opts.BlockSize);
        [M,N] = size(MVx);
        VF(1:M,1:N,j,1) = MVx;
        VF(1:M,1:N,j,2) = MVy;
        VF(:,:,j,3) = reconstruct(img2, VF(:,:,j,1), VF(:,:,j,2), opts.subpixel);
        % invert the sign like (-1*VFx,-1*VFy) for inverse motion compensation
        % MCr = reconstruct(odat(:,:,15), VF(:,:,14,1), VF(:,:,14,2), 0.25);
        % MCrInv = reconstruct(MCr, -1*VF(:,:,14,1), -1*VF(:,:,14,2), 0.25);
        %toc
    end
    multiWaitbar(label, j/(L-1)); 
end

multiWaitbar( label, 'Close' );

end