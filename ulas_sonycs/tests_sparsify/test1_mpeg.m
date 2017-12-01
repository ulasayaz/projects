% tests COMPRESSION RATIO vs. PSNR using FFMPEG

%close all;
%clear all;

fn='cam512.png';
cs=33; % 1/ratio

'start'
%% test variables set here
L=256; % number of frames

videoflag=2; % 1 for rotation, 2 for shifts
compress=1; % 1 for mpeg4 video compression, 0 for raw video (no compression)

% if videoflag==2  % shift
%         init=ones(2,1);
% else % rotation
%         init=ones(1,1);
% end

% f=[0.1,0.5,1,1.5,3,5,10];
% set_units=kron(init, f);
% K=length(set_units);

data=char('100K','400K','700K','1M','1.3M'); % bitrates for ffmpeg % saturates at 100K and 1.3M
%use data(1,:)
K=length(data);
bitrate=data(1,:); %just to initialize

CR=zeros(1,K);
PS_mean=zeros(1,K);
SS_mean=zeros(1,K);

pixshift=1;
%rot=90/L;
%rot=5;
%rot=pixshift; % to rotate as much as 'pixshift'

if videoflag==2  % shift
    units=pixshift*ones(2,1);
else % rotation
    units=pixshift;
end

%%

signal2=Signal2D.make_fromImage(fn); % this scales the image

signal3=Signal3D.make_CyclefromLowerDimSig(signal2,videoflag,L,units);

global data3d;
data3d=signal3.xn;
M=size(data3d,1); N=size(data3d,2);

% MPEG compression with ffmpeg

global directory prefix pre_path video_path;
directory = '~/Desktop/ffmpeg_example/';
prefix = 'images/test_';
pre_path = [directory prefix];

system(['rm ' directory 'images/*']); %deletes everything in 'images/' folder
% write the original images from the data matrix
for i=1:L
    path = [pre_path int2str(i) '.png'];
    imwrite(data3d(:,:,i), path);
    if mod(i,20) == 0
        i
    end
end
% real size of the original video
c_ref=L*file_size(path);

%% First create the uncompressed video
compress=0;
%bitrate='4M';
ffmpeg_create(compress, videoflag, pixshift, bitrate);
%getting the file size of a file
%c_ref=file_size(video_path);
%%

%%% Now compress with different ratios
compress=1;

for k=1:K % loop over several bitrates
    k
    bitrate=data(k,:),
    ffmpeg_create(compress, videoflag, pixshift, bitrate); % also extract the data
    
    CR(k)=file_size(video_path); %file size of the compressed video
    
    % read extracted files
    ext3d=zeros(M,N,L);
    for i=1:L
        x=imread([directory 'extracted/test_' int2str(i)]);
        x=rgb2gray(x);
        x=normalize_frame(x); % is not normalized now
        ext3d(:,:,i)=x;
        if mod(i,20) == 0
            i
        end
    end
    %'psnr of extracted video'
    PS_mean(k)=PSNR3(data3d,ext3d);
    SS_mean(k)=SSIM3(data3d,ext3d);
end

f=c_ref./CR;


%--------
% %% VideoReader % properties: NumberOfFrames, FrameRate
% % read the compressed video created by ffmpeg
% if compress == 1 %with compression    
%     [p,s] = read_video(video_path)
%     'psnr of read ffmpeg video'
% end

%--------
% %% VideoWriter %properties: FrameRate, CompressionRatio, Quality
% %methods: open, close, writeVideo
% global new_path;
% quality=25; %default 75 (0-100)
% write_video(compress,quality);
% 
% % read the compressed video created by VideoWriter
% if compress == 1 %with compression    
%     [p,s] = read_video(new_path)
%     'psnr of read matlab video'
% end

'finish'




        





