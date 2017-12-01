% tests COMPRESSION RATIO vs. PSNR using Matlab Motion JPEG

close all;
%clear all;

fn='cam512.png';
%cs=33; % 1/ratio

'start'
%% test variables set here
L=128; % number of frames

videoflag=2; % 1 for rotation, 2 for shifts
compress=1; % 1 for mpeg4 video compression, 0 for raw video (no compression)


data=[1 10 20 30 40 50]; % quality for VideoWriter
K=length(data);

CR=zeros(1,K);
PS_mean=zeros(1,K);
SS_mean=zeros(1,K);

pixshift=5;
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

% Motion JPEG compression with Matlab

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

%% VideoWriter --properties: FrameRate, CompressionRatio, Quality
%methods: open, close, writeVideo

global new_path;
compress=1;

for k=1:K % loop over several qualities
    
    quality=data(k); %default 75 (0-100)
    write_video(compress,quality); % this creates the 'new_path'
    CR(k)=file_size(new_path); %file size of the compressed video
    
    % read the compressed video created by VideoWriter
    mov = read_video(new_path);
    PS_mean(k)=PSNR3(data3d,mov);
    SS_mean(k)=SSIM3(data3d,mov);
end
 
f=c_ref./CR;
%figure;
%plot(f,PS_mean)

%% TRY
% compress=1;
% quality=0.5;
% write_video(compress,quality);

'finish'




        





