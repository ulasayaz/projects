function [mov] = read_withFFmpeg(video_path,scale)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ffmpeg = '/Applications/ffmpeg/ffmpeg';
directory = '~/Desktop/ffmpeg_example/';

if ~exist('scale','var')
    scale=0.25;
end

% extract frames from the video
extract_command =[ffmpeg ' -i ' video_path ' -q:v 0 -f image2 ' directory 'extracted/test_%d'];
% deletes everything in "/extracted" folder and extracts frames from the video
system(['rm ' directory 'extracted/*']); 
system(extract_command);

[logic,out]=system(['ls -1 ' directory 'extracted | wc -l']); %number of extracted frames
L=sscanf(out, '%u'); % char to double

mov=[];
for i=1:L
    x=imread([directory 'extracted/test_' int2str(i)]);
    x=double(rgb2gray(x));
    %x=normalize_frame(x);
    x=imresize(x,scale);
    mov(:,:,i)=x;
    if mod(i,20) == 0
        i
    end
end
mov=normalize_frame(mov);
system(['rm ' directory 'extracted/*']); 
 
end

