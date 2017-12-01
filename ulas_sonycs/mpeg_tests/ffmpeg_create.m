function [] = ffmpeg_create(compress, videoflag, pixshift, bitrate)
% writes video with ffmpeg function from an image folder and extracts it 
% cannot read a data matrix directly.
% compress = 1 (compressed mpeg4), 0 uncompressed
% videoflag = 1 for rotation, 2 for shifts
% pixshift = amount of shift or rotation
% bitrate= quality of compression [100K to 4M]

global directory prefix pre_path video_path;

% prefix = 'images/test_';
% pre_path = [directory prefix];

ffmpeg = '/Applications/ffmpeg/ffmpeg';

if videoflag==1 %rotate
    mode='rot';
else %shift
    mode='shift';
end

if compress == 0 %no compression
    codec = ' rawvideo';
    ending = '';
else %with compression
    codec = ' mpeg4';
    ending = '_comp';
end

%[] notation concatenates strings
video_path = [directory 'videos/test_' mode num2str(pixshift) ending '.avi'];

if compress == 0
    %no setting for -b:v
    build_command = [ffmpeg  ' -y -framerate 25 -start_number 1 -i ' pre_path '%d.png -vcodec' codec ' ' video_path];
else
    %with -b:v setting
    build_command = [ffmpeg  ' -y -framerate 25 -start_number 1 -i ' pre_path '%d.png -b:v ' bitrate ' -vcodec' codec ' ' video_path];
end

% -crf 10 (0-51 range, smaller better quality) 
% -b:v 500K to change bitrate. -b:v more effective.
% -y forces to overwrite the destination file
%system(['rm -f ' video_path]); % deletes the video file

% create video from the frames
system(build_command);

%% extract frames from the video
extract_command =[ffmpeg ' -i ' video_path ' -q:v 0 -f image2 ' directory 'extracted/test_%d'];
% deletes everything in "/extracted" folder and extracts frames from the video
system(['rm ' directory 'extracted/*']); 
system(extract_command);

end

