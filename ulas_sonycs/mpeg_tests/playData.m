function [] = playData(x)
% creates uncompressed video directly from a data cube and saves in the
% folder 'ffmpeg_example/videos' on the desktop.
% using write_video function (Matlab's VideoWriter class)

% needs the data $x$ to be normalized

global data3d directory;
directory = '~/Desktop/ffmpeg_example/';

data3d = x;
write_video(0); % uncompressed video titled "uncomp.avi"

end

