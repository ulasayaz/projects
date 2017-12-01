function [] = write_video(compress,q)
% writes a video from the data cube 'data3d' with Matlab's VideoWriter 
% writes in the '/ffmpeg/videos/' folder either compressed or uncompressed
% Motion JPEG method with given quality
global new_path directory data3d;

if compress == 0 %no compression
    ext = 'uncomp.avi';
    method='Uncompressed AVI';
else %with compression
    ext= 'comp.avi';
    method='Motion JPEG AVI';
end

clear newVid;
new_path=[directory 'videos/' ext];
newVid=VideoWriter(new_path,method);
newVid.FrameRate=25;

if compress == 1 %with compression
    newVid.Quality=q; %default 75 (0-100)
end

open(newVid);

L=size(data3d,3);
for j=1:L
    writeVideo(newVid,data3d(:,:,j)); %method of VideoWriter
    if mod(j,20) == 0
        j
    end
end
close(newVid);

end

