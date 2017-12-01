function mov = read_video(path)
% reads a video on the 'path' with Matlab's VideoReader class 
% and returns the data cube 'mov'
% ONLY works with compressed video

clear vid;
vid=VideoReader(path); % cannot read uncompressed video
H=vid.Height; W=vid.Width; F=vid.NumberOfFrames;
mov=zeros(H,W,F);
for i=1:F
    x=rgb2gray(read(vid,i));
    mov(:,:,i)=normalize_frame(x);
    if mod(i,20) == 0
        i
    end
end

end

