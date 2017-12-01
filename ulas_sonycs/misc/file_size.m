function [c] = file_size(path)
% returns the size of a file in bytes

%double '' means ' in a string
file=['find ' path ' -type f -print0 | xargs -0 stat -f%z | awk ''{b+=$1} END {print b}''']

[logic,out]=system(file); %executes the command to get the size of the video file
c=sscanf(out, '%u'); % char to double


end

