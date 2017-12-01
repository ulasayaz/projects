% Runs all demo files in 1D/Demos and 2D/Demos.
% Produces all the figures and saves them in 1D/Datasets and 2D/Datasets.


%%%%%%%%%%%%
% 1D Demos.%
%%%%%%%%%%%%
disp('Starting 1D demos.');
disp('');

files = dir('1D/Demos');

for k=1:length(files)
    if (files(k).isdir==0)
    	fprintf('\t Running %s\n',files(k).name(1:end-2));
        eval(files(k).name(1:end-2));
    end
end

disp('');

%%%%%%%%%%%%
% 2D Demos.%
%%%%%%%%%%%%
disp('Starting 2D demos.');
disp('');

files = dir('2D/Demos');

for k=1:length(files)
    if (files(k).isdir==0)
    	fprintf('\t Running %s\n',files(k).name(1:end-2));
        eval(files(k).name(1:end-2));
    end
end
                                       
