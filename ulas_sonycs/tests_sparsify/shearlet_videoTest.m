tic;
fprintf('\n---SLExampleVideoDenoising---\n');
fprintf('loading video... ');

clear;

%%settings
sigma = 30;
scales = 2;
shearLevels = [1 1];
thresholdingFactor = 3;
directionalFilter = modulate2(dfilters('cd','d')./sqrt(2),'c');

%%load data
fn='cam64.png';
%fn='cameraman.bmp';
%fn='barbara.jpg';

cs=40;

L=64; % number of frames

videoflag=2; % 1 for rotation, 2 for shifts

signal2=Signal2D.make_fromImage(fn);

pixshift=1;

%rot=90/L;
%rot=5; 
rot=pixshift; % to rotate as much as 'pixshift'

if videoflag==2  % shift
        units=pixshift*[1;1];
else % rotation
        units=rot;
end

signal3=Signal3D.make_CyclefromLowerDimSig(signal2,videoflag,L,units);
X=signal3.xn;

% ------
%L=3; fn='tennis.avi';
%fn='riksch1.avi';
%fn='tennismen32_fewKB.avi';
%signal3=Signal3D.make_fromVideo(fn,[4,3]);


elapsedTime = toc;
fprintf([num2str(elapsedTime), ' s\n']);
tic;
fprintf('generating shearlet system... ');

%%create shearlets
shearletSystem = SLgetShearletSystem3D(0,size(X,1),size(X,2),size(X,3),scales,shearLevels,0,directionalFilter);

elapsedTime = toc;
fprintf([num2str(elapsedTime), ' s\n']);
tic;
fprintf('decomposition, thresholding and reconstruction... ');

%%decomposition
coeffs = SLsheardec3D(X,shearletSystem);

%%thresholding
c_sort = sort(abs(coeffs(:)),'descend'); %sort in descening order

dim=length(c_sort);
thresh = c_sort(floor(dim/cs));

coeffs = coeffs .* (abs(coeffs) > thresh);

%%reconstruction
Xrec = SLshearrec3D(coeffs,shearletSystem);

elapsedTime = toc;
fprintf([num2str(elapsedTime), ' s\n']);

%%compute psnr
PSNR = SLcomputePSNR(X,Xrec);

fprintf(['PSNR: ', num2str(PSNR) , ' db\n']);

y=normalize_frame(Xrec);
playData(y);