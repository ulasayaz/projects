% test matlab's 3D wavelet transformation using Gunter's code. 

%fn='cam512.png';
fn='cameraman.bmp';
%fn='cam256.png';
%fn='cam64.png';
cs=40; % 1/ratio

'start'
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

% signal3=Signal3D.make_CyclefromLowerDimSig(signal2,videoflag,L,units);

% ------
L=3; fn='tennis.avi';
%fn='riksch1.avi';
%fn='tennismen32_fewKB.avi';
signal3=Signal3D.make_fromVideo(fn,[4,3]);

% ------

% --- 3d sparsification with wavelets:
w3= Wavelet3D_mlab(signal3);
%w3.set_basisname('db1'); 
%w3.set_basisname('db2');
%w3.set_signal(signal3);
%lev=w3.deepestlev; 
%lev=log2(max(size(signal3.xn)));
%lev=w3.deepestlevAnisotropic;
%w3.dec(lev); %wavelet decomposition
w3.dec;

test3=w3.sparseApprox(cs); %compression and decomposition

PSNR3d=signal3.PSNR(test3);
SSIM3d=signal3.SSIM(test3);
mean(PSNR3d)
mean(SSIM3d)

y=normalize_frame(test3.xn);
playData(y);

'finish'