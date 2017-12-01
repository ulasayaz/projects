% this code compares the compressibility of wavelet and 
% curvelet transforms in the globally shifted video sequences
% we change the pixel shift value at every experiment. 

%fn='cam512.png';
%fn='cameraman.bmp';
fn='cam256.png';

cs=33; % 1/ratio

'start'
L=64; % number of frames

videoflag=2; % 1 for rotation, 2 for shifts

signal2=Signal2D.make_fromImage(fn);

F=[0.1,0.5,1,1.5,3,5,10];

K=length(F);
    
PS_mean=zeros(2,K); % for wavelets and curvelets
SS_mean=zeros(2,K);

for k=1:K
    k
    pixshift=F(k);
    
    %rot=90/L;
    %rot=5;
    rot=pixshift; % to rotate as much as 'pixshift'
    
    if videoflag==2  % shift
        units=pixshift*[1;1];
    else % rotation
        units=rot;
    end
    
    signal3=Signal3D.make_CyclefromLowerDimSig(signal2,videoflag,L,units);
    
    % --- 3d sparsification with wavelets:
    w3= Wavelet3D_mlab();
    
    w3.set_basisname('db2'); % may set basis to 'db1'
    w3.set_signal(signal3);
    %lev=w3.deepestlev;
    %lev=log2(max(size(signal3.xn)));
    lev=w3.deepestlevAnisotropic;
    w3.dec(lev); %wavelet decomposition
    
    test3=w3.sparseApprox(cs); %compression and decomposition
    
    PSNR3d=signal3.PSNR(test3);
    SSIM3d=signal3.SSIM(test3);
    PS_mean(1,k)=mean(PSNR3d);
    SS_mean(1,k)=mean(SSIM3d);
    
    % --- 3d sparsification with curvelets
    x=signal3.xn;
    inv_xc= clab3D_compress(x,cs,0); % mode=0 for default, mode=1 modified
    PS_mean(2,k)=PSNR3(x,inv_xc);
    SS_mean(2,k)=SSIM3(x,inv_xc);

end

% plot results 

figure
plot(F,PS_mean(1,:));
hold
plot(F,PS_mean(2,:),'Color','r');
line([0 10],[29.19 29.19]); % constant line of curve2D
figure
plot(F,SS_mean(1,:));
hold
plot(F,SS_mean(2,:),'Color','r');

