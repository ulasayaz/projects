% this code compares the compressibility of wavelet and 
% curvelet transforms in the globally shifted video sequences
% we change the number of frames at every experiment. pixshift=1

%fn='cam512.png';
%fn='cameraman.bmp';
fn='cam256.png';

cs=33; % 1/ratio

'start'

pixshift=1;

videoflag=2; % 1 for rotation, 2 for shifts

signal2=Signal2D.make_fromImage(fn);

Frames=[32,64,128,256];

K=length(Frames);
    
PS_mean=zeros(2,K); % for wavelets and curvelets
SS_mean=zeros(2,K);

for k=1:K
    k
    L=Frames(k);
    
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

figure
plot(Frames,PS_mean(1,:),'Marker','o','LineWidth',2);
hold
plot(Frames,PS_mean(2,:),'Color','r','Marker','s','LineWidth',2);
title('Compressibility with wavelets and curvelets VS. number of frames (cs=33, pixel shift=1)','FontSize',12);
xlabel('number of frames','FontSize',12);
ylabel('PSNR');
%line([0 10],[29.19 29.19]); % constant line of curve2D
figure
plot(Frames,SS_mean(1,:),'Marker','o','LineWidth',2);
title('Compressibility with wavelets and curvelets VS. number of frames (cs=33, pixel shift=1)','FontSize',12);
xlabel('number of frames','FontSize',12);
ylabel('SSIM','FontSize',12);
hold
plot(Frames,SS_mean(2,:),'Color','r','Marker','s','LineWidth',2);

