close all;
figure;

fn='cam512.png';
fn='cameraman.bmp';
%fn='cam256.png';

x=imread(fn);

%x = ReadImage('Daubechies');

%x=double(x);
x=normalize_frame(x);

subplot(2,2,1);
imgray(x); % displays pictures in Wavelab folder
title('original');

cs=33; % 1/ratio

%% CURVELET TRANSFORM
% there are two methods used
inv_xc=curve2D(x,cs,2); %mode=1 usfft, =2  wrapping
% curve2D does not allow allcurvelets=1

% Gunter's class uses fdct_wrapping_GT with different levels
signal2=Signal2D.make_fromImage(fn);
w= Curvelet2_clab(signal2);
w.set_deepestlev(2);
w= Curvelet2_clab(signal2);
w.set_deepestlev(2); w.dec;
test1=w.sparseApprox(cs);

subplot(2,2,2);

imgray(real(test1.xn));
%imgray(real(inv_xc));
title('curvelet, psnr=28.3, ssim=0.81');

%% WAVELET TRANSFORM

c=1;
% %D4, C3, S8 wavelets
% 
% if c==1
%     wave='Daubechies';
% elseif c==2
%     wave='Coiflet';
% elseif c==3
%     wave='Symmlet';
% end

%Z=floor(log2(size(x,1)));
L=3;

%hold;
subplot(2,2,3);

w2= Wavelet2D_mlab(signal2);
w2.set_basisname('db1');
lev=w2.deepestlev; 
%lev=w2.deepestlevAnisotropic;
w2.dec(lev); %wavelet decomposition
test2=w2.sparseApprox(cs);

%inv_xw=wlab2D(x,cs,c);
%imgray(inv_xw);

imgray(real(test2.xn));
title('wavelet, psnr=26.2, ssim=0.76.');

%% FOURIER TRANSFORM

xf = dct2_iv(x);
xf_sort = sort(abs(xf(:)),'descend');

fthresh = xf_sort(floor(numel(x)/cs));
s_xf = xf .* (abs(xf) > fthresh);
inv_xf = dct2_iv(s_xf);
%GrayImage(inf_xf,256);
subplot(2,2,4);
%figure
%AutoImage(inv_xf);
imgray(inv_xf);
title('Fourier, psnr=22.7, ssim=0.52');

%% PSNR and SSIM values

r=zeros(3,1);
%r(1)=psnr(x,inv_xc);
r(1)=signal2.PSNR(test1);
r(2)=signal2.PSNR(test2);
%r(2)=psnr(x,inv_xw);
r(3)=psnr(x,inv_xf);
'order = curvelet, wavelet, fft'
r

sim=zeros(3,1);
%sim(1)=SSIM3(x,inv_xc);
sim(1)=signal2.SSIM(test1);
sim(2)=signal2.SSIM(test2);
%sim(2)=SSIM3(x,inv_xw);
sim(3)=SSIM3(x,inv_xf);
sim

 MSE=1/numel(x)*sum((x(:)-inv_xc(:)).^2);
 max_x=max(abs(x(:)));
 pr=20*log10(max_x/sqrt(MSE));


