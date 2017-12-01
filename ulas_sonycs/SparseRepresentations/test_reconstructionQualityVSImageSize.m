% resolution tests: compressibility as a function of image size

%WavePath;
w=Wavelet2D_mlab();

% fn='parrots.tiff';

fn='house.png';
fig=3;  % total is 3 images
opt.pixsizeX=1000;
resize_method='bicubic';
%resize_method='bilinear';


signal=Signal2D.make_fromImage(fn);
imgsize=signal.size;
signal.crop(min(imgsize));


%% resize image and test sparseApprox
N=signal.size(1);
L=25;
scale=linspace(64,N,L);
f=scale/N;

L=length(f);

%c=33;
c=89;

PSNR=zeros(L,1);
SSIM=zeros(L,1);

for j=1:L
    s2=signal.clone;
    if abs(f(j)-1)>1e-9
        s2.resize(f(j),resize_method);
    end
    w.set_signal(s2);
    w.dec;
    s2=w.sparseApprox(c);
    recon=s2.xn;
    PSNR(j)=s2.PSNR(recon);
    SSIM(j)=s2.SSIM(recon);
    fprintf('.');
    
end

disp(' ');

prepfigure(fig,[],opt);

suptitle([signal.get_signalname,', size=',vec2str(signal.size),', ',...
    'compression(Wavelets) c=', num2str(c,'%3.1f')],14);

subplot(1,2,1);
plot(scale,PSNR);
title('PSNR of reconstruction','fontsize',12); 
xlabel(['image size (',resize_method,')'],'fontsize',12);
ylabel('PSNR');

subplot(1,2,2); 
plot(scale,SSIM); 
xlabel(['image size (',resize_method,')'],'fontsize',12);
ylabel('SSIM');

title('SSIM of reconstruction','fontsize',12); 

waitforbuttonpress ;

%% example of images
ssmall=signal.clone;
ssmall.resize(256);
N2=256;

prepfigure(fig,[],opt);
L2=4;
sd=factor_subplots(L2);
scale2=linspace(64,N2,L2); 
f2=scale2/N2;
signal2=ssmall.clone;
for j=1:4
    s3=signal2.clone;
    if abs(f(j)-1)>1e-9
        s3.resize(f2(j),resize_method);
    end
    subplot(sd(1),sd(2),j);
    s3.graph_signal(false);
end

waitforbuttonpress; 

%% embed smaller picture in larger one

fn='cameraman.bmp';
signal=Signal2D.make_fromImage(fn);
s0=signal;
s0.ison_colorbar=false;
Nold=signal.size;

w=Wavelet2D_mlab();

% -----  sparse approx. of embedded picture:
% padd image to large size
s1=signal.clone;
s1.ison_colorbar=false;
N=2048; s1.set_0padding((N-imgsize(1))/2);
w.set_signal(s1);
w.dec;
c=89;
recon1=w.sparseApprox(c);
recon1.signalname=['Wavelet of 0padded', ' c=', num2str(c)];
% remove padding again
recon1.crop(Nold);
recon1.colormap_active=s0.colormap_active;
recon1.ison_colorbar=false;
PSNR1=s0.PSNR(recon1);
SSIM1=s0.SSIM(recon1);


% -----  sparse approx. of original picture:
w.set_signal(s0);
w.dec;
c=89;
recon0=w.sparseApprox(c);
recon0.colormap_active=s0.colormap_active;
recon0.signalname=['Wavelet', ' c=', num2str(c)];
recon0.ison_colorbar=false;
PSNR0=s0.PSNR(recon0);
SSIM0=s0.SSIM(recon0);


% ---- show results
opt.pixsizeX=1000;
prepfigure(fig, [], opt);
subplot(2,3,1);
s0.graph_signal(false);
xlabel(''); ylabel('');

subplot(2,3,3); 
recon0.graph_signal(false);
ylabel('');
xlabel(['PSNR=', num2str(PSNR0,'%3.1f'),', SSIM=', num2str(SSIM0,'%3.2f')], ...
    'fontsize',12);


subplot(2,3,4);
s0.graph_signal(false);
xlabel(''); ylabel('');
subplot(2,3,5);
s1.graph_signal(false);
xlabel(''); ylabel('');
subplot(2,3,6);
recon1.graph_signal(false);
axis image;
ylabel('');
xlabel(['PSNR=', num2str(PSNR1,'%3.1f'),', SSIM=', num2str(SSIM1,'%3.2f')], ...
    'fontsize',12);

waitforbuttonpress; 

%% tile image

w=Wavelet2D_mlab();
fig=3;

fn='cameraman.bmp';
signal=Signal2D.make_fromImage(fn);
s0=signal;
s0.ison_colorbar=false;
Nold=signal.size;


% ---- sparse approx. of tiled image
s1=s0.clone;
n=2048./s0.N;
s1.tile(n(1),n(2));
s1.signalname=[s0.signalname, ' tiled'];
s1.ison_colorbar=false;

w.set_signal(s1);
w.dec;
c=89;
recon1=w.sparseApprox(c);
recon1.xn=recon1.xn(1:s0.size(1), 1:s0.size(2));
recon1.colormap_active=s0.colormap_active;
recon1.ison_colorbar=false;
recon1.signalname=['Wavelet tiled', ' c=', num2str(c)];

PSNR1=s0.PSNR(recon1);
SSIM1=s0.SSIM(recon1);


% -----  sparse approx. of original picture:
w.set_signal(s0);
w.dec;
c=89;
recon0=w.sparseApprox(c);
recon0.colormap_active=s0.colormap_active;
recon0.signalname=['Wavelet', ' c=', num2str(c)];
recon0.ison_colorbar=false;
PSNR0=s0.PSNR(recon0);
SSIM0=s0.SSIM(recon0);



% ---- show results
opt.pixsizeX=1000;
prepfigure(fig, [], opt);
subplot(2,3,1);
s0.graph_signal(false);
xlabel(''); ylabel('');

subplot(2,3,3); 
recon0.graph_signal(false);
ylabel('');
xlabel(['PSNR=', num2str(PSNR0,'%3.1f'),', SSIM=', num2str(SSIM0,'%3.2f')], ...
    'fontsize',12);


subplot(2,3,4);
s0.graph_signal(false);
xlabel(''); ylabel('');
subplot(2,3,5);
s1.graph_signal(false);
xlabel(''); ylabel('');
subplot(2,3,6);
recon1.graph_signal(false);
axis image;
ylabel('');
xlabel(['PSNR=', num2str(PSNR1,'%3.1f'),', SSIM=', num2str(SSIM1,'%3.2f')], ...
    'fontsize',12);








