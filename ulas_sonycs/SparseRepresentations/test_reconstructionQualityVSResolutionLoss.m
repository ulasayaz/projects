% resolution tests: compressibilty as a function of image details

%WavePath;
w=Wavelet2D_wlab();

fn='parrots.tiff';
fig=2;
opt.pixsizeX=1000;
%resize_method='bicubic';
resize_method='bilinear';


signal=Signal2D.make_fromImage(fn);
imgsize=signal.size;
signal.crop(min(imgsize));

N=signal.size(1);

L=25;
scale=linspace(min(floor(L/2),64),N,L);
f=scale/N;

L=length(f);

c=33;

PSNR=zeros(L,1);
SSIM=zeros(L,1);

for j=1:L
    s2=signal.clone;
    if abs(f(j)-1)>1e-9
        s2.coarsen(f(j),resize_method);
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
plot(f,PSNR);
title('PSNR of reconstruction','fontsize',12); 
xlabel(['resolution loss factor (',resize_method,')'],'fontsize',12);
ylabel('PSNR');

subplot(1,2,2); 
plot(f,SSIM); 
xlabel(['resolution loss factor (',resize_method,')'],'fontsize',12);
ylabel('SSIM');

title('SSIM of reconstruction','fontsize',12); 



ssmall=signal.clone;
ssmall.resize(256);
N2=256;

prepfigure(fig+1,[],opt);
L2=4;
sd=factor_subplots(L2);
scale2=linspace(64,N2,L2); 
f2=scale2/N2;
signal2=ssmall.clone;
for j=1:4
    s3=signal2.clone;
    if abs(f(j)-1)>1e-9
        s3.coarsen(f2(j),resize_method);
    end
    subplot(sd(1),sd(2),j);
    s3.graph_signal(false);
end







