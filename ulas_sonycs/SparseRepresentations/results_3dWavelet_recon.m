% quality of reconstruction TV_IC_MRA3
% FISTA, 500 it

method='BPDN(DR), it=100, \gamma=1, \epsilon=0';
method_short='DouglasRachford recon';
optimizationProblem=0;
decoder='db1 in 3d';
encoder='fft3, subsampling ratio of side frames=8';
meas_compression=8;
signalname='tennis 512x512xL';

framecount=[1,4,8,12,16];
PSNR_reached=fliplr([29.2, 28.7,28.2, 26.4,28.1]);
PSNR_reached_std=fliplr([0.04, 0.04,0.04,0.04,0]);
PSNR_bound=fliplr([37.1,37.0,36.7,35.1,30.4]);
PSNR_bound_std=fliplr([0.04,0.1,0.1, 0.1,0]);

SSIM_reached=fliplr([0.80, 0.78,0.76,0.64,0.86]);
SSIM_reached_std=fliplr([5e-4, 6e-4,3e-4,3e-4,0]);
SSIM_bound=fliplr([0.98, 0.97,0.97,0.96,0.92]);

ppfig.pixsizeX=1000;
prepfigure(1,[],ppfig);

subplot(1,2,1);
errorbar(framecount,PSNR_reached,PSNR_reached_std);
hold all;
errorbar(framecount,PSNR_bound,PSNR_bound_std);
hold off;
legend(method_short,'DT-quality bound','location','best');
xlabel('video length L in frames','fontsize',12);
ylabel('PSNR [dB]','fontsize',12);

subplot(1,2,2);
errorbar(framecount,SSIM_reached,SSIM_reached_std);
hold all;
plot(framecount,SSIM_bound);
hold off;
legend(method_short,'DT-quality bound','location','best');
xlabel('video length L in frames','fontsize',12);
ylabel('SSIM','fontsize',12);

suptitle({[method,', optimization problem ',decoder,...
    ', mean compression c=',num2str(meas_compression,'%3.1f')],...
    ['encoder ',encoder,', video: ',signalname]},14);

