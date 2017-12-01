% quality of reconstruction TV_IC_MRA3
% FISTA, 500 it

method='BPDN(FISTA), it=500, \gamma=0.01, \tau=0.01';
optimizationProblem=3;
encoder='fft2, subsampling ratio of side frames=8';
meas_compression=8;
signalname='tennis 512x512xL';

framecount=[1,3,5,7,9,15];
PSNR_reached=fliplr([37.6, 36.8,36.2, 34.9,32.5,28.1]);
PSNR_reached_std=fliplr([0.6, 0.2,0.2,0.1,0.04,0]);
PSNR_bound=fliplr([42.3,41.5,40.7,39.2,36.5,30.4]);
PSNR_bound_std=fliplr([0.9,0.6,0.5, 0.3,0.2,0]);

SSIM_reached=fliplr([0.98, 0.97,0.97,0.96,0.93,0.86]);
SSIM_reached_std=fliplr([2e-3, 1e-3,1e-3,8e-4,1e-3,0]);
SSIM_bound=fliplr([0.98, 0.98,0.98,0.97,0.96,0.92]);

ppfig.pixsizeX=1000;
prepfigure(1,[],ppfig);

subplot(1,2,1);
errorbar(framecount,PSNR_reached,PSNR_reached_std);
hold all;
errorbar(framecount,PSNR_bound,PSNR_bound_std);
hold off;
legend('FISTA','DT-quality bound','location','best');
xlabel('video length L in frames','fontsize',12);
ylabel('PSNR [dB]','fontsize',12);

subplot(1,2,2);
errorbar(framecount,SSIM_reached,SSIM_reached_std);
hold all;
plot(framecount,SSIM_bound);
hold off;
legend('FISTA','DT-quality bound','location','best');
xlabel('video length L in frames','fontsize',12);
ylabel('SSIM','fontsize',12);

suptitle({[method,', optimization problem ',num2str(optimizationProblem),...
    ', mean compression c=',num2str(meas_compression,'%3.1f')],...
    ['encoder ',encoder,', video: ',signalname]},14);

