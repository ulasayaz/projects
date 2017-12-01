% testing shift invariance of SSIM
% ssim is invariant w.r.t. rescaling of intensities
% BUT is NOT invariant w.r.t. shift of intensities !!!

% in windows x,y of default size of 8x8
% SSIM(x,y)= (2*mux*muy+c1)/(mux^2+muy^2+c2) *
% (2*sigxy+c2)/(sigx^2+sigy^2+c2);
% with small deniminator stabilizers  ci=(ki*L)^2;
% k1=0.01; k2=0.03; L ... dynamic range.
% symmetries of SSI:
%-------------------
% SSIM(f*x, f*y)= SSIM(x,y)
% BUT translation symmetry is missing:  SSIM(x-t, y-t) != SSIM(x,y)

load('testingSSIM.mat');

x=Signal2D(xmat);
xorig=Signal2D(xorigmat);
xopt=Signal2D(xoptmat);

disp(x.stat);


x2=clone(x); xorig2=clone(xorig);
x2.shiftSignalValues(xorig.min);
xorig2.shiftSignalValues(xorig.min);
disp(x2.stat);

SSIM_translation0= xorig.SSIM(x); 
disp(['SSIM(shift=0)=',num2str(SSIM_translation0)]);

SSIM_translation1= xorig2.SSIM(x2); 
disp(['SSIM(shift to [0,1])=',num2str(SSIM_translation1)]);

% first factor of SSIM:
% (second factor containing variances and covariances is shift invariant)
ssim1=@(mux,muy) (2*mux.*muy + 1.e-4)./(mux.^2+muy.^2 + 1e-4);

t=linspace(0,0.49);
figure, plot(t, ssim1(0.49-t,0.5-t),'-r');
xlabel('shift t'); ylabel('ssim factor 1');
title('missing shift invariance of ssim');


