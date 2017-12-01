function [PSN3d_mean, SSIM3d_mean, params]=test_VideoSparsification(fn, params)
% ~(fn,[params]) compare 3d sparsification with 2d-sparsification
% for simulated video with global motion
% fn ... filename of video
% params.c...  compression rate for hardthresholding
% params.L  ... extension in 3d dimension
% e.g.
% fn='tennis.avi';
% params.c=33;
% test_VideoSparsification(fn, params);
%
requireOnly(exist(fn, 'file'),'local','file exists on search path');
if ~exist('params','var')
    params=struct;
end

if ~isfield(params,'c')
    params.c=33;  % compression rate
end
if ~isfield(params,'L')
    params.L=64;
end
if ~isfield(params,'basisname3d')
    params.basisname3d='db1';
end
if ~isfield(params,'fig')
    params.fig=4;
end
if params.fig<2
    params.fig=2;
end

signal3=Signal3D.make_fromVideo(fn,params.L);
idx=1;
signal2=Signal2D(signal3.xn(:,:,idx));
signal2.signalname=[signal3.signalname,', frame ',num2str(idx)];

use_mlab=license('test','Wavelet_Toolbox');

opt.pixsizeX=1000;
videoflag_descr='time lag';
videoflag_type='video stream';


if use_mlab
    w3= Wavelet3D_mlab();
else
    w3= Wavelet3D_wlab();
end
w3.set_basisname(params.basisname3d);
w3.set_signal(signal3);
lev=w3.deepestlev;
%lev=w3.deepestlevAnisotropic;
w3.dec(lev);

test3=w3.sparseApprox(params.c);
%test3.graph_signal;  %3d video reconstructed


PSNR3d=signal3.PSNR(test3);
SSIM3d=signal3.SSIM(test3);

PSNR3d_mean=mean(PSNR3d);
SSIM3d_mean=mean(SSIM3d);



% ---- compare with 2D sparsification:

if use_mlab
    w2= Wavelet2D_mlab(signal2);
else
    w2= Wavelet2D_wlab(signal2);
end
w2.set_basisname('db1');
w2.dec;
test2=w2.sparseApprox(params.c);
test2.signalname=[signal2.get_signalname,' sparse Approx. c=',num2str(params.c)];
test2.colormap_active='gray';
%test2.graph_signal; % 2d recon


PSNR2d=signal2.PSNR(test2);
SSIM2d=signal2.SSIM(test2);

%% results: first figure
%=======================================================================

prepfigure(params.fig,[],opt);

w3class=strrep(class(w3),'_','\_');
w2class=strrep(class(w2),'_','\_');

cellstr={[signal3.signalname,' ', videoflag_type], ...
    ['3d-PSNR=', num2str(PSNR3d_mean,'%3.2f'),', 2d-PSNR=',num2str(PSNR2d,'%3.2f')],...
    ['3d-SSIM=', num2str(SSIM3d_mean,'%3.2f'),', 2d-SSIM=',num2str(SSIM2d,'%3.2f')] ...
    };
suptitle(cellstr, 14);

subplot(2,2,1);
section='1,:,:';
signal3_cut=signal3.cross_section(section);
signal3_cut.graph_signal(false);
titstr=signal3_cut.signalname;
titstr=strrep(titstr,signal2.signalname,'');
title(titstr);

subplot(2,2,2);
w3.fontsize=10;
w3.graph_detail(1,4,false);


%% second part of figure output
%=======================================================================


subplot(2,2,3);
w3.graph_detcoefNorms(false);


subplot(2,2,4);
coeff=w3.ts.make_like();
coeff.replace_signal(w3.C2vec);
coeff.graph_sparsityDefect(false,params.c);


end





