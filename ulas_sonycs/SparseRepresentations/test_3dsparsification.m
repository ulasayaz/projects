function [PSN3d_mean, SSIM3d_mean, params]=test_3dsparsification(fn, params)
% ~(fn,[params]) compare 3d sparsification with 2d-sparsification
% for simulated video with global motion
% fn ... filename of image
% params.c...  compression rate for hardthresholding
% params.videoflag ... 1 is rotation, 2 is shift
% params.L  ... extension in 3d dimension
% params.gen3d ... generation method for 3d signal ('cycle' or 'pan')
% --- e.g. for params.gen3d='cycle' use small images like
% fn='Mondrian.tif'; fn='cameraman.bmp'; 
% --- for params.gen3d='pan' use larger image like
% fn='parrots.tiff'; params.pansize=256*[1,1]; params.origin=[700,300]; params.L=8;
%
% params.c=33;
% params.videoflag=2;  % 1 is rotation, 2 is shift
% test_3dsparsification(fn, params);
%
requireOnly(exist(fn, 'file'),'local','file exists on search path');
if ~exist('params','var')
    params=struct;
end

if ~isfield(params,'videoflag')
    params.videoflag=2;  % shift
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
if ~isfield(params,'units')
    params.units=1;
end
if ~isfield(params,'gen3d')
    params.gen3d='cycle';
end
if ~isfield(params,'pansize')
    params.pansize=[];
end
if ~isfield(params,'origin')
    params.origin=[];
end
if ~isfield(params,'anisotropic')
    % if true, choose deepest possible decomposition level of largest
    % dimensional extension;
    % if false, choose deepest possible level of smallest dimensional
    % extension.
    params.anisotropic=true;
end


signal2=Signal2D.make_fromImage(fn);
ssig2=signal2.size;
use_mlab=license('test','Wavelet_Toolbox');
L=params.L;
if ~use_mlab
    % wlab can only handle dyadic equi-dimensions:
    signal2.crop();  % crop to lower dyadic bound
    if ~isequal(ssig2,signal2.size)
        warning('using wavelab toolbox required cropping of the image');
    end
    L=min(signal2.size);
    params.anisotropic=false;
end

opt.pixsizeX=1000;

% motion units:
units=params.units(:);
if params.videoflag==2  % shift
    if length(params.units)==1
        units=units*ones(2,1);
    end
else % rotation
    if length(units)>1
        units=max(abs(units));
    end
end

f=[0.1,0.5,1,1.5,3,5,10];
set_units=kron(units, f);
K=length(set_units);

PSN3d_mean=zeros(1,K);
SSIM3d_mean=zeros(1,K);

wait_handle = waitbar(0,...
                ['loop over magnitudes of motion (K=',num2str(K),') ...']);
% loop over several magnitudes of motion
for k=1:K
    
    waitbar(k / K);
    units=set_units(:,k);
    
    if params.videoflag==1
        videoflag_descr=['rotation (unit ',num2str(units(1),'%3.1f'),'°)'];
        videoflag_type='rotation (°)';
    else
        videoflag_descr=['global shift (units',vec2str(units,'%3.1f'),')'];
        videoflag_type='global shift (pix)';
    end
    
    % --- 3d sparsification:
    if strcmpi(params.gen3d,'pan')  % slower!
        signal3=Signal3D.make_PanfromLowerDimSig(signal2,L,params.pansize,params.origin,units);       
    else
        signal3=Signal3D.make_CyclefromLowerDimSig(signal2,params.videoflag,L,units);
    end
   
    % prepfigure(1,[],opt); signal3.graph_signal([],false);
    % suptitle([signal3.get_signalname,': video-sequence'],14);
    if use_mlab
        w3= Wavelet3D_mlab();
    else
        w3= Wavelet3D_wlab();
        ssig3=signal3.size;
        signal3.crop();
        if ~isequal(ssig3,signal3.size)
            warning('using wavelab toolbox required cropping of the image');
        end
        L=min(signal2.size);
    end
    w3.set_basisname(params.basisname3d);
    w3.set_signal(signal3);
    
    if  params.anisotropic
        lev=w3.deepestlevAnisotropic;
    else
        lev=w3.deepestlev;
    end
    
    w3.dec(lev);
    if k==1
        signal3_first=signal3;
        w3_first=w3;
    end
    
    
    test3=w3.sparseApprox(params.c);
    if K==1
        test3.graph_signal;  %3d video reconstructed
        % test3.graph_signal(1); title(videoflag_descr,'fontsize',12); % first reconstructed image of video
        
    end
    
    PSNR3d=signal3.PSNR(test3);
    SSIM3d=signal3.SSIM(test3);
    
    PSN3d_mean(k)=mean(PSNR3d);
    SSIM3d_mean(k)=mean(SSIM3d);
    
    rotangles=[0,signal3.rotangles(L)];
    shiftpix=[[0;0],signal3.shiftdist(L,units)];
    
    if params.videoflag==1
        xvals=rotangles;
        xdescr='rot. angle \alpha [°]';
    else
        xvals=(sum(shiftpix.^2,1)).^0.5;
        xdescr='shift [pix]';
    end
    
    
end
close(wait_handle );

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
if K==1
    test2.graph_signal; % 2d recon
end


PSNR2d=signal2.PSNR(test2);
SSIM2d=signal2.SSIM(test2);

%% results: first figure
%=======================================================================
if params.anisotropic
    isotropy_str='anisotropic dec.';
else
    isotropy_str='isotropic dec.';
end

prepfigure(params.fig,[],opt);
if K==1
    
    
    suptitle([signal3.get_signalname,', size=',vec2str(signal3.size),', ',...
        'compression(Wavelets) c=', num2str(params.c,'%3.1f'),', ',videoflag_descr],14);
    
    subplot(2,2,1);
    plot(xvals,signal3.PSNR(test3));
    title(['3d sparsification PSNR  (PSNR(2d)= ',num2str(PSNR2d,'%4.2g'),')'],'fontsize',12);
    xlabel(xdescr,'fontsize',10);
    ylabel('PSNR');
    
    subplot(2,2,2);
    plot(xvals,signal3.SSIM(test3));
    xlabel(xdescr,'fontsize',10);
    ylabel('SSIM');
    
    title(['3d sparsification SSIM (SSIM(2d)= ',num2str(SSIM2d,'%4.3g'),')'],'fontsize',12);
else
    
    norm_units=zeros(1,K);
    for k=1:K
        norm_units(k)=norm( set_units(:,k));
    end
    
    w3class=strrep(class(w3),'_','\_');
    w2class=strrep(class(w2),'_','\_');
    
    suptitle([signal2.signalname,', size=',vec2str(signal3.size),', ',...
        isotropy_str,', ', videoflag_type], 14);
    
    subplot(2,4,1:2);
    semilogx(norm_units, PSN3d_mean);
    line([norm_units(1), norm_units(end)],[PSNR2d,PSNR2d],'color','red');
    legend(['3d: mean PSNR of ',num2str(L),' frames'],['2d: ',w2class],...
        'location','best');
    title(['PSNR (',w3class,'(',w3.basisname,'), compression c=',num2str(params.c,'%3.0f'),')'],'fontsize',12);
    xlabel(videoflag_type,'fontsize',10);
    ylabel('PSNR');
    
    subplot(2,4,3:4);
    semilogx(norm_units, SSIM3d_mean);
    line([norm_units(1), norm_units(end)],[SSIM2d,SSIM2d],'color','red');
    legend(['3d: mean SSIM of ',num2str(L),' frames'],['2d: ',w2class],...
        'location','best');
    title(['SSIM (',w3class,'(',w3.basisname,'), compression c=',num2str(params.c,'%3.0f'),')'],...
        'fontsize',12);
    xlabel(videoflag_type,'fontsize',10);
    ylabel('SSIM');
    
end

subplot(2,4,5);
section='1,:,:';
signal3_first_cut=signal3_first.cross_section(section);
signal3_first_cut.graph_signal(false);
titstr=signal3_first_cut.signalname;
titstr=strrep(titstr,signal2.signalname,'');
title(titstr);

subplot(2,4,6);
w3_first.fontsize=10;
dsig=Signal3D(w3_first.detcoef(1,4));
dsig.signalname='Detail signal 1,4, z=1';
dsig.graph_signal(false,1);
  

subplot(2,4,7);
signal3_last_cut=signal3.cross_section(section);
signal3_last_cut.graph_signal(false);
titstr=signal3_last_cut.signalname;
titstr=strrep(titstr,signal2.signalname,'');
title(titstr);

subplot(2,4,8);
w3.fontsize=10;
dsig=Signal3D(w3.detcoef(1,4));
dsig.signalname='Detail signal 1,4, z=1';
dsig.graph_signal(false,1);



%% second figure
%=======================================================================

prepfigure(params.fig,[],opt);

subplot(2,2,1);
w3_first.graph_detcoefNorms(false);

subplot(2,2,2);
w3.graph_detcoefNorms(false);

subplot(2,2,3);
coeff_first=w3_first.ts.make_like();
coeff_first.replace_signal(w3_first.C2vec);
coeff_first.graph_sparsityDefect(false,params.c);

subplot(2,2,4);
coeff_last=w3.ts.make_like();
coeff_last.replace_signal(w3.C2vec);
coeff_last.graph_sparsityDefect(false,params.c);

end



