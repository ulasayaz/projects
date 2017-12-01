function [PSN2d_mean, params]=test_2dsparsification(signal1, params)
% ~(signal1,[params]) compare 2d sparsification with 1d-sparsification
% for simulated image with global motion.
% the 2d-signal created from singal1 is:
% units=20; signal2=Signal2D.make_CyclefromLowerDimSig(signal1,signal1.length,units);
%
% INPUT:
% signal1 ... 1d signal of type TimeSignal
% params.c...  compression rate for hardthresholding
% params.videoflag ... 1 is rotation, 2 is shift
% params.L  ... extension in 3d dimension
% 
% example
% -----------
% --- define 1-dimensional signal:
% hN=256; signal1=TimeSignal.make_sinus(3,[],hN);
% signal1=TimeSignal.make_triang();
% params.c=33;
% --- choose sparsifying transform:
% params.w2= Wavelet2D_mlab(); params.w2= Wavelet2D_wlab();
% params.w2= Curvelet2_clab(); params.w2=Wavelet2WP_mlab();
% --- call function:
% test_2dsparsification(signal1, params);
%


requireOnly(isa(signal1,'TimeSignal'),'local','signal1 is a 1d signal of type TimeSignal');

use_wavelab=isa(params.w2,'Wavelet2D_wlab');
use_squares =use_wavelab ;

if ~exist('params','var')
    params=struct;
end

if ~isfield(params,'c')
    params.c=33;  % compression rate
end
if ~isfield(params,'L')
    params.L=20;   % minimum for Curvelet2_clab is 17
end

if ~isfield(params,'basisname2d')
    params.basisname2d='db1';
end
if ~isfield(params,'fig')
    params.fig=4;
end

if ~isfield(params,'anisotropic')
    % if true, choose deepest possible decomposition level of largest
    % dimensional extension;
    % if false, choose deepest possible level of smallest dimensional
    % extension.
    params.anisotropic=true;
end



if use_squares
    params.L=signal1.length;
end
L=params.L;


opt.pixsizeX=1000;


units=1;

f=[0.1,0.5,1,1.5,3,5,10,20];
set_units=units*f;
K=length(set_units);

PSN2d_mean=zeros(1,K);

% loop over several magnitudes of motion
wait_handle = waitbar(0,...
                ['loop over magnitudes of motion (K=',num2str(K),') ...']);
            
for k=1:K
    
    waitbar(k / K);
    units=set_units(k);
    
    imageflag_descr=['global shift (units',num2str(units,'%3.1f'),')'];
    imageflag_type='global shift (pix)';
    
    % --- 2d sparsification:
    signal2=Signal2D.make_CyclefromLowerDimSig(signal1,L,units);
    
    % prepfigure(1,[],opt); signal2.graph_signal([],false);
    % suptitle([signal2.get_signalname,': shift_sequence'],14);
    wclass=class(params.w2);
    constructor=str2func(wclass);
    w2= constructor();
    
    if ~isa(w2,'Curvelet2_clab');
        w2.set_basisname(params.basisname2d);
    end
    w2.set_signal(signal2);
    if  params.anisotropic
        lev=w2.deepestlevAnisotropic;
    else
        lev=w2.deepestlev;
    end
    
    %lev=w2.deepestlevAnisotropic;
    w2.dec(lev);
    if params.anisotropic && w2.level~=lev
        params.anisotropic=false;
    end
    
    if k==1
        signal2_first=signal2;
        w2_first=w2;
    end
    
    test2=w2.sparseApprox(params.c);
    if K==1
        test2.graph_signal;  %2d image reconstructed
        % test2.graph_signal(1); title(imageflag_descr,'fontsize',12); % first reconstructed image of image
        
    end
    
    PSNR2d=signal2.PSNR(test2);
    PSN2d_mean(k)=mean(PSNR2d);
    
    shiftpix=[0,signal2.shiftdist(L,units)];
    
    xvals=(sum(shiftpix.^2,1)).^0.5;
    xdescr='shift [pix]';
        
end
close(wait_handle );

% ---- compare with 2D sparsification:

if use_wavelab
    w1= Wavelet1D_wlab(signal1);
else
    w1= Wavelet1D(signal1);
end

w1.set_basisname(params.basisname2d);
w1.dec;
test1=w1.sparseApprox(params.c);
test1.signalname=[signal1.get_signalname,' sparse Approx. c=',num2str(params.c)];
test1.colormap_active='gray';
if K==1
    test1.graph_signal; % 1d recon
end


PSNR1d=signal1.PSNR(test1);


%% results: first figure
%=======================================================================
if params.anisotropic
    isotropy_str='anisotropic dec.';
else
    isotropy_str='isotropic dec.';
end

prepfigure(params.fig,[],opt);

suptitle([signal1.signalname,', size=',vec2str(signal2.size),', ',...
        isotropy_str], 14);
    
subplot(2,4,1:2);
signal1.graph_signal(false);

subplot(2,4,5);
signal2_first.graph_signal(false);
subplot(2,4,6);
try
    w2_first.graph_trafo(false);
catch
end
cbfreeze;
title(w2_first.transformName);
if K>1
    subplot(2,4,7);
    signal2.graph_signal(false);
    subplot(2,4,8);
    try
        w2.graph_trafo(false);
    catch
    end
    title(w2.transformName);
end

subplot(2,4,3:4);
if K==1
    plot(xvals,signal2.PSNR(test2));
    title(['2d sparsification PSNR  (PSNR(1d)= ',num2str(PSNR1d,'%4.2g'),')'],'fontsize',12);
    xlabel(xdescr,'fontsize',10);
    ylabel('PSNR');
    
else
    norm_units=zeros(1,K);
    for k=1:K
        norm_units(k)=norm( set_units(:,k));
    end
    
    w2class=strrep(class(w2),'_','\_');
    w1class=strrep(class(w1),'_','\_');
    
    semilogx(norm_units, PSN2d_mean);
    line([norm_units(1), norm_units(end)],[PSNR1d,PSNR1d],'color','red');
    legend(['2d: mean PSNR of ',num2str(L),' frames'],['1d: ',w1class,'(',...
        w1.basisname,')'],...
        'location','best');
    xlabel(imageflag_type,'fontsize',10);
    ylabel('PSNR');
    
end

title(['PSNR (',w2class,'(',w2.basisname,...
    '), compression c=',num2str(params.c,'%3.0f'),')'],'fontsize',12);

%% second figure
%=======================================================================

prepfigure(params.fig,[],opt);

subplot(2,2,1);
w2_first.graph_detcoefNorms(false);

subplot(2,2,2);
w2.graph_detcoefNorms(false);

subplot(2,2,3);
coeff_first=w2_first.ts.make_like();
coeff_first.replace_signal(w2_first.C2vec);
coeff_first.graph_sparsityDefect(false,params.c);

subplot(2,2,4);
coeff_last=w2.ts.make_like();
coeff_last.replace_signal(w2.C2vec);
coeff_last.graph_sparsityDefect(false,params.c);

end




