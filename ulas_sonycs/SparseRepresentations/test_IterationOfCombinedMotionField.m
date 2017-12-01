% test_IterationOfCombinedMotionField.m

% STEP 0:  initialization
tic;
disp(datestr(now));
p=1;  % lp-norm
if ~exist('L','var')
    L=3;
end
if ~exist('do_reduce_size','var') 
   do_reduce_size=false;
end
if ~exist('signal','var') || isempty(signal) || signal.size(3)~=L
    if ~exist('fn','var')
        fn='riksch1.avi'; %fn='tennis.avi';
    end
    signal=Signal3D.make_fromVideo(fn,L);
end
if do_reduce_size
    signal.resize(min(64,signal.size(1)));
end
if signal.min<0.5
    signal.shiftSignalValues(-1);
end
j=1;
s2=Signal2D(signal.xn(:,:,j)); s2.signalname=signal.signalname;
s2.colormap_active='gray';
toc;

% STEP 1: simulate measurement and preliminary reconstruction
% of each frame separately.
if ~exist('xprelim','var') || ~isequal(xprelim.size,signal.size) ...
        || ~exist('PSNR_prelim','var') || ~exist('SSIM_prelim','var')
    tic;
    [ xprelim,cs0,PSNR_prelim,SSIM_prelim ] = CS_SimPremimRecon( signal );
    toc;    
end

tic;
if ~exist('iter','var')
    iter=2;
end

PSNR=zeros(L,iter+1);
SSIM=PSNR;
PSNR(:,1)=reshape(PSNR_prelim(1:L),[],1);
SSIM(:,1)=reshape(SSIM_prelim(1:L),[],1);

multiWaitbar('Close All');
multiWaitbar('recon-motion-field-iteration',1/(iter+1));
xbest=xprelim.clone;
xset=cell(1,iter);

vf=cell(1,iter);

params.verbose=false;
params.fig=1;
for j=2:iter+1
    % STEP2: motion estimation
    
    sorig=signal;
    opts=struct;
    %-- choose method (OF is best):
    opts.method_estim='OF'; % method_estim='FSBM';
    vf{j-1}=MotionVF.make_fromSignal3(xbest,opts);
    % -- treat occlusions by using the quotient field:
    vftemp=vf{j-1};
    vftemp.set_weightFieldFromSignal3(xbest);
    vftemp.motionfname=[vftemp.motionfname,', it=',num2str(j)];
        
    % STEP4 :  frame composition: partial 3d, total variation
        
    for k=1:L
        params.zref=k;
        [ x, cs ] = CS_MotionCompRecon( sorig, vftemp,params );       
        
        PSNR(k,j)=cs.solutionProps.quality.PSNR(cs.Phi.zref);
        SSIM(k,j)=cs.solutionProps.quality.SSIM(cs.Phi.zref);
        disp(['[j,k]=',vec2str([j,k]),': PSNR=',num2str(PSNR(k,j))]);
        
        xbest.xn(:,:,k)=x.frame(cs.Phi.zref).xn;
        multiWaitbar('recon-motion-field-iteration',j*k/(L*(iter+1)));
    end
        
    xbest.signalname=strrep(xbest.signalname,'preliminary',['it=',num2str(j)]);
    xbest.signalname=strrep(xbest.signalname,['it=',num2str(j-1)],['it=',num2str(j)]);
    
    xset{j-1}=xbest.clone;    
end


multiWaitbar('recon-motion-field-iteration','close');
multiWaitbar('Close All');

disp(datestr(now));
toc;


