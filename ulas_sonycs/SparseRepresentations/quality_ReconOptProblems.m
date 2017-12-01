% quality_ReconOptProblems
% determines quality of reconstruction of a video by varying
% the number of frames used.
% 2 cases implemented: either with or without motion compensation
% where motion compensation uses the original signal instead of
% a preliminary reconstruction.

if ~exist('fn','var')
    % fn='tennis.avi'; % tennis.avi works very well without motion compensation
    fn='riksch1.avi';  % risksch1.avi needs motion compensation
end
if ~exist('framecount','var')
    a=7;
    framecount=2*linspace(0,a,a+1)+1; % framecount=[1,3];
end
if ~exist('compensate_motion','var')
    compensate_motion=true;
end
if ~exist('reduce_size','var')
    reduce_size=false;  % reduces image size in order to accelerate computation
end
if ~exist('max_size','var')
    max_size=64;
end
if ~exist('show_error','var')
    show_error=false;
end
Lfc= length(framecount);

PSNR_reached=cell(1,Lfc);
SSIM_reached=cell(1,Lfc);
SSIM_opt_vals=cell(1,Lfc);
PSNR_opt_vals=cell(1,Lfc);

%% motion estimation
if ~exist('vf','var')
    vf=[];
end
L=max(framecount);
video=Signal3D.make_fromVideo(fn,L);
if video.size(3)<L
   video.repSignal(ceil(L/video.size(3)));
   video=video.frame(1:L);
end
framecount=framecount(framecount<=L);
if reduce_size
    video.resize(max_size);
end
if video.min<0.5
    video.shiftSignalValues(-1+video.min);
end
assert(video.min>0.5,'signal bounded away from 0');
video_size=video.size;

if compensate_motion
    
    if isempty(vf) || vf.size(3)+1<max(framecount) ...
            || ~isequal(vf.size_space,video_size(1:2)) ...
            || ~strcmp(vf.signal3.signalname,video.signalname)
        display('running with motion estimation ...');
        display('(using the original signal (ideal non-realistic case))');
        
        if ~exist('opts','var') || ~isstruct(opts)
            opts=struct;
        end
        %-- choose method (OF is best):
        opts.method_estim='OF'; % method_estim='FSBM';
        vf=MotionVF.make_fromSignal3(video,opts);
        % -- treat occlusions by using the quotient field:
        vf.set_weightFieldFromSignal3(vf.signal3);
    else
        display('using precalculated motion estimation vf');
    end
else
    display('WITHOUT motion estimation !');
end

multiWaitbar('Close All');
multiWaitbar('loop over frame count',0);
disp('PSNR: ');
%% loop over frame counts
for j=1:Lfc
    
    L=framecount(j);
    if L>1
        signal=video.frame(1:L);
        Phi= TV_IC_MRA3(signal);   % improves when c varies over frames !!!
        if compensate_motion
            Phi.motionCF=vf.frame(1:L-1);
            Phi.use_motion=true;
            Phi.set_algo;
        else
            Phi.use_motion=false;
            Phi.set_algo;
        end
    else
        signal=video.frame(1);
        Phi=Wavelet2D_mlab(signal);
    end
    
    
    params=struct;
    params.w1=0.995;
    % best results for Phi= TV_IC_MRA3 and variable compression rates:
    
    
    % -- Curvelets give further improvement BUT take longer
    % Phi.set_transform(Curvelet2_clab());
    % Phi.anweight1=1-1/signal.size(3);
    
    % -- fourier ENCODER:
    if L>1
        B=FrameTrafo3_Partial(signal,HartleyTrafo2());
        % -- wavelet ENCODER:
        %   B=FrameTrafo3_Partial(signal,Wavelet2D_mlab());
        %   B.ftrafo.set_deepestlev(1);
    else
        B=HartleyTrafo2(signal);
        % -- wavelet ENCODER:
        %   B=Wavelet2D_mlab(signal);
        %   B.set_deepestlev(1);
    end
    
    % -- case 2: compression rate varying over frames:
    % params.k0 ... weight of side frames compared to central frame
    % 1/params.p0 ,,, mean compression rate
    % c=zeros(1,L); c0=8; params.p0=1/c0; params.k0=1/2;  % e.g. c0=16; k0=1/16 or 1/8
    
    % -- case 1: compression rate fixed:
    % c=8;
    
    % case 1b: fixed c, but for each frame its own nodes:
    % reduces recon quality of problems 2,3 but not 1. Why?
    c=8*ones(1,L);
    
    %RL=60;nodes=B.nodesOnRadialLines(RL);
    % params.frameNumber=L;  % define for each frame new mask
    % OBSERVATION: quality gets worse for varying mask
    nodes=B.nodesDefault(c,params);
    
    % check, if nodes are ok
    assert(xor(iscell(nodes.data),nodes.fL==L),'#frames is set in nodes');
    
    % create compressive-sensing object:
    cs=CompressiveSensing();
    cs.set_original(signal);
    cs.set_methodActive(4);  % DR:1, FISTA: 2
    if L>=9
        % reduce parameters critical for convergence
        cs.params.PD.sigma=cs.params.PD.sigma/2;
        cs.params.PD.tau=cs.params.PD.tau/2;
        cs.params.PD.maxIters=500;
    end
    
    % simulate measurements with compression c and reconstruct:
    cs.fig=0;
    cs.sim_SampledONS(Phi,B,nodes);
    
    PSNR_reached{j}=  cs.solutionProps.quality.PSNR;
    SSIM_reached{j}= cs.solutionProps.quality.SSIM;   
    
    if isempty(cs.xopt)
        cs.compute_QualityBound();
    end
    SSIM_opt_vals{j}=cs.xorig.SSIM(cs.xopt);
    PSNR_opt_vals{j}=cs.xorig.PSNR(cs.xopt);
        
    
    fprintf([' --*',num2str(L),': ',num2str(max(PSNR_reached{j}))]);        
    multiWaitbar('loop over frame count',j/Lfc);
end
disp(' ');
multiWaitbar('loop over frame count','close');

PSNR_reached_mean=cellfun(@mean,PSNR_reached);
PSNR_reached_max=cellfun(@max,PSNR_reached);
PSNR_reached_min=cellfun(@min,PSNR_reached);
PSNR_reached_std=cellfun(@std,PSNR_reached);

PSNR_bound=cellfun(@max,PSNR_opt_vals);
PSNR_bound_std=cellfun(@std,PSNR_opt_vals);

SSIM_reached_mean=cellfun(@mean,SSIM_reached);
SSIM_reached_max=cellfun(@max,SSIM_reached);
SSIM_reached_min=cellfun(@min,SSIM_reached);
SSIM_reached_std=cellfun(@std,SSIM_reached);

SSIM_bound=cellfun(@max,SSIM_opt_vals);

ppfig.pixsizeX=1000;
prepfigure(1,[],ppfig);

subplot(1,2,1);
opt.outputY=true; opt.format='%3.1f';
if show_error
    errorbar(framecount,PSNR_reached_mean,PSNR_reached_std);
    hold all;
    errorbar(framecount,PSNR_bound,PSNR_bound_std);
    hold off;
    textvec(framecount,PSNR_reached_mean,opt);
else
    plot(framecount,PSNR_reached_max);
    hold all;
    plot(framecount,PSNR_bound);
    hold off;
    textvec(framecount,PSNR_reached_max,opt);
end

legend('PD_{l1}','DT-quality bound','location','best');
xlabel('video length L in frames','fontsize',12);
ylabel('PSNR [dB]','fontsize',12);

subplot(1,2,2);
opt.format='%3.2f';
if show_error
    errorbar(framecount,SSIM_reached_mean,SSIM_reached_std);
    hold all;
    plot(framecount,SSIM_bound);
    hold off;
    textvec(framecount,SSIM_reached_mean,opt);
else
    plot(framecount,SSIM_reached_max);
    hold all;
    plot(framecount,SSIM_bound);
    hold off;
    textvec(framecount,SSIM_reached_max,opt);
end

legend('PD_{l1}','DT-quality bound','location','best');
xlabel('video length L in frames','fontsize',12);
ylabel('SSIM','fontsize',12);

encoder=B.ftrafo.algo.name;
if length(c)>1
    encoder=[encoder,', new nodes for each frame'];
end
meanc=1/mean(1./c);

suptitle({[cs.x.signalname,' (',cs.solver_paramstr,...
    ', it=', num2str(cs.solutionProps.numIters),')'...
    ],...
    ['encoder ',encoder,', mean(c)=',num2str(meanc,'%3.1f')],...
    ['decoder ',cs.Phi.algo.name,' based on orig.']},14);




