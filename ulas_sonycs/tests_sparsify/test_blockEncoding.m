%signal=Signal3D.make_CyclefromLowerDimSig(s2,videoflag,L,motionunit);
L=3; 

%fn='tennis.avi'; 
fn='riksch1.avi';
signal=Signal3D.make_fromVideo(fn,L);

block = [32,64,128,256,512];

K=length(block);

PS_mean=zeros(1,K);
SS_mean=zeros(1,K);

for k=1:K
    k
    blocksize = block(k);
    
    params=struct;
    % best results for Phi= TV_IC_MRA3 and variable compression rates:
    Phi= TV_IC_MRA3(signal);   % improves when c varies over frames !!!
    %Phi= MRA_o_TV3(signal);   % deteriorates when c varies over frames
    
    % -- Curvelets give further improvement BUT take longer
    % Phi.set_transform(Curvelet2_clab());
    % Phi.anweight1=1-1/signal.size(3);
    
    % -- block fourier ENCODER:
    if k == K
        B=FrameTrafo3_Partial(signal,HartleyTrafo2());
    else
        B=FrameTrafo3_Partial_Blocks_U(signal,HartleyTrafo2(),blocksize);
    end
    
    % -- wavelet ENCODER: (do not choose if Phi uses wavelet, too;
    % (mutual coherence increases otherwise).
    %         B=FrameTrafo3_Partial(signal,Wavelet2D_mlab());
    %         B.ftrafo.set_deepestlev(1);
    %         params.w1=0.995;
    
    % -- case 2: compression rate varying over frames:
    % params.k0 ... weight of side frames compared to central frame
    % 1/params.p0 ,,, mean compression rate
    % c=zeros(1,L); c0=8; params.p0=1/c0; params.k0=1/2;  % e.g. c0=16; k0=1/16 or 1/8
    
    % -- case 1: compression rate fixed:
    % c=8;
    
    % case 1b: fixed c, but for each frame its own nodes,
    % same nodes for the blocks in a frame:
    % improves reconstruction quality compared to using the same nodes
    % for each frame!
    c=8*ones(1,L);
    
    %RL=60;nodes=B.nodesOnRadialLines(RL);
    % params.frameNumber=L;  % define for each frame new mask
    % OBSERVATION: quality gets worse for varying mask
    [nodes,c,cmean]=B.nodesDefault(c,params);
    
    % check, if nodes are ok
    assert(xor(iscell(nodes.data),nodes.fL==L),'#frames is set in nodes');
    
    % create compressive-sensing object:
    CS=CompressiveSensing_U(); % uses SampleTrafo_U
    CS.set_original(signal);
    CS.set_methodActive(5);  % DR:1, FISTA: 2
    disp(CS.get_methodActive);
    % setting epsilon >0 ---> improves rescon quality for single image
    % but NOT for multiple frames CS.params.epsilon=
    % 0.01*sqrt(signal.numel);
    
    % simulate measurements with compression c and reconstruct:
    % parameter setting CS.params.epsilon=0.1; PD
    %if L>=9
        % reduce parameters critical for convergence
        %CS.params.PD.sigma=0.05;
        %CS.params.PD.tau=0.05;
        CS.params.PD.maxIters=500;
    %end
    
    CS.fig=3;
    
    CS.sim_SampledONS(Phi,B,nodes);
    
    PS_mean(k) = mean(CS.solutionProps.quality.PSNR);
    SS_mean(k) = mean(CS.solutionProps.quality.SSIM);
    
end

% plot results

N = signal.size(1);

figure
plot(block,PS_mean(:),'LineWidth',2,'Marker','s');
title([fn,', L=',int2str(L),', PSNR vs block size, resolution =', int2str(N)],...
    'FontSize',12);
xlabel('blocksize','FontSize',12);
ylabel('PSNR','Fontsize',12);

figure
plot(block,SS_mean(:),'Color','r','LineWidth',2,'Marker','s');
title([fn,', L=',int2str(L),', SSIM vs block size, resolution =',...
    int2str(N)],'FontSize',12)
xlabel('blocksize','FontSize',12);
ylabel('SSIM','Fontsize',12);




