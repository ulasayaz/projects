% test_combinedMotionField.m

%% initialization
tic;
p=1;  % lp-norm
if ~exist('L','var')
    L=3;
end
if ~exist('fn','var')
    fn='riksch1.avi'; %fn='tennis.avi'; 
end
if ~exist('use_prelim','var')
    use_prelim=false; % use_prelim=true;
end
if ~exist('use_denoised','var')
    use_denoised=false;
end
signal=Signal3D.make_fromVideo(fn,L);
if signal.min<0.5
    signal.shiftSignalValues(-1);
end
j=1;
s2=Signal2D(signal.xn(:,:,j)); s2.signalname=signal.signalname;
s2.colormap_active='gray';
toc;
%% STEP 0: test decoder:
        tic;
    % choose decoder:
    w= TV_IC_MRA3(signal);
    % --- test if analyze and sythesize are inverses of each other:        
        [ok,err1,err2]= w.test_AnalyzeSynthsize(); 
        disp(['Inversion Test 1 Passed=',num2str(ok), ' with error ',num2str(max(err1,err2))]);
        [ok,err1,err2]=w.test_DecRec(); 
        disp(['Inversion Test 2 Passed=',num2str(ok), ' with error ',num2str(max(err1,err2))]);
       
    % --- test if adjointSynthesis is the adjoint of synthesize    
        disp('computing the motion fields takes longer, please wait ...'); 
        [ok,err]= TV_IC_MRA3.test_adjointSynthesis();
        disp(['Adjointness Test 1 passed=',num2str(ok), ' with error ',num2str(err)]);
        toc;
 %% STEP 1: preliminary reconstruction of each frame separately:
% -------------------------------------------------------------------------
        % load('Rikscha15_prelim.mat');
        tic;
        display('running preliminary frame by frame reconstruction ...');
        xprelim=signal.make_like;
        
        xprelim.signalname=['preliminary recon of ',signal.signalname];
        xprelim.xn=zeros(signal.size);
        PSNR_prelim=zeros(1,L);
        SSIM_prelim=zeros(1,L);
        multiWaitbar('Close All');
        multiWaitbar('preliminary recon',0);
        for j=1:L
                       
           sj=signal.frame(j);
           Phi=Wavelet2D_mlab(sj); 
           % Phi2=Curvelet2_clab(sj);
           B=HartleyTrafo2(sj);
       
            if ~exist('c','var')
                c=8; 
            end
            c0=c;        
            params=struct; params.w1=0.995;
            % params.frameNumber=L;
            nodes=B.nodesDefault(c0(1),params);
                 
            cs0=CompressiveSensing(); 
            cs0.set_compression(c0);
            cs0.set_original(sj);
        
            cs0.fig=0;        
            cs0.sim_SampledONS(Phi,B,nodes);   
            PSNR_prelim(j)=cs0.solutionProps.quality.PSNR;
            SSIM_prelim(j)=cs0.solutionProps.quality.SSIM;
            disp(['j=',num2str(j),', PSNR=',num2str(PSNR_prelim(j),'%3.1f')]);
            xprelim.xn(:,:,j)=cs0.x.xn;
            
            multiWaitbar('preliminary recon',j/L);
        end
        multiWaitbar('preliminary recon','close');
        multiWaitbar('Close All');        
        toc;
        
         %% STEP1B denoising preliminary recon
        tic;        
        if use_denoised
             opts.thresh1=0.05;
             [xpDenoised,opts]=xprelim.denoise(Wavelet2D_mlab(),opts);
        end
        toc;
        
        %% STEP2: motion estimation  
        tic;
        
        display('running motion estimation takes longer, please wait ...');
        if use_denoised
            display('(using the denoised reconstruction)');
        elseif use_prelim 
            display('(using the preliminary reconstruction)');
        else
            display('(using the original signal (ideal non-realistic case))');
        end
        if ~exist('opts','var') || ~isstruct(opts)
            opts=struct; 
        end
        %-- choose method (OF is best):
        opts.method_estim='OF'; % method_estim='FSBM';

        if use_denoised
            vf=MotionVF.make_fromSignal3(xpDenoised,opts);
        elseif use_prelim
            vf=MotionVF.make_fromSignal3(xprelim,opts);
        else
            vf=MotionVF.make_fromSignal3(signal,opts);
        end
        
        vf.fig=3;
        vf.show_NormsAnglesDiff(1);
%         vf.show_field_TV();
%         vf.show_norm;    
%         vf.show_occlusions();
%         vf.show_crossings();
        
         % -- treat occlusions by keeping moving pixels
        %vf.show_motionCompensation();
        if use_prelim || use_denoised
            %vf.set_weightField2scalar(1);
            if ~exist('srel','var')
                srel=0.1;
            end
            deltaQ1_PonP=vf.sparsityDefect_delta(srel);
        end
        %vf.test_motionCompensation(1);
        %vf.test_motionCompensation();
        
    
    % -- treat occlusions by using the quotient field:
        vf.set_weightFieldFromSignal3(vf.signal3);
        %vf.show_motionCompensation();
        %vf.test_motionCompensation(1);
        %vf.test_motionCompensation();
        
        occ=vf.occlusions;
        cross=vf.crossings;
        occ_area=zeros(1,size(occ,3));
        cross_area=occ_area;
        L=length(occ_area);
        for j=1:L
            occ_area(j)=sum(reshape(occ(:,:,j),[],1))/numel(occ(:,:,j));
            cross_area(j)=sum(reshape(cross(:,:,j),[],1))/numel(occ(:,:,j));
        end
        
        figure; plot(1:L,occ_area,1:L,cross_area); 
        legend('occlusions','crossings','location','best');
        grid on;
        textopt.outputY=true;
        textvec(1:L,occ_area,textopt);
        xlabel('frame pair number','fontsize',12);
        if ~(use_prelim || use_denoised)
            title(['share of occluded and crossing points: ',fn],'fontsize',12);
        else
            title(['share of occluded and crossing points: ','1. recon of ',fn],'fontsize',12);
        end
        
        
        if ~exist('srel','var')
           srel=0.1;
        end
        if ~(use_prelim || use_denoised)
            delta_OrigonOrig=vf.sparsityDefect_delta(srel);
            diff_OrigOnOrig=vf.sparsityDefect_diff(srel);
        else
            diff_PonP=vf.sparsityDefect_diff(srel);
            delta_PonP=vf.sparsityDefect_delta(srel);
        end
               
        disp(datestr(now));
        toc;
               
        
        %% STEP3: test motion estimation of prelim recon on original
        tic;
        
        L=vf.signal3.size(3)-1;
        
        if use_prelim || use_denoised
                 
            if ~exist('srel','var')
                srel=0.1;
            end
            
            vf.set_signal(xprelim);
            if ~exist('delta_PonP','var')
                diff_PonP=vf.sparsityDefect_diff(srel);
                delta_PonP=vf.sparsityDefect_delta(srel);
            end
            if ~exist('deltaQ1_PonP','var');
                deltaQ1_PonP=NaN*zeros(1,L);
            end
            
            vf.set_signal(signal); 
            
            delta_PonOrig=vf.sparsityDefect_delta(srel);
            diff_OrigOnOrig=vf.sparsityDefect_diff(srel);                            
            %vf.set_signal(xprelim);
                       
            prepfigure(1,[]);
            if exist('delta_OrigonOrig','var') && ...
                    length(delta_OrigonOrig)==L
                plot(...
                    1:L,diff_PonP, ...                    
                    1:L, deltaQ1_PonP,...
                    1:L, delta_PonP,...
                    1:L, delta_PonOrig,...
                    1:L, delta_OrigonOrig,...
                    'LineWidth',2);
                legend(...
                    'uncompensated D of 1.recon', ...
                    'compensated \Delta(V,Q=1) of 1.recon',...
                    'compensated \Delta(V,Q) of 1.recon',...
                    'compensated \Delta(V,Q) of orig using 1.recon',...
                    'compensated \Delta(V,Q) of orig using orig',...
                    'location','best');
            else
                plot(1:L,diff_PonP, ...
                    1:L, delta_PonP,...
                    1:L, deltaQ1_PonP,...
                    1:L, delta_PonOrig,...
                    'LineWidth',2);
                legend('uncompensated D of 1st recon', ...
                    'compensated \Delta(V,Q=1) of 1.recon',...      
                    'compensated \Delta(V,Q) of 1st recon',...                                 
                    'compensated \Delta(V,Q) of orig using recon',...
                    'location','best');
            end
            title({['compression capabilities of difference operators (',fn,')'],
            ['i.e. rel. sparsity defects of differences of consecutive images at c=',...
                num2str(1/srel)]},'fontsize',12);
            ylabel('\sigma_{1/c}(x)_1/dim(x)','fontsize',12);
            xlabel('frame pair number','fontsize',12);
            
            % -- treat occlusions by using the quotient field:            
%             vf.show_motionCompensation();
%             vf.test_motionCompensation(1);
%             vf.test_motionCompensation();
%             title({'use preliminary motion field on original',get_title},...
%                 'fontsize',12);
        end
        toc;
        
        %% STEP4 :  frame composition: partial 3d, total variation
% ------------------------------------------------------
        tic;  
        display('running video reconstruction using motion compensation ...');
        params=struct;
        % best results for Phi= TV_IC_MRA3 and variable compression rates:
        Phi= TV_IC_MRA3(signal);   % improves when c varies over frames !!!
        
        Phi.motionCF=vf;
        
        Phi.use_motion=true;
        Phi.set_algo;
        if exist('zref','var')
            Phi.zref=zref;
        end
        disp(['use_motion=',num2str(Phi.use_motion),', ', class(vf)]);
              
        % -- Curvelets give further improvement BUT take longer
        % Phi.set_transform(Curvelet2_clab());
        % Phi.anweight1=1-1/signal.size(3);
        
        % -- fourier ENCODER:        
        B=FrameTrafo3_Partial(signal,HartleyTrafo2());   
        L=signal.size(3);
        c=8*ones(1,L);
                          
        [nodes,c,cmean]=B.nodesDefault(c,params);        
        
        % check, if nodes are ok
        assert(xor(iscell(nodes.data),nodes.fL==L),'#frames is set in nodes');
           
        % create compressive-sensing object:
        cs=CompressiveSensing();
        cs.set_original(signal);
        cs.set_methodActive(4);  % DR:1, FISTA: 2
        disp(cs.get_methodActive);
        if L>=9
            % reduce parameters critical for convergence
            cs.params.PD.sigma=cs.params.PD.sigma/2;
            cs.params.PD.tau=cs.params.PD.tau/2;
            cs.params.PD.maxIters=500;
        end
            
        if ~exist('fig','var')
            fig=1;
        end
        cs.fig=fig;
        % setting epsilon >0 ---> improves rescon quality for single image but NOT for multiple frames  
        % cs.params.epsilon= 0.01*sqrt(signal.numel);
        
        % simulate measurements with compression c and reconstruct:
        cs.sim_SampledONS(Phi,B,nodes);
        % cs.show_recon;
        disp(datestr(now));
        toc;
