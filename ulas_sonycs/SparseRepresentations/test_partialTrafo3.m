%% ----- create 3d signal:
                        
        if ~exist('L','var')
            L=3;
        end
        if ~exist('fn','var')
            fn='tennis.avi'; %fn='riksch1.avi';
        end
        signal=Signal3D.make_fromVideo(fn,L);
        if signal.min<0.5
            signal.shiftSignalValues(-1);
        end
        j=1;
        s2=Signal2D(signal.xn(:,:,j)); s2.signalname=signal.signalname;
        s2.colormap_active='gray';
        
        % special case: stationary video
        % signal=signal.frame(ceil(L/2)).make_VideofromSingleImage(L);
        

%% CASE 0: comparision with single reconstruction:
% -------------------------------------------------------

        % compare with 2d- reconstuction of individual frames:
        Phi2=Wavelet2D_mlab(s2); 
        % Phi2=Curvelet2_clab(s2);
        B2=HartleyTrafo2(s2);
        B3=Wavelet2D_mlab(s2); B3.set_deepestlev(1);
        
        if ~exist('c','var')
            c=8; 
        end
        c0=c;
        %RL=60;nodes=B.nodesOnRadialLines(RL);
        params=struct; params.w1=0.995;
        % params.frameNumber=L;
        nodes2=B2.nodesDefault(c0(1),params);
        nodes3=B3.nodesDefault(c0(1),params);
         
        cs0=CompressiveSensing();  cs1=CompressiveSensing();
        cs0.set_compression(c0);   cs1.set_compression(c0);
        cs0.set_original(s2);      cs1.set_original(s2);
        % setting epsilon >0 ---> improves rescon quality !!!       
        % cs.paramsRecon.epsilon= 0.01*sqrt(signal.numel);  
        cs0.fig=2;
        cs1.fig=2;
        
        %cs0.set_methodActive(4);  % DR:1, FISTA: 2
        
        cs0.sim_SampledONS(Phi2,B2,nodes2);
        cs1.sim_SampledONS(Phi2,B3,nodes3);
        
        % compare quality results of 2 encoders:
        prepfigure(1);
        subplot(1,2,1);
        cs0.xorig.PSNRquantiles(cs0.x).graph_signal(false);
        hold all;
        cs1.xorig.PSNRquantiles(cs1.x).graph_signal(false);
        hold off;   
        legend([B2.algo.name, ', SSIM=',num2str(cs0.solutionProps.quality.SSIM,'%3.2f')],...
             [B3.algo.name, ', SSIM=',num2str(cs1.solutionProps.quality.SSIM,'%3.2f')],...
             'location','best');               
        subplot(1,2,2);
        cs0.xorig.PSNRquantiles(cs0.x).graph_signal(false);
        hold all;
        cs1.xorig.PSNRquantiles(cs1.x).graph_signal(false);
        hold off;   
        title('partial enlargement','fontsize',12);
        ylim([0,50]); xlim([0.75,1]);
        legend(B2.algo.name,B3.algo.name,'location','best');  
        suptitle('Quality of 2 encoders',14');
        
        %Phi2.dec;             
        %Phi2.graph_sparsityDefect(true,cs0.s1w_Donoho_Asymptotics);

        %% CASE 1: preliminary reconstruction of each frame separately:
% -------------------------------------------------------------------------
        xprelim=Signal3D();
        xprelim.signalname=['preliminary recon of ',signal.signalname];
        xprelim.xn=zeros(signal.size);
        multiWaitbar('Close All');
        multiWaitbar('preliminary recon',0);
        for j=1:L
            
           multiWaitbar('preliminary recon',j/L);
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
        
            cs0.fig=3;        
            cs0.sim_SampledONS(Phi,B,nodes);     
            xprelim.xn(:,:,j)=cs0.x.xn;
        end
        multiWaitbar('preliminary recon','close');
        multiWaitbar('Close All');
        

%% CASE 5: 3-dimensional MRA
% --------------------------------------
               
        
        Phi=Wavelet3D_mlab(signal);       
        
        B=HartleyTrafo3(signal);
       
        c=8; % RL=60;nodes=B.ftrafo.nodesOnRadialLines(RL,c);
        nodes=B.nodesDefault(c);
        cs=CompressiveSensing();
        cs.set_original(signal);
        % setting epsilon >0 ---> improves rescon quality !!!       
        % cs.paramsRecon.epsilon= 0.01*sqrt(signal.numel);  
        cs.sim_SampledONS(Phi,B,nodes);                      
        
        
%% CASE 6: frame pair: partial 3d, total variation
% ---------------------------------------------------
               
        Phi1=MultiResTransform3_Partial(signal,Wavelet2D_mlab());  
        %Phi1=MultiResTransform3_Partial(signal,Curvelet2_clab());  
        Phi2=TV3_Temporal(signal);
        Phi=FrameTrafoPair(signal);
        Phi.set_transforms(Phi1,Phi2);
        Phi.set_synweight1(1);
        % --------------- NEW IDEA: modifying analysis weight
        anweight1=0.5; % 0.25;
        Phi.set_anweight1(anweight1);
        
        B=FrameTrafo3_Partial(signal,HartleyTrafo2());
        params=struct;
        
%         B=FrameTrafo3_Partial(signal,Wavelet2D_mlab());
%         B.ftrafo.set_deepestlev(1); 
%         params.w1=0.995;
        
        c=8; %RL=60;nodes=B.nodesOnRadialLines(RL);
        nodes=B.nodesDefault(c,params);
        cs=CompressiveSensing();
        cs.set_original(signal);
        % setting epsilon >0 ---> improves rescon quality !!!       
        % cs.paramsRecon.epsilon= 0.01*sqrt(signal.numel);  
        % cs.paramsRecon.epsilon=0;
        %  cs.set_APhiNodes(Phi,B,nodes);   
        cs.sim_SampledONS(Phi,B,nodes);
        
        Phi.dec;
        Phi.graph_sparsityDefect(true,33.3);
        
        
%% CASE 7:  frame composition: partial 3d, total variation
% ------------------------------------------------------
          
        params=struct;
        % best results for Phi= TV_IC_MRA3 and variable compression rates:
        Phi= TV_IC_MRA3(signal);   % improves when c varies over frames !!!
        if ~exist('use_motion','var')
            use_motion=false; %use_motion=true;
        end
        if ~exist('mfield','var')
            %mfield=MotionVF(); 
            mfield=MotionQuotientField();
        end
        if use_motion
            if isa(mfield,'MotionQuotientField')
                if signal.min<0.5
                    signal.shiftSignalValues(-1);
                end
            end
            if exist('xprelim','var')
                mfield=mfield.make_fromSignal3(xprelim);
            else
                mfield=mfield.make_fromSignal3(signal);
            end
            Phi.motionCF=mfield;
        end
        Phi.use_motion=use_motion;
        if exist('zref','var')
            Phi.zref=zref;
        end
        disp(['use_motion=',num2str(Phi.use_motion),', ', class(mfield)]);
        
        %Phi= MRA_o_TV3(signal);   % deteriorates when c varies over frames
        
        % -- Curvelets give further improvement BUT take longer
        % Phi.set_transform(Curvelet2_clab());
        % Phi.anweight1=1-1/signal.size(3);
        
        % -- fourier ENCODER:        
        B=FrameTrafo3_Partial(signal,HartleyTrafo2());
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
        
        % case 1b: fixed c, but for each frame its own nodes:
        % improves reconstruction quality compared to using the same nodes 
        % for each frame!
        c=8*ones(1,L);
        
        %RL=60;nodes=B.nodesOnRadialLines(RL);               
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
            fig=2;
        end
        cs.fig=fig;
        % setting epsilon >0 ---> improves rescon quality for single image but NOT for multiple frames  
        % cs.params.epsilon= 0.01*sqrt(signal.numel);
        
        % simulate measurements with compression c and reconstruct:
        cs.sim_SampledONS(Phi,B,nodes);

        %%
        cs.x.graph_totalVar;
        cs.x.graph_motionTV;
        % decompose reconstructed signal x:
        cs.x2C;
        cs.Phi.show_trafo;
        
 
        
%% --- do some additional analysis:
        
        prepfigure(1,[],cs.figopt);
        subplot(2,2,1);
        % check cross similarity between frame pairs
        % 1.) does the reconstruction incorporate temporal changes? 
        cs.x.graph_PSNRcross([],false);
        % 2.) what temporal changes occur in original video stream?
        subplot(2,2,2);
        cs.xorig.graph_PSNRcross([],false);
        % 3.) do the reconstructed frames match best with their synchronous original
        % frame?
        subplot(2,2,3);
        cs.xorig.graph_PSNRcross(cs.x,false);
        
        prepfigure(1,[],cs.figopt);
        subplot(2,2,1);
        % check cross similarity between frame pairs
        % 1.) does the reconstruction incorporate temporal changes? 
        cs.x.graph_SSIMcross([],false);
        % 2.) what temporal changes occur in original video stream?
        subplot(2,2,2);
        cs.xorig.graph_SSIMcross([],false);
        % 3.) do the reconstructed frames match best with their synchronous original
        % frame?
        subplot(2,2,3);
        cs.xorig.graph_SSIMcross(cs.x,false);


        % analyze sparsification achieved by Phi:
        Phi.set_signal(signal);
        Phi.dec;
        if ~exist('cs','var')
            cs=CompressiveSensing();
            cs.set_original(signal);
            if exist('nodes','var')
                cs.set_APhiNodes(Phi,B,nodes);           
            end
        end
        csrate=cs.s1w_Donoho_Asymptotics;
        Phi.graph_sparsityDefect(true,csrate);   
        
        Phi.show_trafo;
        
        thresh=Phi.computeThreshold(csrate);
        Phi.hardthresh1_inplace(thresh);   
        Phi.fig=max(2,Phi.fig);
        Phi.show_trafo;

        
        % check TV of reconstruction and compare to TV of original
        f= TV3_Temporal(signal);
        f.dec;
        f.graph_trafo;
        
        f= TV3_Temporal(cs.x);
        f.dec;
        f.fig=3; f.graph_trafo;
        
        % increasing snweight should decrease TV of video:
        Phi.anweight1=0.1;
        cs.sim_SampledONS(Phi,B,nodes);          
        
        f= TV3_Temporal(cs.x);
        f.dec;
        f.fig=3; f.graph_trafo;
        
