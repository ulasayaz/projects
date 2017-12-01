%% ----- create 3d signal:

        L=16;               
        %s2=Signal2D.make_fromImage('cameraman.bmp');       
        %videoflag=2; motionunit=[1.2,1.4]; %  [0,0];  % [0.1,0.1]
        %signal=Signal3D.make_CyclefromLowerDimSig(s2,videoflag,L,motionunit);
        L=3; fn='tennis.avi'; % fn='riksch1.avi';
        signal=Signal3D.make_fromVideo(fn,L);
        j=1;
        s2=Signal2D(signal.xn(:,:,j)); s2.signalname=signal.signalname;
        s2.colormap_active='gray';
        
        % special case: stationary video
        % signal=signal.frame(ceil(L/2)).make_VideofromSingleImage(L);
        

%% CASE 0: comparison with individual reconstruction:
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
         
        CS0=CompressiveSensing_U();  CS1=CompressiveSensing_U();
        CS0.set_compression(c0); CS1.set_compression(c0);
        CS0.set_original(s2); CS1.set_original(s2);
        % setting epsilon >0 ---> improves rescon quality !!!       
        % CS.paramsRecon.epsilon= 0.01*sqrt(signal.numel);  
        CS0.sim_SampledONS(Phi2,B2,nodes2);
        CS1.sim_SampledONS(Phi2,B3,nodes3);
        
        % compare quality results of 2 encoders:
        prepfigure(1);
        subplot(1,2,1);
        CS0.xorig.PSNRquantiles(CS0.x).graph_signal(false);
        hold all;
        CS1.xorig.PSNRquantiles(CS1.x).graph_signal(false);
        hold off;   
        legend([B2.algo.name, ', SSIM=',num2str(CS0.solutionProps.quality.SSIM,'%3.2f')],...
             [B3.algo.name, ', SSIM=',num2str(CS1.solutionProps.quality.SSIM,'%3.2f')],...
             'location','best');               
        subplot(1,2,2);
        CS0.xorig.PSNRquantiles(CS0.x).graph_signal(false);
        hold all;
        CS1.xorig.PSNRquantiles(CS1.x).graph_signal(false);
        hold off;   
        ylim([0,50]); xlim([0.75,1]);
        legend(B2.algo.name,B3.algo.name,'location','best');  
        suptitle('Quality of 2 encoders',14');
        
        %Phi2.dec;             
        %Phi2.graph_sparsityDefect(true,CS0.s1w_Donoho_AsymptotiCS);


%% CASE 5: 3-dimensional MRA
% --------------------------------------
               
        
        Phi=Wavelet3D_mlab(signal);       
        
        B=HartleyTrafo3(signal);
       
        c=8; % RL=60;nodes=B.ftrafo.nodesOnRadialLines(RL,c);
        nodes=B.nodesDefault(c);
        CS=CompressiveSensing_U();
        CS.set_original(signal);
        % setting epsilon >0 ---> improves rescon quality !!!       
        % CS.paramsRecon.epsilon= 0.01*sqrt(signal.numel);  
        CS.sim_SampledONS(Phi,B,nodes);                      
        
        
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
        CS=CompressiveSensing_U();
        CS.set_original(signal);
        CS.set_methodActive(4);  % DR:1, FISTA: 2
        disp(CS.get_methodActive);
        % setting epsilon >0 ---> improves rescon quality for single image
        % but NOT for multiple frames CS.params.epsilon=
        % 0.01*sqrt(signal.numel);
        
        % simulate measurements with compression c and reconstruct:
        % parameter setting CS.params.epsilon=0.1; PD
        if L>=9
            % reduce parameters critical for convergence
            CS.params.PD.sigma=0.05;
            CS.params.PD.tau=0.05;
            CS.params.PD.maxIters=500;
        end
        % setting epsilon >0 ---> improves rescon quality !!!       
        % CS.paramsRecon.epsilon= 0.01*sqrt(signal.numel);  
        % CS.paramsRecon.epsilon=0;
        %  CS.set_APhiNodes(Phi,B,nodes);   
        CS.sim_SampledONS(Phi,B,nodes);
        
        Phi.dec;
        Phi.graph_sparsityDefect(true,33.3);
        
        
%% CASE 7:  frame composition: partial 3d, total variation
% ------------------------------------------------------
          
        params=struct;
        % best results for Phi= TV_IC_MRA3 and variable compression rates:
        Phi= TV_IC_MRA3(signal);   % improves when c varies over frames !!!
        %Phi= MRA_o_TV3(signal);   % deteriorates when c varies over frames
        
        % -- Curvelets give further improvement BUT take longer
        % Phi.set_transform(Curvelet2_clab());
        % Phi.anweight1=1-1/signal.size(3);
        
        % -- fourier ENCODER: 
        
        % -- case 1: block measurements
        blocksize = 64;
        B=FrameTrafo3_Partial_Blocks_U(signal,HartleyTrafo2(),blocksize); 
        
        % -- case 2: regular measurements       
        % B=FrameTrafo3_Partial(signal,HartleyTrafo2());
           
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
        CS.set_methodActive(4);  % DR:1, FISTA: 2
        disp(CS.get_methodActive);
        % setting epsilon >0 ---> improves rescon quality for single image
        % but NOT for multiple frames CS.params.epsilon=
        % 0.01*sqrt(signal.numel);
        
        % simulate measurements with compression c and reconstruct:
        % parameter setting CS.params.epsilon=0.1; PD
        if L>=9
            % reduce parameters critical for convergence
            CS.params.PD.sigma=0.05;
            CS.params.PD.tau=0.05;
            CS.params.PD.maxIters=500;
        end
        % DR
        CS.params.DR.gamma=0.1;
        CS.params.DR.maxIters=200;
        
        CS.fig=3;
        
        CS.sim_SampledONS(Phi,B,nodes);
        
        blocksize = N;
        
        % show first frame
        figure;
        imgray(CS.x.xn(:,:,1));
        title([fn,', L=',int2str(L),', resolution =',...
            int2str(N),', blocksize =', int2str(blocksize),...
            ', first frame'],'FontSize',12);
        xlabel(['reconstructed ',fn],'FontSize',12);
        
        % show reference frame
        figure; imgray(CS.x.xn(:,:,B.zref));
        title([fn,', L=',int2str(L),', resolution =',...
            int2str(N),', blocksize =', int2str(blocksize),...
            ', reference frame'],'FontSize',12);
        xlabel(['reconstructed ',fn],'FontSize',12);

        
        %% CASE 8: -ulas- frame composition: partial 3d, 
        %{ 
        total variation and motion vector:
        tests an apriori known motion vector (shifted cameraman)
        A new class: Phi = TV_IC_MRA3_U
        is used with properties: motionV (motion vector) and use_motion (logical) 
        %}
% ------------------------------------------------------

        % create a shifted/rotated video

        fn='cameraman.bmp';
        %fn='cam512.png';
        L=10;
        unit=1;
        videoflag=2;

        signal2=Signal2D.make_fromImage(fn);
        %signal2.resize(128);
        signal= Signal_Factory.make_CyclefromSignal2D(signal2,videoflag,L,unit);

        % -- write the video to a file '~/Desktop/ffmpeg_example/'
        % -- change the directory in 'playData.m'
        % y=normalize_frame(signal3.xn); playData(y);
        

        params=struct;
        % best results for Phi= TV_IC_MRA3 and variable compression rates:
        Phi= TV_IC_MRA3_U(signal);   % improves when c varies over frames !!!
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
        % params.frameNumber=L;  % define for each frame new mask
        % OBSERVATION: quality gets worse for varying mask      
        [nodes,c,cmean]=B.nodesDefault(c,params);        
        
        % check, if nodes are ok
        assert(xor(iscell(nodes.data),nodes.fL==L),'#frames is set in nodes');
          
        % create compressive-sensing object:
        CS=CompressiveSensing(); % choosing _U version gives error
                                 % due to use of SampleTrafo_U
        CS.set_original(signal);
        CS.set_methodActive(4);  % DR:1, FISTA: 2
        disp(CS.get_methodActive);
        % setting epsilon >0 ---> improves rescon quality for single image but NOT for multiple frames  
        % CS.params.epsilon= 0.01*sqrt(signal.numel);
        
        % simulate measurements with compression c and reconstruct:
        % parameter setting       
        CS.params.epsilon=0.1;
        % PD
        if L>=9
            % reduce parameters critical for convergence
            CS.params.PD.sigma=0.05;
            CS.params.PD.tau=0.05;
            CS.params.PD.maxIters=500;
        end
        % DR
        CS.params.DR.gamma=0.1;
        CS.params.DR.maxIters=200;
        
        %decide to use motion
        Phi.motionV=unit;
        Phi.use_motion=true;
        
        CS.sim_SampledONS(Phi,B,nodes);
        
        y=normalize_frame(CS.x.xn);
        playData(y);


        %%
        CS.x.graph_totalVar;
        CS.x.graph_motionTV;
        % decompose reconstructed signal x:
        CS.x2C;
        CS.Phi.show_trafo;
        
        
        
%% --- do some additional analysis:
        
        prepfigure(1,[],CS.figopt);
        subplot(2,2,1);
        % check cross similarity between frame pairs
        % 1.) does the reconstruction incorporate temporal changes? 
        CS.x.graph_PSNRcross([],false);
        % 2.) what temporal changes occur in original video stream?
        subplot(2,2,2);
        CS.xorig.graph_PSNRcross([],false);
        % 3.) do the reconstructed frames match best with their synchronous original
        % frame?
        subplot(2,2,3);
        CS.xorig.graph_PSNRcross(CS.x,false);
        
        prepfigure(1,[],CS.figopt);
        subplot(2,2,1);
        % check cross similarity between frame pairs
        % 1.) does the reconstruction incorporate temporal changes? 
        CS.x.graph_SSIMcross([],false);
        % 2.) what temporal changes occur in original video stream?
        subplot(2,2,2);
        CS.xorig.graph_SSIMcross([],false);
        % 3.) do the reconstructed frames match best with their synchronous original
        % frame?
        subplot(2,2,3);
        CS.xorig.graph_SSIMcross(CS.x,false);


        % analyze sparsification achieved by Phi:
        Phi.set_signal(signal);
        Phi.dec;
        if ~exist('CS','var')
            CS=CompressiveSensing_U();
            CS.set_original(signal);
            if exist('nodes','var')
                CS.set_APhiNodes(Phi,B,nodes);           
            end
        end
        CSrate=CS.s1w_Donoho_AsymptotiCS;
        Phi.graph_sparsityDefect(true,CSrate);   
        
        Phi.show_trafo;
        
        thresh=Phi.computeThreshold(CSrate);
        Phi.hardthresh1_inplace(thresh);   
        Phi.fig=max(2,Phi.fig);
        Phi.show_trafo;

        
        % check TV of reconstruction and compare to TV of original
        f= TV3_Temporal(signal);
        f.dec;
        f.graph_trafo;
        
        f= TV3_Temporal(CS.x);
        f.dec;
        f.fig=3; f.graph_trafo;
        
        % increasing snweight should decrease TV of video:
        Phi.anweight1=0.1;
        CS.sim_SampledONS(Phi,B,nodes);          
        
        f= TV3_Temporal(CS.x);
        f.dec;
        f.fig=3; f.graph_trafo;
        
