function [ xprelim,cs0,PSNR_prelim,SSIM_prelim ] = CS_SimPremimRecon( signal )
% STEP 1: simulate measurement and preliminary reconstruction 
% of each frame separately.
% -------------------------------------------------------------------------     
        display('running preliminary frame by frame reconstruction ...');
        xprelim=signal.make_like();
        xprelim.signalname=['preliminary recon of ',signal.signalname];
        xprelim.xn=zeros(signal.size);
        L=signal.size(3);
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
            disp(['j=',num2str(j),', [PSNR,SSIM]=',vec2str([PSNR_prelim(j),SSIM_prelim(j)])]);
            xprelim.xn(:,:,j)=cs0.x.xn;
            multiWaitbar('preliminary recon',j/L);
            
        end
        disp(' ');
        multiWaitbar('preliminary recon','close');             

end

