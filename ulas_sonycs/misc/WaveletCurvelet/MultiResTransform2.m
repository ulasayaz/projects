classdef MultiResTransform2 < MultiResTransform
    % abstract 2-dim. multi-resolution transformations with graphical output
    % implementation classes, e.g.
    % - Wavelet2_mlab
    % - Wavelet2_wlab
    % (analogous to Fourier2)
    %
    % GT Berlin, 2011-2013
    %
    
    properties (Access=public)
       
        
    end
        
    
    %% constructor and commands
    methods
        
        function obj=MultiResTransform2(signal)
            % constructor ~(signal)                     
            if nargin==0 
                signal=Signal2D();
            end
            assert(isa(signal,'Signal2D'),'2d signal');
			obj = obj@MultiResTransform(signal);
			            
            obj.ensure(true,'local','check invariant');        
        end % constructor
       
    end
    
   
    %% queries
    methods               
        
        function v=img(obj)
            % signal
            if ~isempty(obj.ts)
                v=obj.ts.xn;
            else
                v=[];
            end
        end
        
        function L=deepestlevAnisotropic(obj)
            % deepest level in the deepest direction
            L=obj.wmaxlev-obj.wcoarsestlev;
        end
                     
        
    end
    
    methods (Hidden)
        
        function t=detcoef_rangemax(obj,level)
            % t=~() max value of parameter type of the detail coefficients
            t=3;
        end
        
    end
    
    
    %% transforms
    methods
        
        function Ccell=C2cell(obj,lvl)
            % convert transformation result to cell of components
            % must be redefined for those subclasses, where C is not a cell
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || isempty(lvl) || (lvl>=1 && lvl<=obj.level),'local',...
                'level is below computed level');
            if nargin <2
                lvl=[];
            end
            Ccell{1}=obj.appcoef();
            if isempty(lvl)
                for j=1:obj.level
                    Ccell{3*(j-1)+2}=obj.detcoef(j,'h');
                    Ccell{3*(j-1)+3}=obj.detcoef(j,'v');
                    Ccell{3*(j-1)+4}=obj.detcoef(j,'d');
                end
            else
                j=lvl;
                Ccell{2}=obj.detcoef(j,'h');
                Ccell{3}=obj.detcoef(j,'v');
                Ccell{4}=obj.detcoef(j,'d');
            end
        end
        
        function mat=C2graph(obj)
            % mat=~() convert coefficient vector to image
            % will be redefined in child classes
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            mat=reshape(obj.C,[obj.frame_length,obj.frame_length]);
        end
        
        
       
    end
    
    methods (Hidden)
        
        
    end
    
    %% filters
    % cf. superclass FrameTrafo
    methods
        
       
    end
    
    
    %% tests
    
    methods
        
        function w=test1_distortion(obj,omega,sig)
            % test noisefilter deformation of a simple harmonics
            assert(nargin <1 || isnumeric(omega),'omega is angular freq');
            
            % testing the noisefilterfilter
            if nargin ==0 || isempty(omega)
                omega=5;
            end
            if nargin <2 || isempty(sig)
                sig=0;  % noise
            end
            hN=128; del0=[0,0.1,pi/4];
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            suptitle(['Rauschfilterverzerrung einer einfachen Harmonischen - ',...
                'Method : noisefilter\_simple'], obj.fontsize+1);
            
            for j=1:3
                subplot(3,3,j);
                %  choose del0 ~=0 to construct a signal with x[N]~=x[0]
                %  del0=0; del0=0.1; del0=pi/4;
                
                hts=Signal2D.make_sinus(omega,sig,hN,del0(j));
                w=obj.make();
                w.set_signal(hts);
                
                w.dec;
                w.estimateNoise();
                thresh=4.5*w.sigma;
                
                yn=w.hardthresh1(thresh);
                err_max=max(reshape(abs(yn-w.img),[],1));
                err_med=medianNaN(reshape(abs(yn-w.img),[],1));
                
                imagesc(w.img); colormap(w.ts.colormap_active); colorbar;
                title({['Signal (deviation from full period \delta=',...
                    num2str(del0(j),'%3.1f'),')']...
                    },...
                    'fontsize',obj.fontsize);
                
                subplot(3,3,j+3);
                imagescn(yn-w.img); colormap(w.ts.colormap_active);colorbar;
                title({['error[max, median])= ',...
                    vec2str([err_max, err_med],'%3.1e')]...
                    }, 'fontsize',obj.fontsize);
                
                subplot(3,3,j+6);
                imagesc(w.C); colormap(w.ts.colormap_active); colorbar;
                title({['Wavelet Trafo,' w.basisname,],['thresh=',...
                    num2str(thresh)]},...
                    'fontsize',obj.fontsize);
            end
            
        end
        
        function show_test(obj,v,hN)
            % ~(w,[v,hN]) testing transform with star signal of v vertices
            % to examin directional characteristic.
            if ~exist('v','var') || isempty(v)
                v=8;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            hts=Signal2D.make_star(v,hN);
            obj.set_signal(hts);
            obj.set_deepestlev(2);
            obj.figopt.pixsizeX=1000;
            
            obj.dec;
            obj.show_resultd;
        end
        
        
    end
    
    methods (Static)
        
        function best= bestwavelet(signal)
            % select daubechy wavelet with best energy compression at threshold
            assert(isa(signal,'Signal2D'),'class of signal is Signal2D');
            
            if license('test','Wavelet_Toolbox')
                w=Wavelet2D_mlab();
            else
                w=Wavelet2D_wlab();
            end
            
            w.set_signal(signal);
            
            enthresh=w.energy_compressionlevel;
            nmax=20;
            wvn=zeros(1,nmax);
            j=1;
            while (j<=5 || wvn(j-1)<=max(wvn(j-3:j-2))) && j<=nmax
                try   % number j might not be defined in the toolbox used
                    w.set_basisname(['db',num2str(j)]);
                catch
                    break;
                end
                
                w.dec;
                sC=sort(w.C.^2,'descend');
                ennorm=sum(sC);
                enmap=cumsum(sC)/ennorm;
                idx99=find(enmap>=enthresh,1,'first');
                wvn(j)=idx99/numel(sC);
                j=j+1;
            end
            
            nmax=j-1;
            wvn=wvn(1:nmax);
            [~,bestn]=min(wvn(1:nmax));
            if wvn(3)<=wvn(bestn)+1  % take db3 if idx99 differs only by 1
                [~,idx]=sort(wvn,'ascend');
                bestn=idx(2);
            end
            
            best=['db',num2str(bestn)];
            
        end
        
    end
    
    %% graphics
    methods
        
        
        function graph_trafo(obj, open_new, lev)
            % show transform C in a single matrix in open figure window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            obj.requireOnly(nargin<3 || (isnumeric(lev) && lev>0),'local', 'integer');
            
            if ~exist('lev', 'var') || isempty(lev)
                lev=1;
            end
            if ~exist('open_new', 'var') || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            [img_dwt, levelshown]=obj.C2graph(lev);
            if obj.use_repfun
                img_dwt= obj.repfun(img_dwt);
                img_dwt(img_dwt<-10)=-10;
            end
            
            imagesc(img_dwt);
            colormap('default');
            colorbar;
            if obj.use_repfun
                cblabel(func2str(obj.repfun));
            end
            
            tittext=[obj.algo.name,'(',obj.basisname,'), ',num2str(levelshown),...
                ' levels shown',obj.add_signalname];
            tit2=obj.title2;
            if isempty(tit2)
                title(tittext,'fontsize',obj.fontsize);
            else
                title({tittext,[' (',tit2,')']},'fontsize',obj.fontsize);
            end
            
        end
              
        
        function graph_trend(obj,open_new)
            % ~([n]) show trend signal n in open figure window
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            A= obj.appcoef();
            imagesc(A); colormap(obj.ts.colormap_active); colorbar;
            tittext=['Trendsignal ',num2str(obj.level),' ',obj.algo.name,' ',obj.basisname];
            tit2=obj.title2;
            if isempty(tit2)
                title(tittext,'fontsize',obj.fontsize);
            else
                title([tittext,' (',tit2,')'],'fontsize',obj.fontsize);
            end
            
        end
        
        function graph_detail(obj,n,orient,open_new)
            % ~([n,orient])  show fluctuation at level n in open figure window
            % orient for proper wavelets is 'd','h', or 'v'.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin <3 || n<=obj.level,'local',...
                'n is smaller than the highest calculated level');
            if ~exist('n','var')
                % n=length(obj.S)-2;  % last detail signal
                n=1;
            end
            if ~exist('orient','var')
               orient=1;
            end
            if ~exist('open_new', 'var') || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            D=obj.detcoef(n,orient);
            if obj.use_repfun
                D= obj.repfun(D);
                D(D<-10)=-10;
            end
            
            imagesc(D);
            colormap(obj.ts.colormap_active);
            colorbar;
            if obj.use_repfun
                cblabel(func2str(obj.repfun));
            end
            
            if ~ischar(orient)
                orient=obj.detail_idx2str(orient);
            end
            
            tittext=['Detailsignal ',orient,num2str(n),' ',obj.algo.name,' ', obj.basisname];
            tit2=obj.title2;
            if isempty(tit2)
                title(tittext,'fontsize',obj.fontsize);
            else
                title([tittext,' (',tit2,')'],'fontsize',obj.fontsize);
            end
            
        end
        
        
        function [noisy,PSNR_noisy,SSIM_noisy]=graph_test_denoise(obj, signoise, open_new)
            % ~(signoise) add noise, estimate noise, denoise
            obj.requireOnly(nargin<3 || islogical(open_new),'local',...
                'open_new is boolean');
            if isscalar(signoise)
                noisy=Signal2D(obj.img+signoise*randn(size(obj.img)));
            elseif isa(signoise,'Signal2D')
                noisy=signoise;
            end
            
            PSNR_noisy=obj.ts.PSNR(noisy);
            SSIM_noisy=obj.ts.SSIM(noisy);
            
            constructor = str2func(class(obj));
            w=constructor(noisy);           
            w.dec;
            denoised=w.hardthresh1();
            
            PSNR=obj.ts.PSNR(denoised);
            SSIM=obj.ts.SSIM(denoised);
            
            if ~exist('open_new', 'var') || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
                subplot(1,2,1);
                imagesc(noisy.xn);
                axis image;
                colorbar;
                title(['noisy (\sigma=',num2str(signoise),', PSNR=',num2str(PSNR_noisy,'%3.1f'),' dB, SSIM=',...
                    num2str(SSIM_noisy,'%3.2f'),')'],'fontsize',obj.fontsize);
                subplot(1,2,2);
            end
            imagesc(denoised.xn);
            axis image;
            colorbar;
            title(['denoised ',obj.algo.name,' (PSNR=',num2str(PSNR,'%3.1f'),' dB, SSIM=',...
                num2str(SSIM,'%3.2f')],'fontsize',obj.fontsize);
            colormap('gray');
        end
        
        function show_trafo_components(obj,lvl)
            % show components of transform C in subplots of a new window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || (lvl>=1 && lvl<=obj.level),'local',...
                'level is below computed level');
            if nargin <2
                lvl=[];
            end
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            
            cC=obj.C2cell(lvl);
            L=length(cC);
            sd=factor_subplots(L);
            
            tit=[obj.algo.name,'(',obj.basisname,') of ',...
                obj.ts.get_signalname];
            if ~isempty(lvl)
                tit=[tit,', only level ',num2str(lvl)];
            end
            
            suptitle(tit,obj.fontsize);
            for j=1:L
                subplot(sd(1),sd(2),j);
                imagesc(cC{j});
                colorbar;
                if j==1
                    title(['approx. signal (level ',num2str(obj.level),')'],...
                        'fontsize',obj.fontsize);
                else
                    title(['detail signal j=',num2str(j)],'fontsize',obj.fontsize);
                end
            end
            
        end
         
               
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal class';
            ok= isa(obj.ts,'Signal2D');
        end
    end
    
end




