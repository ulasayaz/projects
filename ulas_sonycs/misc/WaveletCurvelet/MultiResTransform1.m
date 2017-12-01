classdef MultiResTransform1 < MultiResTransform
    % abstract 1-dim. multi-resolution transformations with graphical output
    % implementation classes, e.g.
    % - Wavelet1
    % - Wavelet1_wlab
    %
    % GT Berlin, 2011-2013
    %
    
    properties (Access=public)       
        %graphics
        marker %@<string> plot symbol for single points
    end
        
    %% constructor and commands
    methods
        
        function obj=MultiResTransform1(signal)
            % constructor ~(signal)            
            if nargin==0 
                signal=TimeSignal();
            end
            assert(isa(signal,'TimeSignal'),'1d signal');
			obj = obj@MultiResTransform(signal);
			            
            obj.ensure(true,'local','check invariant');
        end % constructor
        
        
        function set_deepestlev(obj, lev)
            %  ~(lev) set deepest level indirectly via wcL (length is 2^L)
            obj.requireOnly(max(obj.length)==0 || (lev>=1 && lev<=obj.wmaxlev-2),'local',...
                ' level restrictions due to signal length');
            wcL_new=min(obj.length_dyadic)-lev;
            if wcL_new~=obj.wcL
                obj.wcL=wcL_new;
                obj.reset_Trafo;
            end
        end
        
        
    end
        
    
    %% queries
    methods
        
        function w2=clone(obj)
            % w2=~() clone object (shallow)
			clone@MultiResTransform(obj);
            w2.marker=obj.marker;
        end
                       
        
        function s=xn(obj)
            % s=~() signal data
            s=obj.ts.xn;
        end
        
       
        
        function LN=wmaxlev(obj)
            % multi-resolution decomposition level
            % allow for different extensions in each dimension
            % (take dimension with largest extension)
            LN=max(0,floor(log2(max(length(obj.xn)))));
            %LN=max(0,min(floor(log2(min(length(obj.xn))))));
        end
               
        
        function L=deepestlev(obj)
            % deepest level chosen by wcoarsestlev
            L=min(obj.wmaxlev,min(obj.length_dyadic)-obj.wcoarsestlev);
        end
        
        function L=level_bottomUp(obj)
            L=min(obj.length_dyadic)-obj.level();
        end
        
        function ss=length(obj)
            % s=~(): sample length with padding inmcluded
            ss=length(obj.xn);
            
        end
        
        function J=length_dyadic(obj)
            % J=~() log2 of length (take dimension with lowest extension)
            J=floor(log2(obj.length));
        end
        
        function J=size_dyadic(obj)
            % J=~() log2 of length (take dimension with lowest extension)%
            % replaces superclass definition by checking only 1 dim
            J=floor(log2(obj.length));
        end             
        
        
    end
    
    methods (Hidden)
        
        function t=detcoef_rangemax(obj,level)
            % t=~() max value of parameter type of the detail coefficients
            t=1;
        end
        
    end
    
    
    %% transforms
    methods
	
         function mat=C2graph(obj)
            % mat=~() convert coefficient vector to graph format            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            mat=obj.C(:);
        end
        
        
    end
    
    
    %% filters
    % cf. superclass FrameTrafo
    methods
        
        
    end
    
    
    %% tests
    
    methods
        
        function show_test(obj,omega,hN)
            % ~([sig,hN]) testing transform with star signal of v vertices
            % to examin directional characteristic.
            if ~exist('omega','var') || isempty(v)
                omega=5;
            end
            if ~exist('hN','var') || isempty(hN)
                hN=256;
            end
            hts=TimeSignal.make_sinus(omega,0,hN);
            obj.set_signal(hts);
            obj.set_deepestlev(2);
            obj.figopt.pixsizeX=1000;
            
            obj.dec;
            obj.show_resultd;
        end
        
        
    end
    
    methods (Static)
        
        function ok= isvalid_framename(wvn)
            ok=ischar(wvn);
        end
        
        function best= bestwavelet(signal)
            % select daubechy wavelet with best energy compression at threshold
            assert(nargin>0 && isa(signal,'TimeSignal'),'class of signal is TimeSignal');
            
            if license('test','Wavelet_Toolbox')
                w=Wavelet1D();
            else
                w=Wavelet1D_wlab();
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
            
            if nargin <3 || isempty(lev)
                lev=obj.level;
            end
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            levelshown=lev;
            
            if obj.use_repfun
                xn_dwt= obj.repfun(obj.C);
                xn_dwt(xn_dwt<-10)=-10;
            else
                xn_dwt=obj.C;
            end
            
            plot(xn_dwt);
            
            if obj.use_repfun
                legend(func2str(obj.repfun),'location','best');
            end
            
            tittext=[obj.algo.name,'(',obj.basisname,'), ',num2str(levelshown),' levels shown'];
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
            plot(A);
            tittext=['Trendsignal ',num2str(obj.level),' ',obj.algo.name,' ',obj.basisname];
            tit2=obj.title2;
            if isempty(tit2)
                title(tittext,'fontsize',obj.fontsize);
            else
                title([tittext,' (',tit2,')'],'fontsize',obj.fontsize);
            end
            
        end
        
        function graph_detail(obj,n,open_new)
            % ~([n])  show fluctuation at level n in open figure window
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin <3 || n<=obj.level,'local',...
                'n is smaller than the highest calculated level');
            if nargin <3
                n=length(obj.L)-2;  % last detail signal
            end
            if nargin <4 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            D=obj.detcoef(n);
            if obj.use_repfun
                D= obj.repfun(D);
                D(D<-10)=-10;
            end
            
            plot(D);
            
            if obj.use_repfun
                cblabel(func2str(obj.repfun));
            end
            
            tittext=['Detailsignal ',obj.algo.name,' ', obj.basisname];
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
            
            if nargin <3 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
                subplot(1,2,1);
                plot(noisy.xn);
                
                title(['noisy (\sigma=',num2str(signoise),', PSNR=',num2str(PSNR_noisy,'%3.1f'),' dB, SSIM=',...
                    num2str(SSIM_noisy,'%3.2f'),')'],'fontsize',obj.fontsize);
                subplot(1,2,2);
            end
            plot(denoised.xn);
            
            title(['denoised ',obj.algo.name,' (PSNR=',num2str(PSNR,'%3.1f'),' dB, SSIM=',...
                num2str(SSIM,'%3.2f')],'fontsize',obj.fontsize);
            
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
                plot(cC{j});
                
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
            ok= isa(obj.ts,'TimeSignal');
        end
    end
    
end




