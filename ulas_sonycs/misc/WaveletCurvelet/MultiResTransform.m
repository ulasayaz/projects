classdef MultiResTransform < FrameTrafo
    % abstract n-dim. multi-resolution transformations with graphical output
    % implementation classes specify dimension.   
    %
    % GT Berlin, 2011-2013
    %
    
    properties (Access=public)
        
        frameSet %@<cellstr> additional basis names to be studied
        
        % computed values
        
    end
    
    properties (SetAccess=protected)
        wvfamily      % @(string) e.g. haar,db,coif
        wvnumber      % (int)
        dwtmode_active %@ (string) DWT extension mode handling border distortion problem
        wcL            %  coarsest level of multi-resolution decomposition (size is 2^wcL)
        
        detail_str2idx, detail_idx2str %@<container> tranlates indices to labels for detail signals
        
    end
    
    
    %% constructor and commands
    methods
        
        function obj=MultiResTransform(signal)
            % constructor ~(signal)            
            if nargin==0 
                signal=SignalClass();
            end            
            obj = obj@FrameTrafo(signal);
                         
            obj.tol=1e-6;                                                            
            obj.set_wcoarsestlev(2);  % keep 2 levels away from bottom most
            
            obj.energy_compressionlevel=0.999;
            obj.repfun=@(x) sign(x).*sqrt(abs(x));
            obj.use_repfun=true;           
            obj.frameSet={'db1','db2','db3'};
            
            obj.frame_issparse=true;  % e.g. wavelets basis is sparse expressed in standard ONB
            
            % do not yet check invariant for abstract class
        end % constructor
        
        
        function set_basisname(obj, wvn)
            % ~(wvn) set multi-resolution name and reset transform
            obj.requireOnly(obj.isvalid_framename(wvn),'local',' wvn is valid frame name');
            p=regexp(wvn,'\d+'); % finds first digit
            if ~isempty(p) && p(1)>1
                obj.wvfamily=wvn(1:p(1)-1);
                obj.wvnumber=str2num(wvn(p(1):end));
            else
                obj.wvfamily=wvn;
                obj.wvnumber=[];
            end
            obj.C=[];
            obj.ensureOnly(~obj.isTrafodone,'local', 'transform reset');
        end
        
        
        function set_algo(obj)
            obj.algo.version=0.9;
            obj.algo.versiondate='5.12.2013';
            obj.algo.name='MRA';
            obj.algo.toolbox=[];
            obj.algo.mfilename=class(obj);
            obj.algo.author='G. Troll';
        end
        
        function set_wcoarsestlev(obj, lev)
            %  ~(lev) set coarsest level of multi-resolution decomposition (size is 2^L)
            obj.requireOnly(min(obj.size)<8 || (lev>=2 && lev<=obj.wmaxlev-1),'local',...
                ' level restrictions due to image size');
            obj.wcL=lev;
        end
        
        
        function set_deepestlev(obj, lev)
            %  ~(lev) set deepest level indirectly via wcL (size is 2^L)
            obj.requireOnly(max(obj.size)==0 || (lev>=1 && lev<=obj.wmaxlev-2),'local',...
                ' level restrictions due to image size');
            wcL_new=min(obj.size_dyadic)-lev;
            if wcL_new~=obj.wcL
                obj.wcL=wcL_new;
                obj.reset_Trafo;
            end
        end
        
        
    end
    
    %% abstract methods
    
    methods (Abstract)
        
        L=level(obj)            % actual level of multi-resolution decomp
        coeff=detcoef(obj,level,type) % detail coefficients of multi-resolution transform
        A= appcoef(obj)  % approximation coefficients
        t=detcoef_rangemax(obj,level) %max value of parameter type of the detail coefficients    		
    end
    
    %% queries
    methods
        
        function yn=adjointSynthesis(obj,x)
            % yn=~(x) apply adjoint of synthesis operator Phi to vector x
            % assuming that the adjoint(Phi)=Psi (orthogonality);
            % result is a vector yn;
            yn=obj.analyze(x);
        end
        
        function w2=clone(obj)
            % w2=~() clone object (shallow)
            constructor = str2func(class(obj));
            w2=constructor();
            w2.set_signal(obj.ts);
            w2.fig=obj.fig;
            w2.figopt=obj.figopt;
            w2.frameSet =obj.frameSet;
            w2.repfun=obj.repfun;
            w2.use_repfun=obj.use_repfun;
            w2.energy_compressionlevel =obj.energy_compressionlevel;
            w2.tol=obj.tol;
            
            w2.wvfamily=obj.wvfamily;
            w2.wvnumber =obj.wvnumber;
            w2.dwtmode_active =obj.dwtmode_active;
            w2.wcL=obj.wcL;
            w2.algo=obj.algo;
            
            w2.C=obj.C;
            w2.frame=obj.frame;
            w2.sigma=obj.sigma;
        end
               
        
        function wvn= basisname(obj)
            % multi-resolution short name
            wvn=[obj.wvfamily,num2str(obj.wvnumber)];
        end
        
        function nv=frame_length(obj)
            % n=~() number of elements of transform
            % must be redefined for some subclasses
            if ~obj.isTrafodone
                obj.dec;
            end
            nv=obj.numelC; %numel(obj.C2vec);
        end
        
        function nv=frame_norm(obj)
            % l2-norm of transform coefficients
            % tests in this form only isometries (e.g. ONS)
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            nv=1;  % norm(obj.C,2)/obj.ts.norm(2) will not work in matlab toolbox for higher dbx
        end
        
        function nv=frame_dim(obj)
            % n=~() dimension of Hilbert space
            nv=obj.ts.numel;  % default is redundancy 1
        end
        
        function nv=detcoefNorms(obj, lvl)
            % nv=~() (rel.) norm of detail coefficients of each sector
            % (from type 1 to detcoef_rangemax) from level 1 to obj.level
            if nargin <2
                nv=cell(obj.level,1);
                totalnorm=norm(obj.C2vec);
                for j=1:obj.level
                    K=obj.detcoef_rangemax(j);
                    for k=1:K
                        nv{j}{k}=norm(reshape(obj.detcoef(j,k),[],1),2)/totalnorm;
                    end
                end
            else
                j=lvl;
                totalnorm=norm(obj.C2vec);
                K=  obj.detcoef_rangemax(j);
                nv=zeros(1,K);
                for k=1:K
                    nv(k)=norm(reshape(obj.detcoef(j,k),[],1),2)/totalnorm;
                end
            end
        end
        
        function nv=appcoefNorm(obj)
            % nv=~() (rel.) norm of detail coefficients of each sector
            % (from type 1 to detcoef_rangemax) from level 1 to obj.level           
            totalnorm=norm(obj.C2vec);
            nv=norm(reshape(obj.appcoef(),[],1),2)/totalnorm;
        end
        
        function nv=coeffNorms(obj)
            % nv=~(lvl) norms of coefficients at computed level
            lvl=obj.level;
            nvtrend=obj.appcoefNorm();
            nv=obj.detcoefNorms(lvl); 
            nv=vertcat(nvtrend,nv(:));
            nv=nv/sum(nv);
        end
        
        function LN=wmaxlev(obj)
            % multi-resolution decomposition level
            % allow for different extensions in each dimension
            % (take dimension with largest extension)
            LN=max(0,floor(log2(max(obj.size))));
            %LN=max(0,min(floor(log2(min(obj.size)))));
        end
        
        function L=wcoarsestlev(obj)
            % L= ~() coarsest level of multi-resolution decomposition (size is 2^L)
            L=obj.wcL;
        end
        
        function L=deepestlev(obj)
            % deepest level chosen by wcoarsestlev
            L=max(1,min(obj.wmaxlev,min(obj.size_dyadic)-obj.wcoarsestlev));
        end
        
        function L=level_bottomUp(obj)
            L=min(obj.size_dyadic)-obj.level();
        end
        
       function ss=size(obj,n)
            % s=~(): sample size
            if nargin ==2
                ss=size(obj.sdata,n);
            else
                ss=size(obj.sdata);
            end
        end
        
        function J=size_dyadic(obj)
            % J=~() log2 of size (take dimension with lowest extension)
            J=floor(log2(obj.size));
        end
               
    end
    
    methods (Hidden)               
        
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
            offset=1;
            if isempty(lvl)
                for j=1:obj.level
                    K=obj.detcoef_rangemax(j);
                    for k=1:K
                        Ccell{offset+k}=obj.detcoef(j,k);
                    end
                    offset=offset+K;
                end
            else
                j=lvl;
                K=obj.detcoef_rangemax(j);
                for k=1:K
                    Ccell{offset+k}=obj.detcoef(j,k);
                end
            end
        end
              
        
        
        function F2=RandnMatTimesDualFrame(obj,m,n)
            % F2=(m,n) multiply random matrix randn(m,n) with dual frame
            % used e.g. for compressive sensing problems
            % Wavelets, curvelets are a Parseval frame (S=Id) and we can avoid
            % computing the frame matrix theta' (synth. op. or analysis op.
            % theta) by using the op. dec:
            % A*S^(-1)*theta'= A*theta', hence
            % theta*A'= dec applied to columns of random matrix A'.
            % the final result is the adjoint (complex transpose).
            obj.require(rem(sqrt(n),1)==0,'local', 'n is square number');
            n1=sqrt(n);
            n2=n1;
            F2=Frame();
            for j=1:m
                obj.set_signal(randn(n1,n2));
                obj.dec;  % apply decomposition to do frame analysis theta
                if j==1
                    F2.set_data(zeros(m,length(obj.C2vec)));
                end
                F2.data(j,:)=obj.C2vec';  % fill rows to implement adjoint
            end
        end
        
        
    end
    
    methods (Hidden)
        
        
    end
    
    %% filters
    % cf. superclass FrameTrafo
    methods
        
        function estimateNoise(obj)
            % estimate noise by measuring the standard deviation of the first
            % diagonal fluctutation of an orthogonal(!) wavelet transform of NoisyImage.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            
            % determine first diagonal detail subsignal which should not
            % contain much of the original signal if the original signal does not fluctuate too wildly
            D1=reshape(obj.detcoef(1,'d'),[],1);
            % theory predicts that Gaussian noise is invariant under orthonormal
            % wavelet transforms, i.e. if the first condition is satisfied,
            % the standard deviation of the first diagonal detail subsignal
            % should be a good approximation to the standard deviation of the
            % added Gaussian noise.
            % obj.sigma=obj.frame_norm*std(D1(:));
            % we use the correction factor obj.frame_norm not in sigma but
            % later as thresh=4*obj.frame_norm*obj.sigma
            
            % Donoho, Johnston (1994)
            % wavelet coefficient of uncorrupted signal are sparse
            % and can be viewed as outliers in the wavelet coeff. of the
            % finest scale:
            obj.sigma=median(abs(D1-median(D1)))/0.6745;
            if obj.sigma<1e-9
                obj.sigma=mean(abs(D1-mean(D1)))/0.6745;
            end
        end
        
        
    end
    
    
    %% tests
    
    methods
                        
        
    end
    
    methods (Static)
        
        function ok= isvalid_framename(wvn)
            ok=ischar(wvn);
        end
                
        
    end
    
    %% graphics
    methods
                
        
        function graph_distribution(obj, open_new,lev)
            % ~() show distribution of frame coefficients in open window;
            % works for any dimension.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || isempty(open_new) || islogical(open_new),'local', 'boolean');
            
            if nargin <2 || isempty(open_new)
               open_new=true;
            end
            if open_new
                 prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            if nargin >=3 && ~isempty(lev)
                dwt=flatten2vec(obj.detcoef(lev));  % for dim>1 collects all filters of level lev               
                levelshown=lev;
            else
                dwt=obj.C2vec;
                levelshown=obj.level;
            end
                       
            HM=HMatrix(dwt);
            HM.graph_distribution(false);
            
            h=get(gca,'Title');
            present_titstr=get(h,'String');
            addtit=['distribution of coefficients ',...
                obj.algo.name,'(',obj.basisname,', level: ',num2str(levelshown),')'];
            new_titstr=present_titstr;
            new_titstr{1}=addtit;
            if ~iscell(present_titstr)
                new_titstr{2}=present_titstr;
            end
            title(new_titstr,'fontsize',12);
            
        end
        
        function graph_normDistrib(obj,open_new)
            % ~() plot norms of trend and detail signal at computed level of MRA.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            lvl=obj.level();
            nv=obj.coeffNorms();
            plot(nv);  
            filter_idx=[1;reshape(sort(cell2mat(obj.detail_str2idx.values)),[],1)+1];
            filter_labels=vertcat({'T'},reshape(obj.detail_idx2str.values,[],1)); 
            set(gca,'XTick',filter_idx);
            set(gca,'XTickLabel',filter_labels);
            ylabel('||:||_2');
            xlabel('coeff. type');
            tittext{1}=['Norms of coeff. ',obj.algo.name,'(',obj.basisname,...
                ', lev=',num2str(lvl),')'];
            tittext{2}=obj.ts.signalname;
            title(tittext,'fontsize',obj.fontsize);
            
        end
        
        function graph_detcoefNorms(obj, open_new)
            % ~() plot norms of detail signals up to computed level of MRA.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            nv=obj.detcoefNorms();            
            L=length(nv);
            legstr=cell(1,L);
            for j=1:L
               semilogy(cell2mat(nv{j}));  
               if j==1
                   hold all;
               end
               legstr{j}=['level ', num2str(j)];
            end
            if ~isempty(obj.detail_str2idx)
                filter_idx=sort(cell2mat(obj.detail_str2idx.values));
                filter_labels=obj.detail_idx2str.values;            
                set(gca,'XTick',filter_idx);
                set(gca,'XTickLabel',filter_labels);
            end
            legend(legstr,'location','best');
            ylabel('||:||_2');
            xlabel('detail signal (filter)');
            tittext{1}=['Norms of detail signals of ',obj.algo.name,'(',obj.basisname,')'];
            tittext{2}=obj.ts.signalname;
            title(tittext,'fontsize',obj.fontsize);
            
        end
                             
               
        function c=graph_energymap(obj, open_new)
            % c=~() energy map of transform in open window returns compression
            % observe: orthogonal multi-resolution transforms only conserve l2-norm
            % l2-norm is essential for quality of reconstruction
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            enthresh=obj.energy_compressionlevel;
            enthresh_str=[num2str(100*enthresh,'%3.2f'),'% point'];
            
            sC=obj.C2vec();  % here (:), but more complicated e.g. for curvlets
            sC=sort(sC.^2,'descend');
            ennorm=sum(sC);
            enmap=cumsum(sC)/ennorm;
            L=length(enmap);
            idx99=find(enmap>=enthresh,1,'first');
            c=prod(obj.N)/idx99;
            
            % just plot curve without idx99, to make it possible to assemble several
            % energy maps for different frames
            semilogx((1:L)/L,enmap); %,idx99/L,enmap(idx99),'or');
            % text(idx99/L,enmap(idx99),num2str(idx99/L,'%3.2g'),'VerticalAlignment','bottom');
            
            xlabel('share of size-ordered coefficients kept');
            ylabel('normed ||.||_2');
            title({[obj.basisname,'-Energy Map of ',...
                obj.ts.get_signalname], ...
                ['redundancy=',num2str(obj.frame_redundancy,'%3.1f'),...
                ', compression=', num2tex(c,'%3.1f','none'),...
                ' at ',enthresh_str]},...
                'fontsize',obj.fontsize);
            ylims=ylim;
            ylim([ylims(1),1.005]);
        end
        
        function graph_sortedCoeff(obj,func,open_new)
            % ~(meth) compare several transforms w.r.t. func in open window
            % e.g. func.f=@(x) x; func.str= '|.|'; func.descr='modulus';
            %      func.f= @(x) cumsum(x.^2); func.desc='||.||_2'; func.descr='energy map';
            %
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.require(isstruct(func) && isa(func.f,'function_handle'),'local',...
                'func.f is function, func.descr string');
            if nargin <3 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            if ~isfield(func,'descr')
                func.descr=func2str(func.f);
            end
            if ~isfield(func,'str')
                func.str=func2str(func.f);
            end
            
            Cs=obj.C2vec;
            Cs=sort(abs(Cs(:)),'descend');
            L=length(Cs);
            x=(1:L)/L;
            y=func.f(Cs);
            y=y/max(y);
            %semilogx(x,y);
            loglog(x,y);
            xlabel('share of size-ordered coefficients kept');
            ylabel(func.str);
            title([func.descr,': MRA ', obj.basisname,', Signal ',...
                obj.ts.signalname],'fontsize',obj.fontsize);
        end
                                      
        
        function show_trafos_selection(obj, frameSet)
            % ~(frameSet) show and compare several multi-resolution transforms
            % in a new window;
            obj.require(nargin <2 || iscellstr(frameSet),'local','frameSet is cell of strings');
            obj.requireOnly(~isempty(obj.sdata),'local','signal image is non-empty');
            obj.requireOnly(obj.fig~=0,'local', 'figure output activated');
            if nargin <2 || isempty(frameSet)
                frameSet=obj.frameSet;
            end
            prepfigure(obj.fig,obj.algo,obj.figopt);
            suptitle(obj.algo.toolbox);
            Lw=length(frameSet);
            sd = factor_subplots( Lw+1 );
            subplot(sd(1),sd(2),1);
            obj.graph_signal(false);
            
            for j=1:Lw
                subplot(sd(1),sd(2),j+1);
                obj.set_basisname(frameSet{j});
                obj.dec;
                obj.graph_trafo(false);
            end
            % keep result of last basisname for further graphical output
            
            
        end
        
        
        function show_resultd(obj)
            % show decpmposition result in a new window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.require(obj.fig~=0,'local', 'figure output activated');
            
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            suptitle(['Toolbox: ',obj.algo.toolbox],obj.fontsize+1);
            
            subplot(2,2,1);
            obj.graph_signal(false);
            subplot(2,2,2);
            obj.graph_trafo(false);
            subplot(2,2,3);
            obj.graph_energymap(false);
            subplot(2,2,4);
            obj.graph_detcoefNorms(false);
%             imagesc(obj.appcoef());
%             title(['coarsest approx. coeff. level ',num2str(obj.level),' (linear scale!)'],...
%                 'fontsize',obj.fontsize);
%             colorbar;
            
        end
        
        
        function show_noisefilter(obj, thresh)
            % ~(thresh) apply noise filter with threshold and show in new win
            % apply e.g. to a test signal created by w= make_SineSignal(5,0.1).
            obj.require(obj.fig~=0,'local', 'figure output activated');
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            
            if nargin < 2
                thresh=[];
            end
            
            [s2, thresh]=obj.hardthresh1(thresh);
            method_str=[obj.basisname,', thresh=',num2str(thresh,'%3.2e')];
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            subplot(1,2,1);
            obj.graph_signal(false);
            
            subplot(1,2,2);
            s2.graph_signal(false);
            title({'MRA-Rauschfilter: ',method_str}, 'fontsize',obj.fontsize);
        end
        
        
        function tit2= title2(obj,levels,disp_dwtmodes)
            % second title line
            if nargin <3
                disp_dwtmodes=true;
            end
            if nargin <2 || isempty(levels)
                levels=obj.level;
            end
            if length(levels)>1
                level_str=[num2str(min(levels)),'...',num2str(max(levels))];
            else
                level_str=num2str(levels);
            end
            tit2=['Level ',level_str];
            if disp_dwtmodes && ~isempty(obj.dwtmode_active)
                tit2=[tit2,', DWTmode=',obj.dwtmode_active];
            end
        end
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal class';
            ok= isa(obj.ts,'SignalClass');
        end
    end
    
end





