classdef MultiResTransform3_Partial<MultiResTransform3
    % 3d partial multiresolution transform
    % applying a 2d MRA to each z-plane separately;
    % used e.g. as one of the two transforms in FrameTrafoPair.
    %
    % Example
    % ========
    %     w=MultiResTransform3_Partial();
    %     signal2=Signal2D.make_fromImage('cameraman.bmp');
    %     L=3; signal3=Signal3D.make_CyclefromLowerDimSig(signal2,1,L);
    %     signal3.play_signal;
    %     w.set_signal(signal3);
    %     w.set_transform(Wavelet2D_mlab());
    %     w.dec;
    %     w.graph_trafo;
    % --- test inversion properties:
    %     w.test_DecRec(); w.test_AnalyzeSynthsize();
    %
    
    properties
        zref   %@<integer> reference z-plane for total variation recon
        anweight1 %<real> weight for transform in z-direction, e.g. for total variation
    end
    
    properties (SetAccess=protected)
        ftrafo    %@<MultiResTransform2>  basic 2d transform
    end
    
    
    %% commands
    methods
        
        function obj=MultiResTransform3_Partial(signal,hftrafo)
            % constructor
            if nargin==0
                signal=Signal3D();
            end
            obj = obj@MultiResTransform3(signal);
            if nargin<2
                if license('test','Wavelet_Toolbox')
                    hftrafo=Wavelet2D_mlab();
                else
                    hftrafo=Wavelet2D_wlab();
                end
            end            
            obj.ftrafo=hftrafo;
            obj.set_signal(signal);  % in order to take further actions
            obj.anweight1=0.5;
            obj.zref=max(1,ceil(signal.size(3)/2));  % mid point z-plane
            
            obj.ensure(true,'local','check invariant');
        end
        
        function set_algo(obj)
            set_algo@MultiResTransform3(obj);
            obj.algo.name='partial MRA3';
        end
        
        function set_signal(obj,hts)
            % ~(hts) set time signal to hts of class SignalClass
            obj.requireOnly(isa(hts,'Signal3D'),'local','ts is a 3d signal');
            obj.requireOnly(~isempty(obj.ftrafo),'local','needs ftrafo');
            if hts.size(3)~=obj.ts.size(3)
                obj.zref=max(1,ceil(hts.size(3)/2)); % mid point z-plane
            end
            set_signal@MultiResTransform3(obj,hts);
            if ~isempty(hts.xn)
                obj.ftrafo.set_signal(Signal2D(hts.xn(:,:,1)));
            end
            
        end
        
        function set_transform(obj, hftrafo1)
            % ~(hftrafo) set basic 2d transform
            obj.requireOnly(isempty(hftrafo1) || isa(hftrafo1,'MultiResTransform2'),...
                'local','needs 2d muliresolutin transform');
            obj.ftrafo=hftrafo1;
            obj.reset_Trafo();
            
            obj.colormap_active='default';
        end
        
        function set_dwtmode(obj, modestr)
            % ~(modestr) border extension mode cf. help dwtmode
            % e.g. 'sym' (symmetric extension), 'zpd' (zero-padding),
            % 'ppd' periodic, 'sp1 ' smooth padding
            obj.dwtmode_active=modestr;
            dwtmode(modestr,'nodisp');
        end
        
        function dec(obj,lev)
            % partial discrete multiresolution analysis
            % apply 2d transform to each z-projection
            obj.requireOnly(obj.ts.numel>0,'local','non-empty sample size');
            
            obj.ftrafo.ts.set_signal(obj.ts.xn(:,:,1));
            if nargin <2
                lev=obj.ftrafo.deepestlev();
            end
            
            L=obj.ts.size(3);
            obj.ftrafo.ts.set_signal(obj.ts.xn(:,:,1));
            obj.ftrafo.dec(lev);
            
            % alternative to object array would be cell array
            %             constructor = str2func(class(obj.ftrafo.C));
            %             obj.C(L)=constructor;
            obj.C=cell(L,1);
            obj.C{1}=obj.ftrafo.C;
            
            for j=2:L
                obj.ftrafo.ts.set_signal(obj.ts.xn(:,:,j));
                obj.ftrafo.dec(lev);
                obj.C{j}=obj.ftrafo.C;
            end
            
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wc)
            % yn=~(cC) % signal reconstruction of decomposition wc,
            % which must be of the same datatype as result of obj.dec.
            
            yn=zeros(obj.ts.size);
            for j=1:obj.ts.size(3)
                yn(:,:,j) = obj.ftrafo.rec(wc{j});
            end
        end
        
    end
    
    
    %% queries
    
    methods
                
        function w2=clone(obj)
            % w2=~() clone object (shallow)
            w2=clone@MultiResTransform3(obj);
            w2.ftrafo=obj.ftrafo.clone;
            w2.zref=obj.zref;
            w2.anweight1=obj.anweight1;
        end
        
        function wvn= basisname(obj)
            % multi-resolution short name
            wvn=obj.ftrafo.basisname;
        end
        
        function ok=isTrafodone(obj)
            % ok=~(): is result of DWT available?
            ok=~isempty(obj.C) && ~isempty(obj.C{1});
        end
        
        function nv=frame_norm(obj)
            % l2-norm of transform coefficients
            % tests in this form only isometries (e.g. ONS)
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            nv=obj.ftrafo.frame_norm;
        end
        
        function L=level(obj)
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            L= obj.ftrafo.level;
        end
        
        function nc=numelC(obj,cC)
            % nc=~([cC]) number of elements in decomposition
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            L=obj.ts.size(3);
            if nargin<2
                cC1=obj.C{1};
            else
                cC1=cC{1};
            end
            nc=L*obj.ftrafo.numelC(cC1);
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose vector x into vector yn
            % apply 2d transform to each z-projection
            L=obj.ts.size(3);
            x=reshape(x,obj.N);
            
            L=obj.ts.size(3);
            cC=reshape(obj.ftrafo.analyze(x(:,:,1)),[],1);
            yn=zeros(length(cC),L);
            yn(:,1)=cC;
            for j=2:L
                yn(:,j)=reshape(obj.ftrafo.analyze(x(:,:,j)),[],1);
            end
        end
        
        function xn= synthesize(obj,y)
            % xn=~(y) synthesize signal from decomposition y, where
            % y must beong to the same data type as result of analyze
            % apply 2d transform to each z-projection
            xn=zeros(obj.ts.size);
            L=obj.ts.size(3);
            y=reshape(y,[],L);
            for j=1:obj.ts.size(3)
                xn(:,:,j) = obj.ftrafo.synthesize(y(:,j));
            end
        end
        
        function n=nnz(obj)
            % n=~() number of non-zero elements in the transformation
            % coefficients
            L=obj.ts.size(3);
            n=zeros(1,L);
            for j=1:L
                n(j)=nnz(obj.ftrafo.C2vec(obj.C{j}));
            end
        end
        
        function [sd,p,rs]=sparsityDefect(obj,rs,cC)
            % d=~(rs) rel. sparsity defect of transformation coefficients
            % for relative sparsity rs in [0,1].
            % a.k.a. lp/error of best s/term approximation 
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || isempty(rs) || (rs>=0 && rs<=1),...
                'local','rs is rel. sparsity'); 
            if nargin<2 || isempty(rs)
                rs=0.1;
            end
            if nargin<3
                cC=obj.C;
            end            
            L=obj.ts.size(3);
            sd=zeros(L,1);
            for j=1:L
                [sd(j),p]=obj.ftrafo.sparsityDefect(rs,cC{j});
            end             
        end
        
        
    end
    
    methods (Hidden)
        
        function t=detcoef_rangemax(obj,level)
            % t=~() max value of parameter type of the detail coefficients
            if nargin <2
                level=[];
            end
            t=obj.ftrafo.detcoef_rangemax(level);
        end
    end
    
    %% filter
    
    methods
        
        function hardthresh1_inplace(obj,thresh)
            % filter all coefficients obj.C
            L=obj.ts.size(3);
            for j=1:L
                obj.ftrafo.C=obj.C{j};
                obj.ftrafo.hardthresh1_inplace(thresh);
                obj.C{j}=obj.ftrafo.C;
            end
        end
        
        function hardthresh2_inplace(obj,thresh1, thresh2)
            % filter with 2 thresholds thresh1<thresh2
            L=obj.ts.size(3);
            for j=1:L
                obj.ftrafo.C=obj.C{j};
                obj.ftrafo.hardthresh2_inplace(thresh1,thresh2);
                obj.C{j}=obj.ftrafo.C;
            end
        end
        
        function softthresh1_inplace(obj,thresh)
            % simple soft filter: reduce all coeff. obj.C by size
            L=obj.ts.size(3);
            for j=1:L
                obj.ftrafo.C=obj.C{j};
                obj.ftrafo.softthresh1_inplace(thresh);
                obj.C{j}=obj.ftrafo.C;
            end
        end
        
        
    end
    
    %% transforms
    methods
        
        function coeff=detcoef(obj,lvl,type)
            % coeff=~(level[,type]) extract detail coefficients
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            L=obj.ts.size(3);
            obj.ftrafo.C=obj.C{1};
            coeff0=obj.ftrafo.detcoef(lvl,type);
            coeff=zeros(size(coeff0,1),size(coeff0,2),L);
            coeff(:,:,1)=coeff0;
            for j=2:L
                obj.ftrafo.C=obj.C{j};
                coeff(:,:,j)=obj.ftrafo.detcoef(lvl,type);
            end
        end
        
        
        function coeff= appcoef(obj)
            % coeff=~() extract approximate coefficients from wavelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            L=obj.ts.size(3);
            obj.ftrafo.C=obj.C{1};
            coeff0=obj.ftrafo.appcoef();
            coeff=zeros(size(coeff0,1),size(coeff0,2),L);
            coeff(:,:,1)=coeff0;
            for j=2:L
                obj.ftrafo.C=obj.C{j};
                coeff(:,:,j)=obj.ftrafo.appcoef();
            end
        end
        
        function Cvec=C2vec(obj,cC)
            % convert transformation result to vector form
            % using C2vec of obj.ftrafo.
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            L=obj.ts.size(3);
            if nargin<2
                cC=obj.C;
            end
            cC1=obj.ftrafo.C2vec(cC{1});
            sn=numel(cC1);
            Cvec=zeros(L*numel(cC1),1);
            for j=1:L
                Cvec((j-1)*sn+1:j*sn)=obj.ftrafo.C2vec(cC{j});
            end
        end
        
        function cC=vec2C(obj,vC)
            % cC=~(vC) convert frame coefficients to vector form
            L=length(obj.C);
            cC=cell(L,1);
            for j=1:L
                cC{j}=obj.ftrafo.vec2C(vC((j-1)*L+1:j*L));
            end
        end
        
        function [mat, levelshown]=C2graph(obj, lvl)
            % mat=~() convert coefficient vector to format suitable for
            % graphical output
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || lvl<=obj.level,'local', 'smaller than max. #levels');
            
            if nargin < 2 || isempty(lvl)
                lvl=obj.ftrafo.level;
            end
            levelshown=lvl;
            L=obj.ts.size(3);
            img=obj.ftrafo.C2graph(lvl);
            mat=zeros(size(img,1),size(img,2),L);
            for j=1:L
                obj.ftrafo.C=obj.C{j};
                mat(:,:,j)= obj.ftrafo.C2graph(lvl);
            end
            
        end
        
        
    end
    
    
    %% statics
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=MultiResTransform3_Partial();
        end
        
    end
    
    %% graphics
    
    methods
        
        function new_titstr=modified_title(obj,z)
            h=get(gca,'Title');
            present_titstr=get(h,'String');
            if isempty(present_titstr)
                present_titstr=obj.ts.signalname;
            end
            addtit=obj.algo.name;
            if nargin >1
                addtit=[addtit,', z=',num2str(z)];
            end
            
            new_titstr=present_titstr;
            if iscell(present_titstr)
                new_titstr{1}=[addtit,', ',present_titstr{1}];
            else
                new_titstr=[addtit,', ',present_titstr];
            end
        end
        
        function graph_trafo(obj, open_new, z)
            % ~([open_new, z) show projection to z-plane of transform C
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if ~exist('z', 'var') || isempty(z)
                z=obj.zref;
            end
            if ~exist('open_new', 'var') || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            obj.ftrafo.C=obj.C{z};
            
            lev=obj.ftrafo.level();
            obj.ftrafo.graph_trafo(false,lev);
            
            title(obj.modified_title(z),'fontsize',12);
            
        end
        
        function show_trafo(obj)
            % show transform C  in open figure window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            suptitle(obj.modified_title(),14);
            
            lvl=obj.level;
            L=obj.ts.size(3);
            sd=factor_subplots(L);
            maxAbsC=zeros(L,1);
            
            for j=1:L
                subplot(sd(1),sd(2),j);
                obj.ftrafo.C=obj.C{j};
                M1=obj.ftrafo.C2graph(lvl);
                maxAbsC(j)=max(abs(M1(:)));
                
                imagesc(obj.repfun(M1));
                colormap('default');
                colorbar;
                if ~isempty(obj.repfun)
                    cblabel_str=func2str(obj.repfun);
                    cblabel(cblabel_str);
                end
                maxstr=num2str(maxAbsC(j),'%3.1e');
                if j==obj.zref
                    title(['z=',num2str(j),', MRA_{xy}=',obj.ftrafo.algo.name,...
                        ' (',obj.ftrafo.basisname,'), max=',maxstr],...
                        'fontsize',12);
                else
                    title(['z=',num2str(j),', ',obj.algo.name,', max=',...
                        maxstr],'fontsize',12);
                end
            end
            
        end
        
        function graph_sparsityDefect(obj,open_new,c)
            % ~() sparsity defect (lp-error of best s-term approx) of transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            if nargin<3
                c=10;
            end
            [sd,p]=obj.sparsityDefect(1/c);
            s=TimeSignal(sd);
            s.signalname={['rel. sparsity defect (c=',num2str(c),...
                ', norm p=',num2str(p),', z_0=',num2str(obj.zref),')'],...   
                obj.ts.get_signalname};
            s.marker='.';
            s.graph_signal(false);
            grid on;
            xlabel('frame (z)');   
            ylabel('\sigma_{1/c}(x)_p');
            
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 3d';
            ok= isa(obj.ts,'Signal3D');
            if ok
                descr='transform is 2d multiresolution';
                ok=isa(obj.ftrafo,'MultiResTransform2');
            end
        end
    end
    
end



