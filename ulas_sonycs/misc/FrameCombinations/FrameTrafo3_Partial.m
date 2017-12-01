classdef FrameTrafo3_Partial<FrameTrafo
    % 3d partial frame transform
    % applying a 2d frame transforms to each z-plane separately;
    % used e.g. as one of the two transforms in FrameTrafoPair.
    %
    % Example
    % ========
    %     w=FrameTrafo3_Partial();
    %     signal2=Signal2D.make_fromImage('cameraman.bmp');
    %     L=3; signal3=Signal3D.make_CyclefromLowerDimSig(signal2,1,L);
    %     signal3.play_signal;
    %     w.set_signal(signal3);
    %     w.set_transform(WalshHadamard2());
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
        ftrafo    %@<FrameTrafo>  basic 2d transform
    end
    
    
    %% commands
    methods
        
        function obj=FrameTrafo3_Partial(signal,hftrafo)
            % constructor
            if nargin==0
                signal=Signal3D();
            end
            obj = obj@FrameTrafo(signal);
            if nargin<2
                hftrafo=FrameTrafo();
            end
            obj.ftrafo=hftrafo;
            obj.repfun=hftrafo.repfun;
            if ~obj.ts.isemptydata
                obj.frameDepth=obj.ts.size(3);
                obj.ftrafo.ts=Signal2D(signal.xn(:,:,1));
            end
            obj.zref=max(1,ceil(signal.size(3)/2));  % mid point z-plane
            obj.anweight1=0.5;
            
        end
        
        function set_algo(obj)
            set_algo@FrameTrafo(obj);
            obj.algo.name='partial 3d';
        end
        
        function set_signal(obj,hts)
            % ~(hts) set time signal to hts of class SignalClass
            obj.requireOnly(isa(hts,'Signal3D'),'local','ts is a 3d-signal');
            set_signal@FrameTrafo(obj,hts);
            if ~isempty(hts.xn)
                obj.ftrafo.set_signal(Signal2D(hts.xn(:,:,1)));
            end
            obj.zref=max(1,ceil(hts.size(3)/2));  % mid point z-plane
        end
        
        function set_transform(obj, hftrafo1)
            % ~(hftrafo) set basic 2d transform
            obj.requireOnly(isempty(hftrafo1) || isa(hftrafo1,'FrameTrafo'),...
                'local','needs 2d muliresolutin transform');
            obj.ftrafo=hftrafo1;
            obj.repfun=hftrafo1.repfun;
            obj.reset_Trafo();
            
            obj.colormap_active='default';
        end
        
        
        function dec(obj)
            % partial discrete multiresolution analysis
            % apply 2d transform to each z-projection
            obj.requireOnly(obj.ts.numel>0,'local','non-empty sample size');
            
            obj.ftrafo.ts.set_signal(obj.ts.xn(:,:,1));
            
            L=obj.ts.size(3);
            obj.frameDepth=L;
            
            obj.ftrafo.ts.set_signal(obj.ts.xn(:,:,1));
            obj.ftrafo.dec();
            
            % alternative to object array would be cell array
            %             constructor = str2func(class(obj.ftrafo.C));
            %             obj.C(L)=constructor;
            obj.C=cell(L,1);
            obj.C{1}=obj.ftrafo.C;
            
            for j=2:L
                obj.ftrafo.ts.set_signal(obj.ts.xn(:,:,j));
                obj.ftrafo.dec();
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
        
        function bn= basisname(obj)
            % bn=~() basis name
            bn=obj.ftrafo.basisname();
        end
        
        function ok=isTrafodone(obj)
            % ok=~(): is result of DWT available?
            ok=~isempty(obj.C) && ~isempty(obj.C{1});
        end
        
        function nv=frame_length(obj)
            % n=~() dimension of transform space
            % must be redefined for some subclasses
            nv=obj.ts.size(3)*obj.ftrafo.frame_length;
        end
        
        function nv=frame_dim(obj)
            % n=~() dimension of Hilbert space
            % must be redefined for some subclasses
            nv=obj.ts.size(3)*obj.ftrafo.frame_dim;
        end
        
        function r=frame_redundancy(obj)
            % redundancy of the frame used (1 for orthogonal wavelets
            % coeff, max. ~7.8 for 2nd gen. curvelets)
            r=obj.ftrafo.frame_length/obj.frame_dim;
        end
        
        function nv=frame_norm(obj)
            % l2-norm of transform coefficients
            % tests in this form only isometries (e.g. ONS)
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            nv=obj.ftrafo.frame_norm;
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
            nc=L*obj.ftrafo.numel(cC1);
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
        
        function set_nodesDefault(obj,p,params)
            % ~(p,[params]) set default nodes
            obj.requireOnly(~isempty(obj.ftrafo),'local','transform exists');
            if nargin <2 || isempty(p)
                p=[];
            end
            if ~exist('params','var')
                params=struct;
            end
            obj.nodes=obj.nodesDefault(p,params);
        end
        
        function [nodes,c,cmean]=nodesDefault(obj,c,params)
            % [nodes,c,cmean]=~(c,[params]) default nodes using compression rate(s) c
            % if c is a vector and any if its elements is zero, c will
            % be computed s.t. its harmonic mean will be 1/params.p0
            % and the c(each side frame)=1/params.k0 * c(reference frame).
            
            obj.requireOnly(~isempty(obj.ftrafo),'local','transform exists');
            obj.requireOnly(all(c>1) || all(c==0),'local','compression rates in ]1,Inf[');
            obj.requireOnly(isscalar(c) || length(c)==obj.ts.size(3),'local','c scalar or length is #frames');
            if nargin <2 
                c=[];
            end
            if ~exist('params','var')
                params=struct;
            end
            if ~isfield(params,'p0')
                params.p0=1/8;  % mean undersampling rate
            end
            if ~isfield(params,'k0')
                params.k0=1/4;  % weight of side frames compared to central frame
            end
            if ~isfield(params,'frameNumber')
                params.frameNumber=1;
            end
            dim=obj.signal_dim;
            if ~isfield(params,'zref') || dim <3
                if dim==3
                    params.zref=max(1,ceil(obj.ts.size(3)/2)); % central frame
                else
                    params.zref=1;
                end
            end
            p=1./c;
            if length(p)>params.frameNumber
                params.frameNumber=length(p);
            end
            Lc=params.frameNumber;
            if any(~isfinite(p))  % compute default values for probabilities for each frame
                % c0=1/p0 is harmonic mean: 1/c0=(L-1)/L * k0/c1 + 1/L * 1/c1            
                k0=params.k0;
                p0=params.p0;
                p=k0*p0*Lc/(k0*(Lc-1)+1) * ones(1,Lc);
                p(params.zref)=p(params.zref)/k0;
            end
            
            c=1./p;   
            cmean=1/mean(p);  % harmonic mean
            if Lc==1
                nodes=obj.ftrafo.nodesDefault(c,params);
                nodes.fL=obj.ts.size(3);
            else  % separate nodes for each frame
                nodes= SampledNodes([],obj.ts.N,obj.ftrafo);
                nodes.set_data(cell(Lc,1));
                for j=1:Lc
                    hnodes=obj.ftrafo.nodesDefault(c(j),params);
                    nodes.set_cell(j,hnodes.data);
                end
            end
            
        end % nodesDefault
        
        
        function yn=sample(obj,nodes)
            % Y=~() measure transform ftrafo at at nodes
            % each z-plane will be sampled at the same nodes.
            cC=obj.C2vec;
            L=obj.ts.size(3);
            LN=nodes.numel;
            yn=zeros(LN,1);
            LC=length(cC)/L;
            offset=0;
            
            for j=1:L   % sample each z-plane of cC:
                cCj=cC((j-1)*LC+1:j*LC);
                Lj=nodes.numel(j);
                yn(offset+1: offset+Lj)=cCj(nodes.get(j));
                offset=offset+Lj;
            end
            
        end
        
        function cC=embedSample(obj,nodes,yn)
            % Y=~(nodes,cC) embed samples yn into frame coefficients;
            % ~ is right inverse of sample;
            cC=zeros(obj.frame_length,1);
            % set values at nodes
            L=obj.ts.size(3);
            
            % assume that even for separateNodesForEachFrame
            % each cell has the same length
            % assert(length(unique(cellfun(@numel,nodes.data)))==1;
            LN=nodes.numel;
            assert(length(yn)==LN,'yn has right length');
            LC=length(cC)/L;
            offset=0;
            
            for j=1:L   % sample each z-plane of cC:
                cCj=zeros(LC,1);
                Lj=nodes.numel(j);
                cCj(nodes.get(j))=yn(offset+1:offset+Lj);
                cC((j-1)*LC+1:j*LC)=cCj;
                offset=offset+Lj;
            end
        end
        
        
    end
    
    methods (Hidden)
        
        
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
        
        function cC=dec2graph(obj,cC)
            % post process decomposition obj.dec for graphical output
            cC=reshape(cC,obj.ts.N);
            L=obj.ts.size(3);
            for z=1:L
                cC(:,:,z)= obj.ftrafo.dec2graph(cC(:,:,z));
            end
        end
        
        function [mat, levelshown]=C2graph(obj, lvl)
            % mat=~([lvl]) convert coefficient vector to format suitable for
            % graphical output
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            
            if nargin < 2
                lvl=[];
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
        
        
    end
    
    %% graphics
    
    methods
        
        function new_titstr=modified_title(obj,z)
            h=get(gca,'Title');
            present_titstr=get(h,'String');
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
            obj.ftrafo.graph_trafo(false);
            
            title(obj.modified_title(z),'fontsize',12);
            
        end
        
        function show_trafo(obj)
            % show transform C  in open figure window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            suptitle(obj.modified_title(),14);
            
            L=obj.ts.size(3);
            sd=factor_subplots(L);
            
            for j=1:L
                subplot(sd(1),sd(2),j);
                obj.ftrafo.C=obj.C{j};
                imagesc(obj.repfun(obj.ftrafo.C2graph()));
                colormap('default');
                colorbar;
                if ~isempty(obj.repfun)
                    cblabel_str=func2str(obj.repfun);
                    cblabel(cblabel_str);
                end
                if j==obj.zref
                    title(['z=',num2str(j),', MRA_{xy}=',obj.ftrafo.algo.name,...
                        ' (',obj.ftrafo.basisname,')'],...
                        'fontsize',12);
                else
                    title(['z=',num2str(j),', ',obj.algo.name],'fontsize',12);
                end
            end
            
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 3d';
            ok= isa(obj.ts,'Signal3D');
            if ok
                descr='client transform is FrameTrafo too';
                ok=isa(obj.ftrafo,'FrameTrafo');
            end
        end
    end
    
end


