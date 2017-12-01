classdef SampleTrafo < FrameTrafo
    % sampling of an orthogonal basis transform (cf. invariant)
    % given in the form of 1 implicit operator (ftrafo).
    % This defines the encoder of a compressive sensing situation.
    % GT Stuttgart, 2014
    %
    % example:
    % signal=Signal2D.make_fromImage('cameraman.bmp');
    % signal.resize(0.25);
    % B=IdTransform(signal);
    % c=4; % compression rate
    % params.P=0; params.iter=1;
    % nodes=B.nodesPDF(c, params);
    % f=SampleTrafo(signal);
    % f.set_transform(B,nodes);
    % [Xi,Yi]=meshgrid(1:f.ts.size(2), 1:f.ts.size(1));
    % z=f.interp(Xi,Yi);
    %
    
    properties
        
        ftrafo   %@<FrameTrafo> (encoder) orthogonal basis transform to be sampled
        
    end
    
    properties (SetAccess=protected)
        
        nodes   %@<SampledNodes> nodes of sample
        L       %@<integer> #sample lines
        
    end
    
    %% constructor and commands
    methods
        
        function obj=SampleTrafo(signal)
            % constructor obj=~(signal)
            assert(nargin<1 || isa(signal,'SignalClass'),'needs SignalClass');
            if nargin<1
                signal=[];
            end
            obj = obj@FrameTrafo(signal);
        end
        
        function set_transform(obj, hftrafo, hnodes)
            % ~(hnodes) set nodes which are to be sampled
            obj.requireOnly(~isempty(hnodes),'local','nodes exist');
            obj.requireOnly(isempty(hftrafo) || isa(hftrafo,'FrameTrafo'),...
                'local','needs a FrameTrafo');
            obj.requireOnly(~isempty(obj.ts),'local','signal is set');
            obj.ftrafo=hftrafo;
            obj.nodes=hnodes;
            obj.ftrafo.ts=obj.ts;
            obj.dec;
            if isa(hftrafo,'IdTransform')
                obj.colormap_active=obj.ts.colormap_active;
            else
                obj.colormap_active='default';
            end
        end
        
        function set_signal(obj,hts)
            % ~(hts) set time signal to hts of class SignalClass
            obj.requireOnly(isa(hts,'SignalClass'),'local','ts is a 2d signal');
            obj.ts=hts;
            obj.ftrafo.ts=hts;
            obj.reset_Trafo;
            % do not check invariant because caller may be constructor
            obj.ensureOnly(~obj.isTrafodone,'local', ...
                'decomposition is reset');
        end
        
        function reset_Trafo(obj)
            % reset trafo (e.g. after signal has changed)
            if ~isempty(obj.ftrafo)
                if isempty(obj.nodes) || max(obj.nodes.get(1))>numel(obj.ftrafo.C2vec)
                    obj.set_nodesDefault();
                end
                obj.ftrafo.C=[];
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
        
        function nodes=nodesDefault(obj,p,params)
            % nodes=~(p,[params]) default nodes
            obj.requireOnly(~isempty(obj.ftrafo),'local','transform exists');
            if nargin <2 || isempty(p)
                p=[];
            end
            if ~exist('params','var')
                params=struct;
            end
            nodes=obj.ftrafo.nodesDefault(p,params);
        end
        
        function set_nodes(obj, hnodes)
            % ~(hnodes) set nodes which are to be sampled
            obj.requireOnly(isa(hnodes,'SampledNodes')','local','nodes type ok');
            obj.nodes=hnodes;
        end
        
        function dec(obj)
            % decomposition (analysis) yields frame coeffs. obj.C
            obj.requireOnly(~isempty(obj.ftrafo),'local', 'ftrafo set');
            obj.ftrafo.dec;
            obj.C=obj.sample_nodes();
        end
        
        function yn=rec(obj,cC)
            % signal reconstruction from frame coefficients C
            obj.requireOnly(~isempty(obj.ftrafo),'local', 'ftrafo set');
            X=obj.ftrafo.embedSample(obj.nodes,cC);
            yn=obj.ftrafo.rec(X);
        end
        
        function yn=recSynthesis(obj,cC)
            % xn=~(C) reconstruct signal (vector) from frame coefficients C
            % is a left inverse of theta=obj.data'
            yn=obj.analyze(cC);
        end
        
        function yn= synthesize(obj,x)
            % yn=~(x) decompose x and sample at nodes
            % SampleTrafo.synthesize uses dec instead of rec
            x = reshape(x,obj.ftrafo.N);
            obj.ftrafo.set_content(x);
            obj.ftrafo.dec;
            yn=obj.sample_nodes();
        end
        
        function x=recAnalysis(obj,yn)
            % x=~(yn) left inverse of analyze
            x=obj.synthesize(yn);
        end
        
        function xn= analyze(obj,cC)
            % yn=~(x) decompose x and sample at nodes
            % SampleTrafo.analyze uses rec instead of dec
            
            % embed samples into image of correct size
            X=obj.ftrafo.embedSample(obj.nodes,cC);
            % rescontruct
            xn=real(obj.ftrafo.synthesize(X));
        end
        
        function yn=sample_nodes(obj)
            % Y=~() measure transform ftrafo at at nodes
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(~isempty(obj.nodes),'local','nodes exist');
            %cC=obj.ftrafo.C2vec;
            %assert(max(obj.nodes)<=numel(cC) && min(obj.nodes)>=1,...
            %    'nodes within range');
            % yn=cC(obj.nodes);
            yn=obj.ftrafo.sample(obj.nodes);
            obj.ensureOnly(isvector(yn),'local', 'returns vector');
        end
        
        function y=sim_measurement(obj)
            % y=~([L]) simulate a measurement vector y by sampling nodes
            
            obj.requireOnly(~isempty(obj.ts),'local','needs signal to process');
            obj.requireOnly(~isempty(obj.nodes),'local','needs nodes to sample');
            
            obj.dec;  % apply sampling transform
            y=obj.sample_nodes;  % sample nodes in in image of sampling transform
            if ~isempty(obj.SNR) && obj.SNR>0
                % add SNR [dB] noise to measurement
                obj.sigma = std(y(:))*10^(-obj.SNR/20);
                y=y+obj.sigma*randn(size(y));
            end
        end
        
        function ss=measurement2graph(obj,yn)
            % ss=~(yn) transforms measurement vector to output form, i.e.
            % in this class embeds it in the original signal
            % yn ... measurement vector, e.g. output of sim_measurement
            obj.requireOnly(~isempty(obj.ts),'local','needs signal to process');
            obj.requireOnly(~isempty(obj.nodes),'local','needs nodes to sample');
            obj.requireOnly(isvector(yn) && numel(obj.nodes)==length(yn),...
                'local', 'yn is measurement vector');
            
            M=obj.ftrafo.embedSample(obj.nodes,yn);
            M(M==0)=NaN;
            %             M=NaN*zeros(obj.ts.N);
            %             M(obj.nodes)=yn;
            
            M=obj.ftrafo.dec2graph(M);  % post processing e.g. shift for fft
            ss=obj.ts.make();
            ss.set_signal(M);
            ss.set_automaticStyle;
            ss.repfun=obj.ftrafo.repfun;
            ss.colormap_active=obj.colormap_active;
            
        end
        
        
    end
    
    %% queries
    methods
        
        function ok=isTrafodone(obj)
            % ok=~(): is frame decomposition of signal available?
            ok=~isempty(obj.ftrafo) && ~isempty(obj.ftrafo.C);
        end
        
        function F=get_frame(obj)
            F=obj.ftrafo.frame;
        end
        
        function bn= basisname(obj)
            % bn=~() basis name
            bn=obj.ftrafo.basisname;
        end
        
        function str=transformName(obj)
            % name of frame transform
            str=['Projections (',obj.basisname,')'];
        end
        
        function nv=frame_length(obj)
            % n=~() dimension of transform space
            % must be redefined for some subclasses
            nv=obj.ftrafo.frame_length;
        end
        
        function nv=frame_norm(obj)
            % l2-norm of transform coefficients
            % tests in this form only isometries (e.g. ONS)
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            nv=obj.ftrafo.frame_norm;
        end
        
        function nv=frame_dim(obj)
            % n=~() dimension of Hilbert space
            % must be redefined for some subclasses
            nv=obj.nodes.numel;
        end
    end
    
    %% filters
    methods
        
        function val=interp(obj,Xi,Yi,Zi, method)
            % [signal2,BSE]=~(Xi,Yi,Zi) interpolation sampled signal
            % e.g. [Xi,Yi]=meshgrid(1:obj.ts.size(2), 1:obj.ts.size(1));
            % val ... interpolation values at (Xi,Yi)
            obj.requireOnly(~isempty(obj.nodes),'local','nodes are set');
            if ~exist('method','var') || isempty(method)
                method='linear';
            end
            
            N=squeeze(obj.ts.N);
            z=obj.ts.xn(obj.nodes.get(1));
            if length(N)==3
                [nodesY,nodesX,nodesZ]=ind2sub(obj.ts.size,obj.nodes.get(1));
                F=TriScatteredInterp(nodesX,nodesY,nodesZ, z,method);
                val=F(Xi,Yi,Zi);
            elseif length(N)==2
                [nodesY,nodesX]=ind2sub(obj.ts.size,obj.nodes.get(1));
                F=TriScatteredInterp(nodesX,nodesY,z,method);
                val=F(Xi,Yi);
            elseif length(N)==1
                F=TriScatteredInterp(obj.nodes.get(1),z,method);
                val=F(Xi);
            end
            
            
        end
        
    end
    
    %% static methods
    methods (Static)
        
    end
    
    %% graphics
    methods
        
        function graph_trafo(obj, open_new)
            % ~() show frame of transform (computed by computeFrame) in open window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            obj.ftrafo.graph_trafo(false);
            h=get(gca,'Title');
            present_titstr=get(h,'String');
            addtit='Encoder';
            new_titstr{1}=addtit;
            if ~iscell(present_titstr)
                new_titstr{2}=present_titstr;
            elseif length(present_titstr)>1
                new_titstr{2}=present_titstr{2};
            end
            title(new_titstr,'fontsize',12);
            
        end
        
        function graph_distribution(obj, open_new, yn)
            % ~() show sampling projections of encoding transform
            obj.requireOnly(~isempty(obj.nodes),'local', 'needs nodes (e.g. via sim_measurement');
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2
                open_new=true;
            end
            if nargin <3
                yn=[];  % optional measurement vector
            end
            if open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            if isempty(yn)
                ts2=obj.ts;
                %  SampleTrafo.synthesize uses dec instead of rec
                % so yn contains the coefficients of the forward transform
                yn= obj.synthesize(obj.ts.xn);  % sideeffect!
                obj.ts=ts2;
            end
            
            ss=obj.measurement2graph(yn);
            ss.colormap_freeze=~open_new;
            ss.graph_signal(false);
            %imagescnan(obj.ftrafo.repfun(M));axis image;rmaxis
            %cblabel(func2str(obj.ftrafo.repfun));
            
            c= obj.ts.numel/numel(yn);
            nf=obj.nodes.frameNumberUnique;
            
            SNR_str=[num2str(obj.SNR,'%3.1f'),' dB noise'];
            if isempty(obj.SNR)
                SNR_str='no noise';
            end
            title({['Projections ',obj.ftrafo.basisname],...
                ['#node frames=',num2str(nf),...
                ', c=',num2str(c,'%3.2f'),...
                ', ', SNR_str]}, 'fontsize', 12);
            
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal exists';
            ok= ~isempty(obj.ts) && isnumeric(obj.img); % && ~isvector(obj.img);
            if ok
                descr='frame transform exists';
                ok=~isempty(obj.ftrafo);
            end
            if ok
                descr='frame transform is (orthonormal) basis transform';
                % check dimension of transform result
                ok=isempty(obj.ftrafo.C) || isequal(numel(obj.ftrafo.C),obj.ts.numel);
            end
        end
    end
    
end



