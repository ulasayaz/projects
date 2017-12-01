classdef IdTransform < FrameTrafo
    % identy frame transform (any dimension): stub class
    % to avoid handling the special case of missing transforms by if cases.
    % GT Stuttgart, 2014
    %
    % Example
    %
    % f=IdTransform();
    % res=f.test_framefeatures(); display(res);
    %
    
    properties
        dwtmode_active %@ (string) extension mode
    end
    
    
    %% constructor and commands
    methods
        
        function obj=IdTransform(signal)
            assert(nargin<1 || isa(signal,'SignalClass'),'needs SignalClass');
            % constructor
            if nargin==0
                signal=Signal2D();
            end
            obj = obj@FrameTrafo(signal);
            obj.dwtmode_active='bicubic';
        end
        
        function set_algo(obj)
            set_algo@FrameTrafo(obj);
            obj.algo.name='identity';
        end
        
        function set_dwtmode(obj, modestr)
            % ~(modestr) interpolation mode
            obj.requireOnly(ismember(modestr,{'bilinear','bicubic'}),'local',...
                'mode is admitted');
            obj.dwtmode_active=modestr;
        end
        
        function dec(obj)
            % decomposition (analysis) yields frame coeffs. obj.C
            % identity transform
            obj.C=obj.ts.xn;
        end
        
        function yn=rec(obj,cC)
            % signal reconstruction from frame coefficients C
            % identity transform
            yn=cC;
        end
        
        function yn=recSynthesis(obj,cC)
            % xn=~(C) reconstruct signal (vector) from frame coefficients C
            % is a left inverse of theta=obj.data'
            % identity transform
            yn=cC;
        end
        
        function yn= synthesize(obj,x)
            % yn=~(x) decompose x
            yn = x;
        end
        
        function x=recAnalysis(obj,yn)
            % x=~(yn) left inverse of analyze
            x=yn;
        end
        
        function xn= analyze(obj,cC)
            % yn=~(x) decompose x
            xn=cC;
        end
        
        
    end
    
    %% queries
    methods
        
        function bn= basisname(obj)
            % bn=~() basis name
            bn='natural basis';
        end
        
        function nv=frame_length(obj)
            % n=~() dimension of transform space
            nv=obj.ts.numel;
        end
        
        function nv=frame_dim(obj)
            % n=~() dimension of Hilbert space
            nv=obj.ts.numel;
        end
        
        
    end
    
    %% transforms
    methods
        
        function [signal2,BSE]=gridApprox(obj,c,p)
            % [signal2,BSE]=~(c,p) best grid approximation at compression c_s
            % signal2 ... signal approximation by coarsening(reduction+interpolation)
            % BSE ... lp-error of best s-term approximation (sparsity defect)
            obj.requireOnly(c>1,'local','compression rate >1');
            if ~exist('c','var') || isempty(c)
                c=16;
            end
            if ~exist('p','var') || isempty(p)
                p=2;  % lp-norm
            end
            method=obj.dwtmode_active;
            
            %signal2=obj.ts.clone();
            %signal2.coarsen(1/c,'bicubic');  % imresize introduces shift!
            
            signal2=obj.ts.make_like();
            signal2.signalname=[obj.ts.signalname,', interpolated (',method,...
                ', c=',num2str(c,'%3.1f'),')'];
            
            %% Sampling using nearest neighbour
            
            sz = obj.ts.N;
            
            indx = round(1:sqrt(c):sz(1));
            indy = round(1:sqrt(c):sz(2));
            Z =obj.ts.xn(indx,indy);
            % figure;imshow(uint8(Z));
            
            %% Simple Interpolation using Meshgrid and Interp2 fn.
            [X,Y] = meshgrid(indx,indy);
            [XI,YI] = meshgrid(1:sz(1),1:sz(2));
            if strcmp(method,'bilinear'),
                out = interp2(X,Y,Z,XI,YI,'linear',0); % performing bilinear interpolation
                % border handling
                out(:,sz(2)-ceil(sqrt(c))+2) = out(:,sz(2)-ceil(sqrt(c))+1);
                out(sz(2)-ceil(sqrt(c))+2,:) = out(sz(2)-ceil(sqrt(c))+1,:);
            elseif strcmp(method,'bicubic'),
                out = interp2(X,Y,Z,XI,YI,'cubic',0); % performing bicubic interpolation
                % border handling
                out(:,sz(2)-ceil(sqrt(c))+2) = out(:,sz(2)-ceil(sqrt(c))+1);
                out(sz(2)-ceil(sqrt(c))+2,:) = out(sz(2)-ceil(sqrt(c))+1,:);
            end
            signal2.replace_signal(out);
            
            if nargout>=2   % ||W-F(W)||/||W||
                BSE=sum(abs(obj.ts.xn(:)-signal2.xn(:)).^p)^(1/p);   % needs original values
            end
        end
        
        function [nodes,mask]=nodesPDF(obj,c, params)
            % nodes=~(c) nodes sampled with compression rate c
            obj.requireOnly(c>1,'local','compression rate in ]0,1[');
            obj.requireOnly(~obj.ts.isemptydata,'local','non-empty signal');
            DN = obj.ts.N; % data Size
            params.iter=1;
            p = 1/c; % undersampling factor
            % uniform density
            pdf= p*ones(DN);
            assert(abs(sum(pdf(:))-p*obj.ts.numel)<1,'definition of pdf');
            % generate sampling mask:
            mask = genSampling(pdf,params);
            % compute nodes:
            nodes=SampledNodes(find(mask),obj.ts.N,obj);
        end
        
        function nodes=nodesDefault(obj,c,params)
            % nodes=~() default nodes
            if nargin <2 || isempty(c)
                c=4;
            end
            if ~exist('params','var')
                params=struct;
            end
            nodes=obj.nodesPDF(c,params);
        end
        
    end
    
    %% graphics
    methods
        
        function cC=dec2graph(obj,cC)
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            cC=reshape(cC,obj.ts.N);
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal exists';
            ok= isa(obj.ts,'SignalClass');
        end
    end
    
    
    
end



