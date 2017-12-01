classdef FourierTrafo < FrameTrafo
    % N-dim. Fourier transform with graphical output;
    %
    % GT Stuttgart, 2014
    %
    
    properties (Access=public)
        fourierMethod  %@<struct> .str, .descr
    end
    
    properties (SetAccess=private)
        orderingActive
        orderingSchemes
        ok_license     %@ (bool) all necessary toolbox licenses exist
    end
    
    
    %% constructor and commands
    methods
        
        function obj=FourierTrafo(signal)
            % constructor
            obj.requireOnly(nargin==0 || isempty(signal) || isa(signal,'SignalClass'),...
                'local', ' signal data type ok');
            if nargin==0
                signal=SignalClass();
            end
            obj.ts=signal;
            obj.C=[];
            obj.repfun=@(F) log(1+abs(F));
            
            obj.fourierMethod.str='fft';
            obj.fourierMethod.descr='fast fourier transform';
            
            obj.set_algo;
            obj.fig=1;
            obj.use_repfun=true;
            obj.fontsize=12;
            
            obj.ok_license= true;
            
        end % constructor
        
        function set_algo(obj)
            obj.algo.version=0.9;
            obj.algo.versiondate='14.3.2014';
            obj.algo.name='Fourier';
            obj.algo.toolbox='matlab';
            obj.algo.mfilename=class(obj);
            obj.algo.author='G. Troll';
        end
        
        
        function set_0padding(obj,padd)
            % ~(padd) symm. 0-padding (left and right) by padd;
            
            obj.requireOnly(all(arrayfun(@obj.isnatural0,padd)),'local',...
                'padd is pair of non-negative integer');
            obj.ts.set_0padding(padd);
            
            obj.reset_Trafo;
            obj.ensureOnly(~obj.isTrafodone,'local','Fourier trafo reset');
            
        end
        
        function reset_Trafo(obj)
            % reset fourier trafo (e.g. after signal has changed)
            obj.C=[];
        end
        
        
        function dec(obj)
            % ~() decomposition whose inverse is rec, i.e. rec(dec)=Id
            % yielding an ONB (check with test_framefeatures() );
            % i.e. operation computeFrame by applying the linear dec operation
            % on the canonical basis (e_i) will yield an ONB (but is no
            % longer involution).
            obj.require(obj.ts.numel>0,'local','non-empty sample size');
            
            % n-dim-transform
            obj.C= fftn(obj.ts.xn);
            % normalize to get an orthogonal transform 
            obj.C=obj.C/sqrt(numel(obj.C));
            
        end
        
        function yn=rec(obj,cC)
            % ~(cC) signal reconstruction from frame coefficients C
            % should satisfy rec(dec)=Id
            
            % n-dim-transform
            yn=real(ifftn(cC));
            % normalize to get an orthogonal transform 
            yn=yn*sqrt(numel(yn));
            
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose vector x using an ONB; 
                       
            yn= fftn(reshape(x,obj.ts.N));
            yn=yn/sqrt(numel(yn));% to get an orthogonal transform
        end
        
        function xn= synthesize(obj,y)
            % yn=~(x) synthesize as inverse operation to analyze
                       
            xn= ifftn(reshape(y,obj.ts.N));
            xn=xn*sqrt(numel(xn)); % to get an orthogonal transform
        end
        
        
    end
    
    %% queries
    
    methods
        
        function bn= basisname(obj)
            % bn=~() basis name
            bn='fftn';
        end    
        
        function ok=isTrafodone(obj)
            % ok=~(): is result of a Fourier trafoi available?
            ok=~isempty(obj.C);
        end
        
        function nv=frame_length(obj)
            % n=~() number of elements of transform
            nv=obj.ts.numel;
        end
        
        function nv=frame_dim(obj)
            % n=~() dimension of Hilbert space
            nv=obj.ts.numel;
        end
        
        
        
        %% transforms
        function mat=C2graph(obj)
            % mat=~() convert coefficient vector to image
            % will be redefined in child classes
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            mat=obj.C;
        end
        
        function cC=vec2C(obj,vC)
            % convert vector to form usable for reconstruction rec
            cC=reshape(vC,obj.ts.N);
        end
        
        function [nodes,mask]=nodesPDF(obj,c,params)
            % nodes=~(c) nodes sampled with compression rate c
            obj.requireOnly(isscalar(c) && c>1,'local','compression rate(s) in ]0,1[');
            obj.requireOnly(~obj.ts.isemptydata,'local','non-empty signal');
            if nargin<3 || ~isfield(params,'P')
                params.deg=5; % Variable density polymonial degree
            end
            if ~isfield(params,'iter')
                params.iter=10;
            end
            params.box=[-1,1];  % centered coordinates (origin is center)
            DN = obj.ts.N; % data Size
            p = 1/c;% [0.25];% undersampling factor
            
            nodes=SampledNodes([],obj.ts.N,obj);
            
            % generates the sampling PDF concentrated at center
            pdf = genPDF(DN,p,params);
            
            % generate sampling mask:
            mask = genSampling(pdf,params);
            % shift back to non-centered Fourier coord:
            mask = fftshift(mask);
            % compute nodes:
            nodes.set_data(find(mask));
        end
        
        function nodes=nodesDefault(obj,p,params)
            % nodes=~() default nodes
            if nargin <2 || isempty(p)
                p=4;
            end
            if ~exist('params','var')
                params=struct;
            end
            nodes=obj.nodesPDF(p,params);
        end
        
        
    end
    
    methods (Hidden)
        
        function v=xn(obj)
            % signal vector
            v=obj.ts.xn;
        end
        
    end
    
    methods
        %% coordinate transformations
        
        
        
    end
    
    
    %% tests
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=FourierTrafo();
        end
        
    end
    
    %% graphics
    methods
        
        function tit2= title2(obj)
            % second title line
            tit2=[];
            if obj.padd>0
                tit2=['0-padding: ',num2str(obj.padd)];
            end
            if ~isempty(tit2)
                tit2=[tit2,', '];
            end
            tit2=[tit2, 'winfun=',obj.ts.winfun2str];
        end
        
        function cC=dec2graph(obj,cC)
            % post process decomposition obj.dec for graphical output
            ss=obj.ts.N;
            if numel(cC)==prod(ss)
                cC= fftshift(reshape(cC,ss));  % e.g. Fourier3 needs complete ss
            elseif numel(cC)==ss(1)*ss(2)
                cC= fftshift(reshape(cC,ss(1:2)));               
            end
        end
        
        
        function graph_distribution(obj, open_new)
            % ~() show distribution of frame coefficients in open window
            obj.require(~isempty(obj.C),'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            img_dwt=abs(obj.C2graph());
            HM=HMatrix(img_dwt);
            HM.graph_distribution(false);
            
            h=get(gca,'Title');
            present_titstr=get(h,'String');
            addtit=['distribution of moduli of coefficients (',obj.basisname,')'];
            new_titstr=present_titstr;
            new_titstr{1}=addtit;
            if ~iscell(present_titstr)
                new_titstr{2}=present_titstr;
            end
            title(new_titstr,'fontsize',12);
            
        end
        
        function graph_winfun(obj, open_new)
            % show window function winfun in the open window
            obj.requireOnly(isa(obj.ts,'Signal2D'),'local','signal exists');
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo);
            end
            obj.ts.graph_winfun;
        end
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is SignalClass';
            ok= isa(obj.ts,'SignalClass');
        end
    end
    
end




