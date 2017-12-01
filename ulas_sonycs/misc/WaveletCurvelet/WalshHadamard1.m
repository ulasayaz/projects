classdef WalshHadamard1 < FrameTrafo
    % 1-dim. Walsh-Hadamard transformations with graphical output;
    %
    % Example:
    % Example
    % 1.) define a signal, e.g.
    %     signal=TimeSignal.make_sinus(3);
    % 2.) create object
    %     f=WalshHadamard1(signal);
    %
    % 3.) apply transform
    %     f.dec;
    % 4.) show result
    %     f.graph_trafo;
    %
    % Test frame features:
    % res=f.test_framefeatures(); display(res);
    %
    % GT Stuttgart, 2014
    %
    
    properties (Access=public)
        
    end
    
    properties (SetAccess=private)
        orderingActive %@<integer> ordering of frame vectors
        orderingSchemes %@<struct> available ordering schemes
        ok_license     %<logical> all necessary toolbox licenses exist
    end
    
    
    %% constructor and commands
    methods
        
        function obj=WalshHadamard1(signal)
            % constructor
            obj.requireOnly(nargin==0 || isempty(signal) || isa(signal,'SignalClass'),...
                'local', ' signal belongs to class SignalClass');
            if nargin==0
                signal=TimeSignal();
            end
            obj.ts=signal;
            obj.C=[];
            obj.repfun=@(F) F;
            
            obj.set_algo;
            obj.fig=1;
            obj.use_repfun=true;
            obj.fontsize=12;
            obj.set_orderingSchemes();
            obj.orderingActive=1;
            
            obj.ok_license= license('test','signal_toolbox');
            
            obj.ensure(isa(obj.ts,'TimeSignal'),'local','2D signal created');
            
        end % constructor
        
        function set_algo(obj)
            obj.algo.version=0.9;
            obj.algo.versiondate='14.3.2014';
            obj.algo.name='Walsh-Hadamard1';
            obj.algo.toolbox='matlab';
            obj.algo.mfilename=class(obj);
            obj.algo.author='G. Troll';
        end
        
        function set_orderingActive(obj,id)
            % set ordering scheme
            obj.requireOnly(nargin>1 && id>=1 && id <=length(obj.orderingSchemes),...
                'local', 'adimissible value');
            obj.orderingActive=id;
        end
        
        function set_orderingSchemes(obj)
            % ~() set data of ordering schemes for Walsh-Hadamard basis
            % the order by sequency orders the basis w.r.t to increasing resolution,
            % i.e. frequency. 
            field1 = 'id';
            value1 = {1,2,3};
            
            field2 = 'name';
            value2 = {'sequency','Hadamard','dyadic'};
            
            
            obj.orderingSchemes=struct(field1,value1,field2,value2);
            
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
            
            n=2^obj.ts.size_dyadic;
            % 1d-transform
            obj.C= fwht(obj.ts.xn,n,obj.orderingSchemes(obj.orderingActive).name);
            % normalize to get an orthogonal transform (otherwise
            % involution)
            obj.C=obj.C*sqrt(numel(obj.C));
            
        end
        
        function yn=rec(obj,cC)
            % ~(cC) signal reconstruction from frame coefficients C
            % should satisfy rec(dec)=Id
            
            n=2^obj.ts.size_dyadic;
            
            % 1d-transform
            yn=ifwht(cC,n,obj.orderingSchemes(obj.orderingActive).name);
            % normalize to get an orthogonal transform (otherwise
            % involution)
            yn=yn/sqrt(numel(yn));
            
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x using an ONB
            n=2^obj.ts.size_dyadic;
            
            % 21-transform
            yn= fwht(x,n,obj.orderingSchemes(obj.orderingActive).name);
            yn=yn*sqrt(numel(yn));% to get an orthogonal transform
        end
        
        function xn= synthesize(obj,y)
            % yn=~(x) synthesize as inverse operation to analyze
            n=2^obj.ts.size_dyadic;
            
            % 1d-transform
            xn= ifwht(y,n(1),obj.orderingSchemes(obj.orderingActive).name);
            xn=xn/sqrt(numel(xn)); % to get an orthogonal transform
        end
        
        
    end
    
    %% queries
    
    methods
        
        function bn= basisname(obj)
            % bn=~() basis name
            bn=['fwht (ordering= ',obj.orderingSchemes(obj.orderingActive).name,')'];
        end
        
        function h=get_orderingActive(obj)
            % h=~() get ordering scheme
            h=obj.orderingSchemes(obj.orderingActive);
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
        
        function nodes=nodesDefault(obj,c,params)
            if nargin <2 || isempty(c)
                c=4;
            end
            if nargin <3
                params=struct;
            end
            try
                nodes=nodesExpRnd(obj,c);
            catch                
                nodes=obj.nodesPDF(c,params);
            end
        end
        
        function [nodes,mask]=nodesPDF(obj,c,params)
            % nodes=~(c) nodes sampled with compression rate c
            obj.requireOnly(c>1,'local','compression rate in ]0,1[');
            obj.requireOnly(~obj.ts.isemptydata,'local','non-empty signal');
            if nargin <3
                params=struct;
            end
            params.box=[0,1];  % origin is left (not center)
            [nodes,mask]=nodesPDF@FrameTrafo(obj,c,params);
        end
        
        function nodes=nodesExpRnd(obj,c)
            % nodes=~(c) compute random sample of nodes which is
            % exponentially distributed with highest density at index
            % origin; 
            obj.require(~obj.ts.isemptydata,'local',' signal needs data');
            obj.requireOnly(c>=1,'local','compression rate >1');
            
            L=obj.ts.length; 
            K=L/c;
            mu=L/10; 
            tol=0.01;
            ctol=max(1,tol*L);            
            itmax0=20;
            d=Inf;
            it=0;
            f=1;
            while d > ctol && it<=itmax0
                f=f*2;
                nodes=unique(max(1,min(L,round(exprnd(mu,[1,f*ceil(K)])))));
                Ln=min(length(nodes),round(K));
                nodes=nodes(1:Ln);
                d=abs(Ln-K);
            end
            nodes=SampledNodes(nodes,obj.ts.N,obj);
            obj.ensureOnly(max(nodes.data)<=numel(obj.ts.xn) && min(nodes.data)>=1,'local',...
                'nodes within range');
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
            w=WalshHadamard1();
        end
        
    end
    
    %% graphics
    methods
        
        function tit2= title2(obj)
            % second title line
            tit2=['ordering ',obj.orderingSchemes(obj.orderingActive).name];
        end
        
        function graph_trafo(obj,open_new)
            %  show coefficients of transformed signal in open figure window
            obj.require(obj.isTrafodone,'local','transform result is set');
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo);
            end
            
            cC=obj.C;
            if obj.use_repfun
                cC= obj.repfun(cC);
                cC(cC<-10)=-10;
            end
            plot(obj.dec2graph(cC)) ;
            
            tittext=['Transformation: ',obj.basisname,', Signal: ',...
                obj.ts.signalname];
            title(tittext,'fontsize',12);
            
        end
        
        function graph_distribution(obj, open_new)
            % ~() show distribution of frame coefficients in open window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            cC=obj.C2vec();
            HM=HMatrix(cC);
            HM.graph_distribution(false);
            
            h=get(gca,'Title');
            present_titstr=get(h,'String');
            addtit=['distribution of coefficients: ',...
                obj.basisname];
            new_titstr{1}=addtit;
            if ~iscell(present_titstr)
                new_titstr{2}=present_titstr;
            elseif length(present_titstr)>1
                new_titstr{2}=present_titstr{2};
            end
            title(new_titstr,'fontsize',12);
            
        end
        
        
    end
    
    %% invariant
    methods
       function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is TimeSignal';
            ok= isa(obj.ts,'TimeSignal');  
        end
    end
    
end



