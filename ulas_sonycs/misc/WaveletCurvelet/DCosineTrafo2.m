classdef DCosineTrafo2 < FrameTrafo
    % 2-dim. discrete cosine transformations with graphical output;
    %
    % Example:
    % Example
    % 1.) define a signal, e.g.
    %     signal=Signal2D.make_fromImage('Mondrian.tif');
    %     signal=Signal2D.make_fromImage('cameraman.bmp');
    %     signal=Signal2D.make_zernike(48);
    % 2.) create object
    %     f=DCosineTrafo2(signal);
    %
    % 3.) apply transform
    %     f.dec;
    % 4.) show result
    %     f.graph_trafo;
    %
    % Test frame features:
    %     res=f.test_framefeatures(); display(res);   
    %    [ok,err1,err2]= f.test_DecRec(); disp(ok);
    %    [ok,err1,err2]= f.test_AnalyzeSynthsize(); disp(ok);
    %
    % GT Stuttgart, 2014
    %
    
    properties (Access=public)
        
    end
    
    properties (SetAccess=private)       
        ok_license     %<logical> all necessary toolbox licenses exist
    end
    
    
    %% constructor and commands
    methods
        
        function obj=DCosineTrafo2(signal)
            % constructor
            obj.requireOnly(nargin==0 || isempty(signal) || isa(signal,'Signal2D'),...
                'local', ' signal is 2D');
            if nargin==0
                signal=Signal2D();
            end
            obj.ts=signal;
            obj.C=[];
            obj.repfun=@(F) log(abs(F));
            
            obj.set_algo;
            obj.fig=1;
            obj.use_repfun=true;
            obj.fontsize=12;           
            
            obj.ok_license= license('test','signal_toolbox');
            
            obj.ensure(isa(obj.ts,'Signal2D'),'local','2D signal created');
            
        end % constructor
        
        function set_algo(obj)
            obj.algo.version=0.9;
            obj.algo.versiondate='14.3.2014';
            obj.algo.name='discr. cosine';
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
                       
            obj.C= dct2(obj.ts.xn);  
            % normalize to get an orthogonal transform 
            obj.C=obj.C/sqrt(numel(obj.C));            
            
        end
        
        function yn=rec(obj,cC)
            % ~(cC) signal reconstruction from frame coefficients C
            % should satisfy rec(dec)=Id
            
            n=2.^obj.ts.size_dyadic;
            if isvector(cC)
                cC=reshape(cC,n);
            end
            % 2d-transform, first along columns, then along rows:
             % normalize to get an orthogonal transform 
            yn=idct2(cC);
            yn=yn*sqrt(numel(yn));
            
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x using an ONB
           
            % 2d-transform
            yn= dct2(reshape(x,obj.ts.N));
            yn=yn/sqrt(numel(yn));% to get an orthogonal transform     
        end
        
        function xn= synthesize(obj,y)
            % yn=~(x) synthesize as inverse operation to analyze            
            xn= idct2(reshape(y,obj.ts.N));
            xn=xn*sqrt(numel(xn)); % to get an orthogonal transform
        end
        
        
    end
    
    %% queries
    
    methods
        
         function bn= basisname(obj)
            % bn=~() basis name
            bn='dct2';
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
        
        
        function cC=vec2C(obj,vC)
            % convert vector to form usable for reconstruction rec
            cC=reshape(vC,obj.ts.N);
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
        
        function [nodes,mask]=nodesPDF(obj,c,params)
            % nodes=~(c) nodes sampled with compression rate c
            obj.requireOnly(c>1,'local','compression rate in ]0,1[');
            obj.requireOnly(~obj.ts.isemptydata,'local','non-empty signal');
            if nargin<3 || ~isfield(params,'P')
                params.deg=5; % Variable density polymonial degree
            end
            if ~isfield(params,'iter')
                params.iter=10;
            end
            params.box=[0,1];  % origin is in upper left corner
            DN = obj.ts.N; % data Size
            p = 1/c;% [0.25];% undersampling factor                        
            % generates the sampling PDF concentrated at center
            pdf = genPDF(DN,p,params);
            % generate sampling mask:
            mask = genSampling(pdf,params);            
            % compute nodes:
            nodes=SampledNodes(find(mask),obj.ts.N,obj);
        end
        
        function [nodes,mask,res]=nodesOnRadialLines(obj,RL,c)
            % ~([RL,c]) nodes on RL radial lines yielding a compression rate c;   
            % if c is given, RL will be adapted (by iteration) to reach c;
            % nodes are used to subsample the transform with
            % emphasis on the region around index (1,1).
                
            function d= compression_deviation(RL1)
                mask = LineMask(round(RL1),min(n));
                nodes=find(mask); % transform vector of indices of non-zeros
                d=c-obj.ts.numel/numel(nodes);
            end
            
            if nargin <2 || isempty(RL)
                RL=100; % start value
            end
            n=obj.ts.size;
            nodes=[];
            res=struct;
            res.RL=RL;
            
            % find number RL of lines corresponding to compression c
            % using radial lines starting from center of image (faster
            % but yielding the same c).
            if nargin>2 && ~isempty(c)
                options = optimset('TolX',0.01);
                [res.RL,res.d,res.flag,res.output]=fzero(@compression_deviation,RL,options);
                res.RL=round(res.RL);
            end
            
            % compute RL radial lines in one quadrant:
            mask=zeros(n);
            angles=linspace(0,pi/2,res.RL);
            x1=1;
            y1=1;
            for j=1:res.RL
                x2=n(2);
                y2=tan(angles(j))*(n(2)-1)+1;
                if y2>n(1)
                    x2=(n(1)-1)/tan(angles(j))+1;
                    y2=n(1);
                end
                x2=max(1,min(n(2),x2));
                y2=max(1,min(n(1),y2));
                [x y]=bresenham(x1,y1,x2,y2);
                idx=sub2ind(size(mask),y,x);
                mask(idx)=1;
            end
            nodes=SampledNodes(find(mask),obj.ts.N,obj);
        end
        
        function cC=dec2graph(obj,cC)
            % cC=~(cC) post process decomposition obj.dec for graphical output
            cC=reshape(cC, obj.ts.N);
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
            w=DCosineTrafo2();
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
        end
        
        function graph_trafo(obj,open_new)
            %  show coefficients of transformed signal in open figure window
            obj.require(obj.isTrafodone,'local','transform result is set');
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo);
            end
            
            img=obj.C;
            if obj.use_repfun
                img= obj.repfun(img);
                img(img<-10)=-10;
            end
            imagesc(obj.dec2graph(img)) ; colorbar;
            cblabel(func2str(obj.repfun));
            
            tittext=['Transformation: ',obj.basisname, obj.add_signalname];
            title(tittext,'fontsize',12);
            
        end
        
        function graph_distribution(obj, open_new)
            % ~() show distribution of frame coefficients in open window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            img_dwt=obj.C2vec();
            HM=HMatrix(img_dwt);
            HM.graph_distribution(false);
            
            h=get(gca,'Title');
            present_titstr=get(h,'String');
            addtit=['distribution of coefficients: ',...
                obj.basisname];
            
            new_titstr=[addtit,present_titstr];
            title(new_titstr,'fontsize',12);
            
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



