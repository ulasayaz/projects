classdef TV3_Temporal<FrameTrafo
    % total variation of a 3d signal along 3rd (temporal) dimension
    % This class can be used e.g. as one of the frame transforms of FrameTrafoPair.
    %
    % GT Stuttgart, 2014
    %
    % Example
    % =========
    % 1.) define a signal, e.g.
    %     signal2=Signal2D.make_fromImage('cameraman.bmp');
    %     L=32; signal=Signal3D.make_CyclefromLowerDimSig(signal2,1,L);
    %
    % 2.) create object
    %     f= TV3_Temporal(signal);
    %
    % 3.) apply transform
    %     f.dec;
    % 4.) show result
    %     f.graph_trafo;
    %
    % Test frame features:
    %     ok=f.test_DecRec; disp(ok);
    %     res=f.test_framefeatures(); display(res);
    %
    %
    
    properties (Access=public)
       
    end
    
    properties (SetAccess=protected)
       
    end
    
    
    %% constructor and commands
    methods
        
        function obj=TV3_Temporal(signal)
            % constructor
            obj.requireOnly(nargin==0 || isempty(signal) || isa(signal,'Signal3D'),...
                'local', ' signal data type ok');
            if nargin==0
                signal=Signal3D();
            end
            obj.ts=signal;
                       
            obj.set_algo;
            obj.fig=1;
            obj.use_repfun=true;
            obj.fontsize=12;                       
            
        end % constructor
        
        function set_algo(obj)
            %~() tags of algorithm
            obj.algo.version=0.9;
            obj.algo.versiondate='28.4.2014';
            obj.algo.name='temporal TV';   
            obj.algo.toolbox=[];
            obj.algo.mfilename=class(obj);
            obj.algo.author='G. Troll';
        end               
        
        function reset_Trafo(obj)
            % reset fourier trafo (e.g. after signal has changed)
            obj.C=[];
        end
        
        function tv=totalVar(obj)
            % tv=~()  total variation along z-axis
            obj.require(obj.ts.numel>0,'local','non-empty sample size');
            if ~obj.isTrafodone
                obj.dec;
            end
            tv=Signal2D(sum(abs(obj.C),3));            
        end
        
        function dec(obj)
            % ~() decomposition whose inverse is rec, i.e. rec(dec)=Id
            % yielding an ONB (check with test_framefeatures() );
            % i.e. operation computeFrame by applying the linear dec operation
            % on the canonical basis (e_i) will yield an ONB (but is no
            % longer involution).
            obj.require(obj.ts.numel>0,'local','non-empty sample size');
            
            % 3d-transform
            obj.C= diff(obj.ts.xn,1,3);            
            
        end
        
        function yn=rec(obj,cC)
            % ~(cC) signal reconstruction from frame coefficients C
            % should satisfy rec(dec)=Id
            
            % 3d-transform
            L=obj.ts.size(3);
            yn=zeros(obj.ts.size);
            % reconstrict first plane (j=1) from initial values of signal
            yn(:,:,1)=obj.ts.xn(:,:,1);
            for j=2:L
                yn(:,:,j)=yn(:,:,j-1)+cC(:,:,j-1);
            end            
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose vector x using an ONB
            
            yn= diff(reshape(x,obj.ts.N),1,3);           
        end
        
        function xn= synthesize(obj,y)
            % yn=~(x) synthesize as inverse operation to analyze
            L=obj.ts.size(3);
            M=obj.ts.size;
            xn=zeros(obj.ts.size);
            y=reshape(y,M(1),M(2),L-1);
            % first plane j=1  not reconstructible without initial values
            for j=2:L
                xn(:,:,j)=xn(:,:,j-1)+y(:,:,j-1);
            end       
            
        end
        
        
    end
    
    %% queries
    
    methods
        
        function bn= basisname(obj)
            % bn=~() basis name
            bn='TV';
        end    
        
        function ok=isTrafodone(obj)
            % ok=~(): is result of a Fourier trafoi available?
            ok=~isempty(obj.C);
        end
        
        function nv=frame_length(obj)
            % n=~() number of elements of transform
            nv=obj.ts.size(1)*obj.ts.size(2)*(obj.ts.size(3)-1);
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
            w=TV3_Temporal();
        end
        
    end
    
    %% graphics
    methods
                
                
        function graph_trafo(obj, open_new)
            % ~([open_new, z) show projection to z-plane of transform C 
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin <3 || obj.isnatural(z),...
                'local','z is permissible z coordinate');
            if ~exist('open_new', 'var') || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            tv=obj.totalVar();
            tv.signalname=[obj.algo.name,obj.add_signalname]; 
            tv.graph_signal(false);
            cblabel(func2str(obj.repfun));           
            
        end
        
        function graph_projection(obj, open_new,z)                         
            % ~([open_new, z) show projection to z-plane of transform C 
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin <3 || obj.isnatural(z),...
                'local','z is permissible z coordinate');
            if ~exist('z', 'var') || isempty(z)
                z=1;
            end
            if ~exist('open_new', 'var') || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            imagesc(obj.repfun(obj.C(:,:,z))); 
            colorbar;
            cblabel(func2str(obj.repfun));
            tittext=[obj.algo.name,', z=',num2str(z),obj.add_signalname];
            title(tittext,'fontsize',12);
            
        end
        
        function show_trafo(obj)
            % ~() show signal reconstructed from transform in new window
            prepfigure(obj.fig,obj.algo,obj.figopt);
            L=obj.ts.size(3);
            
            sd=factor_subplots(L-1);
            suptitle([obj.algo.name,obj.add_signalname],14);
            for j=1:L-1
                subplot(sd(1),sd(2),j);
                obj.graph_projection(false,j);
                title(['z=',num2str(j)],'fontsize',12);
            end
        end
        
         
         
        
       
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is Signal3D';
            ok= isa(obj.ts,'Signal3D');
        end
    end
    
end



