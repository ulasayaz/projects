classdef Fourier1 < FourierTrafo
    % 1-dim. Fourier transformations with graphical output;
    %
    % Example:
    % Example
    % 1.) define a signal, e.g.
    %     signal=TimeSignal.make_sinus(3);    
    % 2.) create object
    %     f=Fourier1(signal);
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
    
    properties 
        
    end
    
    
    %% constructor and commands
    methods
        
        function obj=Fourier1(signal)
            % constructor
            assert(nargin==0 || isempty(signal) || isa(signal,'TimeSignal'),...
                ' signal is 1D');
            if nargin==0
                signal=TimeSignal();
            end
            
            obj = obj@FourierTrafo(signal);
            
        end % constructor
        
    end
    
    %% queries
    
    methods
        
        function bn= basisname(obj)
            % bn=~() basis name
            bn='fft';
        end                
       
        %% transforms
        
        function [nodes,mask]=nodesExpRnd(obj,c)
            % nodes=~(c) compute random sample of nodes which is
            % exponentially distributed with highest density at origin of
            % centered Fourier coordinates.
            obj.require(~obj.ts.isemptydata,'local',' signal needs data');
            obj.requireOnly(c>=1,'local','compression rate >1');
            
            % define mask and nodes in centered Fourier coordinates:
            L=floor(obj.ts.length/2); 
            mu=L/5; 
            cnodes=round(L/c);
            nodes=unique(max(1,min(L,round(exprnd(mu,[1,cnodes]+1)))));
            mask=zeros(L,1);
            mask(nodes)=1;
            mask=[flipud(mask);mask];
            
            % shift back to non-centered Fourier coordinates
            mask=fftshift(mask);            
            nodes=SampledNodes(find(mask),obj.ts.N,obj);
            
            obj.ensureOnly(max(nodes.data)<=numel(obj.ts.xn) && min(nodes.data)>=1,'local',...
                'nodes within range');            
        end
        
        
        
        
    end
    
   
    
    %% tests
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Fourier1();
        end
        
    end
    
    %% graphics
    methods
            
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




