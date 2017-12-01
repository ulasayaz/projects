classdef Fourier2 < FourierTrafo
    % 2-dim. Fourier transform with graphical output;
    %
    % Example:
    % Example
    % 1.) define a signal, e.g.
    %     signal=Signal2D.make_zernike(48);
    %     signal=Signal2D.make_fromImage('cameraman.bmp');
    %
    % 2.) create object
    %     f=Fourier2(signal);
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
        
        function obj=Fourier2(signal)
            % constructor
            assert(nargin==0 || isempty(signal) || isa(signal,'Signal2D'),...
                 ' signal is 2D');
            if nargin==0
                signal=Signal2D();
            end
            obj = obj@FourierTrafo(signal);
            
        end % constructor
       
    end
    
    %% queries
    
    methods
        
        function bn= basisname(obj)
            % bn=~() basis name
            bn='fft2';
        end    
        
        %% transforms
        
        function [nodes,mask,res]=nodesOnRadialLines(obj,RL,c)
            % ~([RL,c]) nodes on RL radial lines yielding a compression rate c;   
            % if c is given, RL will be adapted (by iteration) to reach c;
            % nodes are used to subsample the transform with
            % emphasis on low frequencies
            
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
            if nargin>2 && ~isempty(c)
                options = optimset('TolX',0.01);
                [res.RL,res.d,res.flag,res.output]=fzero(@compression_deviation,RL,options);
                res.RL=round(res.RL);
            end
            
            compression_deviation(res.RL);
            nodes=SampledNodes(nodes,obj.ts.N,obj);
            
        end
        
        
    end
    
    
    %% tests
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Fourier2();
        end
        
    end
    
    %% graphics
    methods
         
        function graph_trafo(obj,open_new)
            % show |fft2| without additional data in open figure window
            % used together with obj.dec
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            imagesc(obj.dec2graph(abs(obj.repfun(obj.C)))); colorbar;
            cblabel(func2str(obj.repfun));
            tittext=['|',obj.basisname,'|',obj.add_signalname];
            title(tittext,'fontsize',12);
        end
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is Signal2D';
            ok= isa(obj.ts,'Signal2D');
        end
    end
    
end




