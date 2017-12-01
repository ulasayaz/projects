classdef Fourier3 < FourierTrafo
    % 3-dim. Fourier transform with graphical output;
    %
    % Example:
    %{
    %1.) define a signal, e.g.
        signal2=Signal2D.make_fromImage('cameraman.bmp');
        L=32; signal=Signal_Factory.make_CyclefromSignal2D(signal2,1,L);
    
    %2.) create object
        f=Fourier3(signal);
    
    %3.) apply transform
        f.dec;
    %4.) show result
        f.graph_trafo;
    
    %Test frame features:
    res=f.test_framefeatures(); display(res);
    %}    
    %GT Stuttgart, 2014
    
    
    properties 
      
    end
   
   
    %% constructor and commands
    methods
        
        function obj=Fourier3(signal)
            % constructor
            assert(nargin==0 || isempty(signal) || isa(signal,'Signal3D'),...
                 ' signal is 3D');
            if nargin==0
                signal=Signal3D();
            end
            obj = obj@FourierTrafo(signal);
            
        end % constructor
       
    end
    
    %% queries
    
    methods
        
        function bn= basisname(obj)
            % bn=~() basis name
            bn='fft3';
        end    
        
        %% transforms
        
       function [nodes,mask]=nodesOnRadialLines(obj,L)
            % ~([L]) nodes are distributed on L radial lines
            % nodes are used to subsample the transform with
            % emphasis on the region around (0,0) in the
            % fftshift coordinate system.
           
            if nargin <2 || isempty(L)
                L=150;
            end
            n=obj.ts.size;
            mask=zeros(n);
           
            % center of lines
            x0=ceil(n(1)/2);
            y0=ceil(n(2)/2);
            z0=ceil(n(3)/2);
            
            Ls=max(1,ceil(3*n/L));
            
            % define a planar grid on 3 planes of the cuboid:
            [Y1,X1,Z1]=meshgrid(1:Ls(1):n(1),1:Ls(2):n(2),n(3));
            [Y2,X2,Z2]=meshgrid(1:Ls(1):n(1),n(2),1:Ls(3):n(3));
            [Y3,X3,Z3]=meshgrid(n(1),1:Ls(2):n(2),1:Ls(3):n(3));
            
            X2=squeeze(X2); Y2=squeeze(Y2); Z2=squeeze(Z2);
            X3=squeeze(X3); Y3=squeeze(Y3); Z3=squeeze(Z3);
            
            % connect points on planar grids with center and beyond:
            for m=1:size(X1,1)
                for n=1:size(X1,2)
                    x1e=X1(m,n); y1e=Y1(m,n); z1e=Z1(m,n);
                    x1s=max(1,2*x0-x1e); y1s=max(1,2*y0-y1e); z1s=max(1,2*z0-z1e);
                    [x,y,z]=bresenham3([x1s,y1s,z1s],[x1e,y1e,z1e]);
                    idx=sub2ind(size(mask),x,y,z);
                    mask(idx)=1;
                end
            end
            for m=1:size(X2,1)
                for n=1:size(X2,2)        
                    x1e=X2(m,n); y1e=Y2(m,n); z1e=Z2(m,n);
                    x1s=max(1,2*x0-x1e); y1s=max(1,2*y0-y1e); z1s=max(1,2*z0-z1e);
                    [x,y,z]=bresenham3([x1s,y1s,z1s],[x1e,y1e,z1e]);
                    idx=sub2ind(size(mask),x,y,z);
                    mask(idx)=1;
                end
            end
            for m=1:size(X3,1)
                for n=1:size(X3,2)
                    x1e=X3(m,n); y1e=Y3(m,n); z1e=Z3(m,n);
                    x1s=max(1,2*x0-x1e); y1s=max(1,2*y0-y1e); z1s=max(1,2*z0-z1e);
                    [x,y,z]=bresenham3([x1s,y1s,z1s],[x1e,y1e,z1e]);
                    idx=sub2ind(size(mask),x,y,z);
                    mask(idx)=1;
                end
            end
                
            % shift to non-centered Fourier space
            mask=fftshift(mask);            
            nodes=SampledNodes(find(mask),obj.ts.N,obj);
        end
        
        
    end
    
    
    %% tests
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Fourier3();
        end
        
    end
    
    %% graphics
    methods
         
        function graph_trafo(obj, open_new, z)
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
            
            imagesc(obj.dec2graph(abs(obj.repfun(obj.C(:,:,z))))); 
            colorbar;
            cblabel(func2str(obj.repfun));
            tittext=['|',obj.basisname,'|, z=',num2str(z),obj.add_signalname];
            title(tittext,'fontsize',12);
            
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






