classdef MRA_o_TV3 < MultiResTransform3_Partial
    % Frame transform composition of any 3d partial multiresolution transform
    % in xy planes combined with total variation along z axis;
    % The 2d MRA is applied to each z-plane separately;
    % This class can be used as a decoder of data type MultiResTransform3.
    %
    % Example
    % ========
    % --- create signal:
    %     signal3=Signal3D.make_fromVideo('tennis.avi',[1,3]);
    %     % signal3.play_signal;
    %     signal3.graph_totalVar;
    % --- create transform object:
    %     w= MRA_o_TV3(signal3);   
    % --- choose 2d multi-resolution analysis:
    %     w.set_transform(Wavelet2D_mlab());
    % OR:  w.set_transform(Curvelet2_clab());
    % --- perform decomposition (analysis):
    %     w.dec;
    %     w.show_trafo;
    % --- test inversion properties:
    %     ok=w.test_DecRec(); disp(ok);
    %     ok=w.test_AnalyzeSynthsize(); disp(ok);
    %
    
    properties
        
    end
    
    
    %% commands
    methods
        
        function obj= MRA_o_TV3(signal,hftrafo)
            % constructor
            if nargin==0
                signal=Signal3D();
            end
            if nargin<2 || isempty(hftrafo)
                if license('test','Wavelet_Toolbox')
                    hftrafo=Wavelet2D_mlab();
                else
                    hftrafo=Wavelet2D_wlab();
                end
            end
            obj = obj@MultiResTransform3_Partial(signal,hftrafo);            
           
            obj.ensure(true,'local','check invariant');
        end
        
        function set_algo(obj)
            set_algo@MultiResTransform3(obj);
            obj.algo.name='(MRA_{xy} \circ (D_z,I)';
        end

        function dec(obj,lev)
            % partial discrete multiresolution analysis
            % apply transform to each z-projection;
            % the 2d MRA gets a weight of L to yield
            % an approx. isometry (which is important fo the convergence
            % of the CS-solvers).
            obj.requireOnly(obj.ts.numel>0,'local','non-empty sample size');
            
            obj.ftrafo.ts.set_signal(obj.ts.xn(:,:,1));
            if nargin <2
                lev=obj.ftrafo.deepestlev();
            end
            
            L=obj.ts.size(3);
            obj.ftrafo.ts.set_signal(obj.ts.xn(:,:,1));
            obj.ftrafo.dec(lev);
            
            
            obj.C=cell(L,1);
            
            
             % TV: operating on signal planes:
            for j=2:obj.zref
                obj.C{j-1}=obj.anweight1*...
                    (obj.ts.xn(:,:,j)-obj.ts.xn(:,:,j-1));
            end
            for j=obj.zref+1:L
                obj.C{j}=obj.anweight1*...
                    (obj.ts.xn(:,:,j)-obj.ts.xn(:,:,j-1));
            end
            obj.C{obj.zref}=(1-obj.anweight1)*obj.ts.xn(:,:,obj.zref);
                                    
            % MRA applied to each obj-C-cell, which is a matrix!            
            for j=1:L
                obj.ftrafo.ts.set_signal(obj.C{j});
                obj.ftrafo.dec(lev);
                obj.C{j}=obj.ftrafo.C;  
            end
           
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wc)
            % yn=~(cC) % signal reconstruction of decomposition wc,
            % which must be of the same datatype as result of obj.dec.
            
            L=obj.ts.size(3);
            yn=zeros(obj.ts.size);
            h=yn;
            
            % undo MRA:
            for j=1:L
                h(:,:,j) = obj.ftrafo.rec(wc{j});
            end
            
            cCref0=h(:,:,obj.zref)/(1-obj.anweight1);
            cCref=cCref0;
            
            % undo TV
            for j=obj.zref+1:L
                yn(:,:,j)=h(:,:,j)/obj.anweight1+cCref;                
                cCref=yn(:,:,j);
            end
            
            cCref=cCref0;
            for j=obj.zref-1:-1:1
                yn(:,:,j)=cCref-h(:,:,j)/obj.anweight1;                
                cCref=yn(:,:,j);
            end
            yn(:,:,obj.zref)=cCref0;
                        
        end
        
    end
    
    
    %% queries
    
    methods
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose vector x into vector yn
            % apply 2d transform to each z-projection
            x=reshape(x,obj.N);
            
            L=obj.ts.size(3);                        
            
            obj.C=cell(L,1);
            cC=cell(L,1);
            
            % TV: operating on signal
            for j=2:obj.zref
                cC{j-1}=obj.anweight1*...
                    (x(:,:,j)-x(:,:,j-1));
            end
            for j=obj.zref+1:L
                cC{j}=obj.anweight1*...
                    (x(:,:,j)-x(:,:,j-1));
            end
            cC{obj.zref}=(1-obj.anweight1)*x(:,:,obj.zref);
                                    
            % MRA
            yn=zeros(obj.ftrafo.numelC,L);            
            for j=1:L
                yn(:,j)=reshape(obj.ftrafo.analyze(cC{j}),[],1);                 
            end            
        end
        
        function xn= synthesize(obj,y)
            % xn=~(y) synthesize signal from decomposition y, where
            % y must beong to the same data type as result of analyze
            % apply 2d transform to each z-projection
            xn=zeros(obj.ts.size);
            L=obj.ts.size(3);
            y=reshape(y,[],L);
                        
            h=xn;            
            % undo MRA:
            for j=1:L
                h(:,:,j) = obj.ftrafo.synthesize(y(:,j));
            end
            
            cCref0=h(:,:,obj.zref)/(1-obj.anweight1);
            cCref=cCref0;
            
            % undo TV
            for j=obj.zref+1:L
                xn(:,:,j)=h(:,:,j)/obj.anweight1+cCref;                
                cCref=xn(:,:,j);
            end
            
            cCref=cCref0;
            for j=obj.zref-1:-1:1
                xn(:,:,j)=cCref-h(:,:,j)/obj.anweight1;                
                cCref=xn(:,:,j);
            end
            xn(:,:,obj.zref)=cCref0;                       
        end
        
    end
    
    
    %% transforms
    methods
        
        
    end
    
    
    %% statics
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w= MRA_o_TV3();
        end
        
    end
    
    %% graphics
    
    methods
        
    end
    
    %% invariant
    
    
end







