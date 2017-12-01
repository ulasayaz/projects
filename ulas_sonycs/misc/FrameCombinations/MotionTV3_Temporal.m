classdef MotionTV3_Temporal<TV3_Temporal
    % total variation along a motion vector field of a 3d signal
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
    %     f= MotionTV3_Temporal(signal);
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
        motion  %@<MotionVF> motion vector field
    end
    
    
    %% constructor and commands
    methods
        
        function obj=MotionTV3_Temporal(signal)
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
            ss=signal.N;
            obj.motion= MotionVF(zeros([ss(1),ss(2),min(0,ss(3)-1),2]));
            
        end % constructor
        
        function set_algo(obj)
            %~() tags of algorithm
            obj.algo.version=0.9;
            obj.algo.versiondate='28.4.2014';
            obj.algo.name='motion total variation';
            obj.algo.toolbox=[];
            obj.algo.mfilename=class(obj);
            obj.algo.author='G. Troll';
        end
        
        function set_motion(obj,vf,M)
            % ~(vf,[]) set motion vector field
            obj.requireOnly(length(size(vf))==4 || (nargin >2 && ismatrix(vf)),...
                'local','vf dimensions ok');
            obj.requireOnly(ismember(obj.ts.numel,[vf.size(1),vf.size(1)*vf.size(2)]),...
                'local', 'match with signal');
            obj.requireOnly(ismember(obj.ts.size(3)+1,[vf.size(2),vf.size(min(length(size(vf)),3))]),...
                'local', 'match with signal');
            if nargin <3
                M=[];
            end
            obj.motion.set_motionVF(obj,vf,M);
        end
        
        function v=delta(obj,j)
            % v=~(j) difference frame(j+1)-frame(j) using motion 
            v=obj.ts.xn(:,:,j);
            w=obj.ts.xn(:,:,j-1);
            idxvals=~isnan(obj.motion(:,j));
            M=obj.ts.size(1);
            N=obj.ts.size(2);
            % subtract signal values of corresponding pixels:
            v=reshape(v((1:M*N)-obj.motion(idxvals,j))-w(idxvals),M, N);
        end
        
        function dec(obj)
            % ~() decomposition whose inverse is rec, i.e. rec(dec)=Id
            % yielding an ONB (check with test_framefeatures() );
            % i.e. operation computeFrame by applying the linear dec operation
            % on the canonical basis (e_i) will yield an ONB (but is no
            % longer involution).
            obj.require(obj.ts.numel>0,'local','non-empty sample size');
            
            % 3d-transform
            hN=obj.ts.size;
            obj.C=zeros(hN(1),hN(2),hN(3)-1);
            for j=2:hN(3)
                obj.C(:,:,j)= obj.delta(j);
            end
            
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
            bn='motion TV';
        end
        
        
        
        %% transforms
        
        
        
    end
    
    methods (Hidden)
        
        
        
    end
    
    methods
        %% coordinate transformations
        
        
        
    end
    
    
    %% tests
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=MotionTV3_Temporal();
        end
        
    end
    
    %% graphics
    methods
        
        
        
        
    end
    
    %% invariant
    
    
end




