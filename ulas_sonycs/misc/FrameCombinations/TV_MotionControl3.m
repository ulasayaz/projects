classdef TV_MotionControl3 < TV_IC_MRA3
    % combine TV_IC_MRA3 with a motion control
    % where TV_IC_MRA3 was a frame transform composition of any 3d partial multiresolution transform
    % in one reference frame combined with total variation of original signal along z axis;
    % This class can be used as a decoder of data type MultiResTransform3.
    %    
    % Example
    % ========
    % --- create signal:
    %     signal3=Signal3D.make_fromVideo('tennis.avi',[1,3]);
    %     % signal3.play_signal;
    %     signal3.graph_totalVar;
    % --- create transform object:
    %     w= TV_MotionControl3(signal3);
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
    
    properties (SetAccess=protected)
        motionVF  %@<matrix> size (M*N,L-1) motion vector field
    end
    
    
    
    %% commands
    methods
        
        function obj= TV_MotionControl3 (signal,hftrafo)
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
            obj = obj@TV_IC_MRA3(signal,hftrafo);
            
            obj.ensure(true,'local','check invariant');
        end
        
        function set_algo(obj)
            set_algo@MultiResTransform3(obj);
            obj.algo.name='(MRA_{xy}\circ \pi_1,Motion_z)';
        end
        
        function set_motionVF(obj,vf)
            % ~(vf2d) set the motionVF operating on vectors using a 2d VF
            obj.requireOnly(isnumeric(vf),'local','vf is matrix');
            ss1=size(vf);
            ss=obj.ts.size;
            obj.requireOnly((ss1(1)==ss(1)*ss(2) && ss1(2)==ss(3)-1) || ...
                (isequal(ss1(1:3),ss-[0,0,1]) && ss1(4)==2),'local',...
                '1d or 2d vector field operating on frames');
            if ss1(1)==ss(1)*ss(2)
                obj.motionVF=vf;
            else % convert 2d to 1d:
                obj.motionVF=obj.motion_sub2ind(vf);
            end                                    
            
        end
        
        function dec(obj,lev)
            % partial discrete multiresolution analysis
            % apply transform to each z-projection;
            % the 2d MRA gets a weight of L to yield
            % an approx. isometry (which is important fo the convergence
            % of the CS-solvers).
            obj.requireOnly(obj.ts.numel>0,'local','non-empty sample size');
            
            % instead of applying (1-anweight1) to result of ftrafo (which might be
            % cell), apply it to signal of obj.ftrafo and use linearity:
            obj.ftrafo.ts.set_signal((1-obj.anweight1)*obj.ts.xn(:,:,obj.zref));
            if nargin <2
                lev=obj.ftrafo.deepestlev();
            end
            
            L=obj.ts.size(3);
            obj.C=cell(L,1);
            
            % MRA applied only to reference plane
            obj.ftrafo.dec(lev);
            obj.C{obj.zref}=obj.ftrafo.C;
            
            % TV: operating on signal planes:
            for j=2:obj.zref
                obj.C{j-1}=obj.anweight1*obj.delta(j);
            end
            for j=obj.zref+1:L
                obj.C{j}=obj.anweight1*obj.delta(j);
            end
            
            
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wc)
            % yn=~(cC) % signal reconstruction of decomposition wc,
            % which must be of the same datatype as result of obj.dec.
            
            L=obj.ts.size(3);
            yn=zeros(obj.ts.size);
            
            % undo MRA in reference frame:
            yn(:,:,obj.zref)=obj.ftrafo.rec(wc{obj.zref})/(1-obj.anweight1);
            cCref=yn(:,:,obj.zref);
            
            % undo TV
            for j=obj.zref+1:L
                yn(:,:,j)=obj.add(wc{j}/obj.anweight1,cCref);
                cCref=yn(:,:,j);
            end
            
            cCref=yn(:,:,obj.zref);
            for j=obj.zref-1:-1:1
                yn(:,:,j)=obj.add(cCref,-wc{j}/obj.anweight1);
                cCref=yn(:,:,j);
            end
            
        end
        
    end
    
    
    %% queries
    
    methods
                
        
        function v=delta(obj,j)
            % difference between consecutive frame using the motion vector
            % field.
            v=obj.ts.xn(:,:,j);
            w=obj.ts.xn(:,:,j-1);
            idxvals=~isnan(obj.motionVF(:,j));
            M=obj.ts.size(1);
            N=obj.ts.size(2);
            % subtract signal values of corresponding pixels:
            v=reshape(v((1:M*N)-obj.motionVF(idxvals,j))-w(idxvals),M, N);
        end
        
        function v=add(obj,a,b,j)
            % integrate delta values
            v=a;
            idxvals=~isnan(obj.motionVF(:,j));
            v=v(idxvals)+b(idxvals);
            M=obj.ts.size(1);
            N=obj.ts.size(2);
            % undo shift introduced by motion vector field
            v=reshape(v((1:M*N)+obj.motionVF(idvals,j)),M,N);
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose vector x into vector yn
            x=reshape(x,obj.N);
            
            L=obj.ts.size(3);
            % MRA on refrence frame:
            cC=reshape(obj.ftrafo.analyze((1-obj.anweight1)*x(:,:,obj.zref)),[],1);
            numel_signal=obj.ts.size(1)*obj.ts.size(2);
            numel_C=length(cC);
            yn=zeros(numel_C+(L-1)*numel_signal,1);
            
            % TV: operating on signal
            offset=0;
            for j=2:obj.zref
                yn(offset+1:offset+numel_signal)=reshape(obj.anweight1*...
                    (obj.ts.xn(:,:,j)-obj.ts.xn(:,:,j-1)),[],1);
                offset=offset+numel_signal;
            end
            % assign MRA of reference frame:
            yn(offset+1:offset+numel_C)=cC;
            offset=offset+numel_C;
            
            for j=obj.zref+1:L
                yn(offset+1:offset+numel_signal)=reshape(obj.anweight1*...
                    (obj.ts.xn(:,:,j)-obj.ts.xn(:,:,j-1)),[],1);
                offset=offset+numel_signal;
            end
            
        end
        
        function xn= synthesize(obj,y)
            % xn=~(y) synthesize signal from decomposition y, where
            % y must beong to the same data type as result of analyze
            % apply 2d transform to each z-projection
            xn=zeros(obj.ts.size);
            L=obj.ts.size(3);
            
            % undo MRA of reference frame:
            numel_signal=obj.ts.size(1)*obj.ts.size(2);
            numel_C=length(y)-(L-1)*numel_signal;
            frame_size=[obj.ts.size(1), obj.ts.size(2)];
            offset0=(obj.zref-1)*numel_signal;
            
            xn(:,:,obj.zref)=reshape(...
                obj.ftrafo.synthesize(y(offset0+1:offset0+numel_C)/(1-obj.anweight1)),...
                frame_size);
            
            % undo TV:
            cCref=xn(:,:,obj.zref);
            offset=offset0+numel_C;
            for j=obj.zref+1:L
                xn(:,:,j)=cCref+reshape(y(offset+1:offset+numel_signal)/obj.anweight1,frame_size);
                offset=offset+numel_signal;
                cCref=xn(:,:,j);
            end
            
            cCref=xn(:,:,obj.zref);
            offset=offset0;
            for j=obj.zref-1:-1:1
                xn(:,:,j)=cCref-reshape(y(offset-numel_signal+1:offset)/obj.anweight1,frame_size);
                offset=offset-numel_signal;
                cCref=xn(:,:,j);
            end
            
        end
        
    end
    
    %% test
    
    methods 
        
        function vf2d= simulate_motionVF(obj)
            % vf2d=~() simulates 2d motion vector field
            
            ss=obj.ts.size;
            rectangle=ceil(0.25*[ss(1), ss(2)]);
            vf2d=zeros(ss(1),ss(2),ss(3)-1,2);
            offset=0;
            for j=1:ss(3)-1;
                vf2d(offset+1:offset+rectangle(1),j,j,1)=NaN;
                vf2d(j,offset+1:rectangle(2)+offset,j,1)=NaN;               
                vf2d(j+1:rectangle(1)+j,j+1:rectangle(2)+j,j,1)=ones(rectangle);
                
                vf2d(offset+1:offset+rectangle(1),j,j,2)=NaN;
                vf2d(j,offset+1:offset+rectangle(2),j,2)=NaN;
                vf2d(j+1:rectangle(1)+j,j+1:rectangle(2)+j,j,2)=ones(rectangle);
                offset=offset+1;
            end
        end
        
        function s3=apply_motionVF(obj, signal2d)
            % s3=~(signal2d) apply the motion vector field to a 2d signal
            obj.requireOnly(isa(signal2d,'Signal2D') && ~iempty(obj.motionVF)...
                && size(obj.motionVF,1)==signal2d.numel,...
                'local','size matches with motion vectro field');
            M=signal2d.size(1);
            N=signal2d.size(2);
            L=size(obj.motionVF,2)+1;
            s3=Signal3D(zeros(M,N,L));
            s3(:,:,1)=signal2d;
            
            for j=1:L-1
                s3.xn(:,:,j+1)=obj.add(s3.xn(:,:,j),obj.motionVF(:,j));
            end            
            
        end
        
        
        
        
    end
    
    %% transforms
    methods
        
        function vf1d=motion_sub2ind(obj,vf2d)
            % v=~(sub) convert 2d motion to 1d motion;
            % is inverse of motion_ind2sub;
            obj.require(isequal(size(vf2d), [obj.ts.size-[0,0,1],2]),...
                'local','argument is a 2 dim. vector field operating on L-1 frames');
            ss=obj.ts.size;
            M=ss(1);
            N=ss(2);
            L=ss(3);
            
            [plane2,plane1]=meshgrid(1:M,1:N);
            idx1=(1:M*N)';
            vf1d=zeros(M*N,L-1);
            for j=1:L-1
                shifted1=plane1+vf2d(:,:,j,1);
                shifted2=plane2+vf2d(:,:,j,2);
                vf1d(:,j)=sub2ind([M,N],shifted1(:),shifted2(:))-idx1;
            end
            
        end
        
        function vf2d=motion_ind2sub(obj,vf1d)
            % v=~(sub) convert 2d motion to 1d motion;
            % is inverse of motion_sub2ind;
            obj.require(isequal(size(vf1d), [obj.ts.size(1)*obj.ts.size(2),obj.ts.size(3)-1]),...
                'local','argument is a 1 dim. vector field operating on L-1 frames');
            ss=obj.ts.size;
            M=ss(1);
            N=ss(2);
            L=ss(3);
           
            idx1=(1:M*N)';
            vf2d=zeros(M,N,L-1,2);
            [plane2,plane1]=meshgrid(1:M,1:N);
            
            for j=1:L-1
                [P1,P2]=ind2sub([M,N],vf1d(:,j)+idx1);
                vf2d(:,:,j,1)=reshape(P1,M,N)-plane1;
                vf2d(:,:,j,2)=reshape(P2,M,N)-plane2;
            end
            
        end
        
    end
    
    
    %% statics
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w= TV_MotionControl3 ();
        end
        
    end
    
    %% graphics
    
    methods
        
        
    end
    
    %% invariant
    
    
end



