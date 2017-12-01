classdef TV_IC_MRA3 < MultiResTransform3_Partial
    % Frame transform composition of any 3d partial multiresolution transform
    % in one reference frame combined with total variation of original signal along z axis;
    % This class can be used as a decoder of data type MultiResTransform3.
    %
    % Example
    % ========
    %{
    % --- create signal:
        signal3=Signal3D.make_fromVideo('tennis.avi',[1,3]);
        % signal3.play_signal;
        signal3.graph_totalVar;
    % --- create transform object:
        w= TV_IC_MRA3(signal3);
    % --- choose 2d multi-resolution analysis:
        w.set_transform(Wavelet2D_mlab());
    % OR:  w.set_transform(Curvelet2_clab());
    % --- perform decomposition (analysis):
        w.dec;
        w.show_trafo;
    % --- test inversion properties:
        ok=w.test_DecRec(); disp(ok);
        ok=w.test_AnalyzeSynthsize(); disp(ok);
    % --- test if adjointSynthesis is the adjoint of synthesize           
        [ok,err]= TV_IC_MRA3.test_adjointSynthesis();
        disp(['Test passed=',num2str(ok), ' with error ',num2str(err)]);
    %}
    
    properties
        motionCF  %@<MotionCompensation> motion compensation field
        use_motion %@<logical>
    end
    
    
    
    %% commands
    methods
        
        function obj= TV_IC_MRA3 (signal,hftrafo)
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
            obj.use_motion=false;
            
            obj.ensure(true,'local','check invariant');
        end
        
        function set_algo(obj)
            set_algo@MultiResTransform3(obj);
            if obj.use_motion
                obj.algo.name='(MRA_{xy}\circ \pi_1, \Delta_z), motion comp.';
            else                
                obj.algo.name='(MRA_{xy}\circ \pi_1,D_z)';
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
            w2=1; %w2=(1-obj.anweight1);
            w1=1; %w1=obj.anweight1;
            
            obj.ftrafo.ts.set_signal(w2*obj.ts.xn(:,:,obj.zref));
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
                if ~obj.use_motion
                    obj.C{j-1}=w1*...
                        (obj.ts.xn(:,:,j)-obj.ts.xn(:,:,j-1));
                else
                    obj.C{j-1}=w1*...
                        (obj.motionCF.moveBackward(j-1)-obj.ts.xn(:,:,j-1));
                end
            end
            for j=obj.zref+1:L
                if ~obj.use_motion
                    obj.C{j}=w1*...
                        (obj.ts.xn(:,:,j)-obj.ts.xn(:,:,j-1));
                else
                    obj.C{j}=w1*...
                        (obj.ts.xn(:,:,j)-obj.motionCF.moveForward(j-1));
                end
            end
            
            
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wc)
            % yn=~(wc) % signal reconstruction of decomposition wc,
            % which must be of the same datatype as result of obj.dec.
            
            L=obj.ts.size(3);
            yn=zeros(obj.ts.size);
            w2=1; %w2=(1-obj.anweight1);
            w1=1; %w1=obj.anweight1;
            
            % undo MRA in reference frame:
            yn(:,:,obj.zref)=obj.ftrafo.rec(wc{obj.zref})/w2;
            cCref=yn(:,:,obj.zref);
            
            % undo TV
            for j=obj.zref+1:L
                if ~obj.use_motion
                    yn(:,:,j)=wc{j}/w1+cCref;
                else
                    yn(:,:,j)=wc{j}/w1+obj.motionCF.moveForward(j-1,cCref);
                end
                cCref=yn(:,:,j);
            end
            
            cCref=yn(:,:,obj.zref);
            for j=obj.zref-1:-1:1
                if ~obj.use_motion
                    yn(:,:,j)=cCref-wc{j}/w1;
                else
                    yn(:,:,j)=obj.motionCF.moveBackward(j,cCref)-wc{j}/w1;
                end
                cCref=yn(:,:,j);
            end
            
        end
        
    end
    
    
    %% queries
    
    methods
        
        function nc=numelC(obj,cC)
            % nc=~([cC]) number of elements in decomposition
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            L=obj.ts.size(3);
            if nargin<2
                cC1=obj.C{obj.zref};
            else
                cC1=cC{obj.zref};
            end
            nc=obj.ftrafo.numelC(cC1)+(L-1)*obj.ts.size(1)*obj.ts.size(2);
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose vector x into vector yn
            % Psi=Psi_1 o D_t, where D is the extended difference op. with IC
            x=reshape(x,obj.N);
            w2=1; %w2=(1-obj.anweight1);
            w1=1; %w1=obj.anweight1;
            
            L=obj.ts.size(3);
            % MRA on reference frame:
            cC=reshape(obj.ftrafo.analyze(w2*x(:,:,obj.zref)),[],1);
            numel_signal=obj.ts.size(1)*obj.ts.size(2);
            numel_C=length(cC);
            yn=zeros(numel_C+(L-1)*numel_signal,1);
            
            % TV: operating on signal (difference op.)
            offset=0;
            for j=2:obj.zref
                if ~obj.use_motion
                    yn(offset+1:offset+numel_signal)=reshape(w1*...
                        (x(:,:,j)-x(:,:,j-1)),[],1);
                else
                    yn(offset+1:offset+numel_signal)=reshape(w1*...
                        (obj.motionCF.moveBackward(j-1,x(:,:,j))-x(:,:,j-1)),...
                        [],1);
                end
                offset=offset+numel_signal;
            end
            % assign MRA of reference frame:
            yn(offset+1:offset+numel_C)=cC;
            offset=offset+numel_C;
            
            for j=obj.zref+1:L
                if ~obj.use_motion
                    yn(offset+1:offset+numel_signal)=reshape(w1*...
                        (x(:,:,j)-x(:,:,j-1)),[],1);
                else
                    yn(offset+1:offset+numel_signal)=reshape(w1*...
                        (x(:,:,j)-obj.motionCF.moveForward(j-1,x(:,:,j-1))),...
                        [],1);
                end
                offset=offset+numel_signal;
            end
            
        end
        
        function xn= synthesize(obj,y)
            % xn=~(y) synthesize signal from decomposition y, where
            % y must belong to the same data type as result of analyze
            % apply 2d transform to each z-projection
            % Phi=I_t o Phi_1, where I_t is the integration op. (left
            % inverse of D_t).
            xn=zeros(obj.ts.size);
            L=obj.ts.size(3);
            w2=1; %w2=(1-obj.anweight1);
            w1=1; %w1=obj.anweight1;
            
            % undo MRA of reference frame:
            numel_signal=obj.ts.size(1)*obj.ts.size(2);
            numel_C=length(y)-(L-1)*numel_signal;
            frame_size=[obj.ts.size(1), obj.ts.size(2)];
            offset0=(obj.zref-1)*numel_signal;
            
            if ~obj.ftrafo.isTrafodone
                obj.ftrafo.dec;
            end
            
            xn(:,:,obj.zref)=reshape(...
                obj.ftrafo.synthesize(y(offset0+1:offset0+numel_C)/w2),...
                frame_size);
            
            % undo TV (i.e. integrate):
            cCref=xn(:,:,obj.zref);
            offset=offset0+numel_C;
            for j=obj.zref+1:L
                if ~obj.use_motion
                    xn(:,:,j)=cCref+reshape(y(offset+1:offset+numel_signal)/w1,frame_size);
                else
                    xn(:,:,j)=obj.motionCF.moveForward(j-1,cCref)+reshape(y(offset+1:offset+numel_signal)/w1,frame_size);
                end
                offset=offset+numel_signal;
                cCref=xn(:,:,j);
            end
            
            cCref=xn(:,:,obj.zref);
            offset=offset0;
            for j=obj.zref-1:-1:1
                if ~obj.use_motion
                    xn(:,:,j)=cCref-reshape(y(offset-numel_signal+1:offset)/w1,frame_size);
                else
                    xn(:,:,j)=obj.motionCF.moveBackward(j,cCref)-reshape(y(offset-numel_signal+1:offset)/w1,frame_size);
                end
                offset=offset-numel_signal;
                cCref=xn(:,:,j);
            end
            
        end
        
        function yn=adjointSynthesis(obj,x)
            % yn=~(x) apply adjoint of synthesis operator Phi to vector x
            if ~obj.use_motion
                yn=obj.adjointSynthesis_NoMotionComp(x);
            else
                yn=obj.adjointSynthesis_MotionComp(x);
            end
        end
        
        function yn=adjointSynthesis_NoMotionComp(obj,x)
            % xn=~(y) apply adjoint of synthesis operator Phi to vector y
            % must tbe redefined because of lack of orthogonality
            % i.e. adjoint(Phi) ~= Psi ;
            % this function returns vector xn
            % the following is no longer true: xn=obj.analyze(y);
            % \Phi=I_{t,\mu} \circ \Phi_1
            % => \Phi^*= \Psi_1 \circ I^*_{t,\mu}
            % TEST: test passed!
            % w= TV_IC_MRA3(signal3);
            % y= w.analyze(w.ts.xn);
            % w.use_motion=false;
            % err=dot(reshape(w.synthesize(y),[],1),y(:))-dot(reshape(w.adjointSynthesis(y),[],1),y(:));
            % disp(err);
            %
            
            x=reshape(x,obj.N);
            L=obj.ts.size(3);
            w2=1; %w2=(1-obj.anweight1);
            w1=1; %w1=obj.anweight1;
            
            numel_signal=obj.ts.size(1)*obj.ts.size(2);
            if ~obj.ftrafo.isTrafodone
                obj.dec;
            end
            numel_C=length(obj.ftrafo.C);
            yn=zeros(numel_C+(L-1)*numel_signal,1);
            
            % adjoint of integration: operating on signal
            % (while removing weight to avoid divergence):
            
            % 1) for j<zref: yn= -sum_{i=1}^n x_i
            offset=0;
            x0=zeros(numel_signal,1);
            for j=1:obj.zref-1
                x0=x0-reshape(x(:,:,j)*w1,[],1);
                yn(offset+1:offset+numel_signal)=x0;
                offset=offset+numel_signal;
            end
            offset_zref=offset;
            x0_zref=-x0; % i.e. +sum_{i=1}^{zref-1} x_i
            
            % 2) for j>zref: yn= sum_{i=n}^L x_i
            offset=numel_C+(L-2)*numel_signal;
            x0=zeros(numel_signal,1);
            for j=L:-1:obj.zref+1;
                x0=x0+reshape(x(:,:,j)*w1,[],1);
                yn(offset+1:offset+numel_signal)=x0;
                offset=offset-numel_signal;
            end
            % final x0 is sum_{i=zref+1}^L x_i
            
            % assign value for reference frame:
            % y_zref= sum_{i=1}^L x_i
            offset=offset_zref;
            x0_zref=x0_zref+x0+reshape(x(:,:,obj.zref)*w2,[],1);
            % MRA on reference frame % (while removing weight to avoid divergence):
            cC=reshape(obj.ftrafo.analyze(x0_zref),[],1);
            yn(offset+1:offset+numel_C)=cC;
            
        end
        
        function yn=adjointSynthesis_MotionComp(obj,x)
            % xn=~(y) apply adjoint of synthesis operator Phi to vector y
            % must tbe redefined because of lack of orthogonality
            % i.e. adjoint(Phi) ~= Psi ;
            % this function returns vector xn
            % the following is no longer true: xn=obj.analyze(y);
            % \Phi=I_{t,\mu} \circ \Phi_1
            % => \Phi^*= \Psi_1 \circ I^*_{t,\mu}
            % TEST: [ok, err]= TV_IC_MRA3.test_adjointSynthesis();
            %        disp(err);
            %
            
            x=reshape(x,obj.N);
            L=obj.ts.size(3);
            w2=1; %w2=(1-obj.anweight1);
            w1=1; %w1=obj.anweight1;
            z0=obj.zref;
            
            numel_signal=obj.ts.size(1)*obj.ts.size(2);
            if ~obj.ftrafo.isTrafodone
                obj.dec;
            end
            numel_C=length(obj.ftrafo.C);
            yn=zeros(numel_C+(L-1)*numel_signal,1);
            
            % adjoint of integration: operating on signal
            % (while removing weight to avoid divergence):
            
            % 1) for j<zref: yn= -sum_{i=1}^n x_i
            offset=0;
            if z0>1
                x0=-x(:,:,1)*w1;
                yn(offset+1:offset+numel_signal)=x0(:);
                offset=offset+numel_signal;
            else
                x0=0;
            end
            for j=2:z0-1
                x0=obj.motionCF.moveBackwardAdj(j-1,x0)-x(:,:,j)*w1;
                yn(offset+1:offset+numel_signal)=x0(:);
                offset=offset+numel_signal;
            end
            offset_zref=offset;
            x0_zref=-x0; % i.e. +sum_{i=1}^{zref-1} x_i
            
            % 2) for j>zref: yn= sum_{i=n}^L x_i
            offset=numel_C+(L-2)*numel_signal;
            if z0<L
                x0=x(:,:,L)*w1;
                yn(offset+1:offset+numel_signal)=x0(:);
                offset=offset-numel_signal;
            else
                x0=0;
            end
            for j=L-1:-1:z0+1;
                x0=obj.motionCF.moveForwardAdj(j,x0)+x(:,:,j)*w1;
                yn(offset+1:offset+numel_signal)=x0(:);
                offset=offset-numel_signal;
            end
            % final x0 is sum_{i=zref+1}^L x_i
            
            % assign value for reference frame:
            % y_zref= sum_{i=1}^L x_i
            offset=offset_zref;
            if z0>1
                x0_zref=obj.motionCF.moveBackwardAdj(z0-1,x0_zref);
            end
            x0_zref=x0_zref+x(:,:,z0)*w2;
            if z0<L
                x0_zref=x0_zref+obj.motionCF.moveForwardAdj(z0,x0);
            end
            % MRA on reference frame % (while removing weight to avoid divergence):
            cC=reshape(obj.ftrafo.analyze(x0_zref),[],1);
            yn(offset+1:offset+numel_C)=cC;
            
        end
        
        function coeff=detcoef(obj,lvl,type)
            % coeff=~(level[,type]) extract detail coefficients
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.ftrafo.C=obj.C{obj.zref};
            coeff=obj.ftrafo.detcoef(lvl,type);
        end
        
        function coeff= appcoef(obj)
            % coeff=~() extract approximate coefficients from wavelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.ftrafo.C=obj.C{obj.zref};
            coeff=obj.ftrafo.appcoef();
        end
        
    end
    
    %% filter
    
    methods
        
        function hardthresh1_inplace(obj,thresh)
            % filter all coefficients obj.C
            L=obj.ts.size(3);
            for j=1:L
                if j~=obj.zref
                    obj.C{j}(abs(obj.C{j})<thresh)=0;
                else
                    obj.ftrafo.C=obj.C{j};
                    obj.ftrafo.hardthresh1_inplace(thresh);
                    obj.C{j}=obj.ftrafo.C;
                end
            end
        end
        
        function hardthresh2_inplace(obj,thresh1, thresh2)
            % filter with 2 thresholds thresh1<thresh2
            L=obj.ts.size(3);
            for j=1:L
                if j~=obj.zref
                    obj.C{j}(abs(obj.C{j})<thresh2)=0;
                else
                    obj.ftrafo.C=obj.C{j};
                    obj.ftrafo.hardthresh2_inplace(thresh1,thresh2);
                    obj.C{j}=obj.ftrafo.C;
                    
                end
            end
        end
        
        function softthresh1_inplace(obj,thresh)
            % simple soft filter: reduce all coeff. obj.C by size
            L=obj.ts.size(3);
            for j=1:L
                if j~=obj.zref
                    obj.C{j}(abs(obj.C{j})<thresh)=0;
                    % now also reduce other coefficients
                    obj.C{j}=obj.C{j}-sign(obj.C{j})*thresh;
                else
                    obj.ftrafo.C=obj.C{j};
                    obj.ftrafo.softthresh1_inplace(thresh);
                    obj.C{j}=obj.ftrafo.C;
                end
            end
        end
        
        
    end
    
    %% transforms
    methods
        
        function Cvec=C2vec(obj,cC)
            % convert transformation result to vector form
            % using C2vec of obj.ftrafo.
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin<2
                cC=obj.C;
            end
            M=obj.ts.N;
            L=length(cC);
            
            % frame zref is transformed by ftrafo
            cC1=obj.ftrafo.C2vec(cC{obj.zref});
            numel_C=numel(cC1);
            numel_signal=M(1)*M(2);
            Cvec=zeros((L-1)*numel_signal+numel_C,1);
            offset=0;
            for j=1:L
                if j~=obj.zref  % these frames are result of difference op.
                    Cvec(offset+1:offset+numel_signal)=reshape(cC{j},[],1);
                    offset=offset+numel_signal;
                else  % this frame is result of ftrafo
                    Cvec(offset+1:offset+numel_C)=cC1;
                    offset=offset+numel_C;
                end
            end
        end
        
        function cC=vec2C(obj,vC)
            % cC=~(vC) convert frame coefficients to vector form
            L=length(obj.C);
            
            numel_signal=obj.ts.size(1)*obj.ts.size(2);
            numel_C=length(vC)-(L-1)*numel_signal;
            
            cC=cell(L,1);
            offset=0;
            for j=1:L
                if j~=obj.zref  % these frames are result of difference op.
                    cC{j}=vC(offset+1:offset+numel_signal);
                    offset=offset+numel_signal;
                else
                    cC{j}=obj.ftrafo.vec2C(vC(offset+1:offset+numel_C));
                    offset=offset+numel_C;
                end
            end
        end
        
    end
    
    
    %% statics
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w= TV_IC_MRA3 ();
        end
        
        function [ok,err]= test_adjointSynthesis(use_motion,mfield,s3)
            % ok=~([use_motion,mfield,s3]) tests if adjointSynthesis is the adjoint of synthesize
            % with or without motion compensation.
            % mfield ... (empty) motion field, e.g. mfield=MotionVF();
            % s3 ... Signal3D
            % e.g. 
            % [ok, err]=TV_IC_MRA3.test_adjointSynthesis(true,MotionQuotientField());
            %
            if ~exist('use_motion','var') || isempty(use_motion)
                use_motion=true;
            end
            if ~exist('mfield','var') || isempty(mfield)               
               mfield=MotionVF();    
            end
            if ~exist('s3','var') || isempty(mfield)               
                L=3; fn='tennis.avi'; 
                s3=Signal3D.make_fromVideo(fn,L);  
                s3.shiftSignalValues(-1); % shift awy from 0 intensities for quotient field
            end
            
            if ~isequal(mfield.size_space,s3.size(1:2))
                mfield=mfield.make_fromSignal3(s3);
            end
            
            w = TV_IC_MRA3(s3);            
            w.motionCF=mfield;
            w.use_motion=use_motion;
            y= w.analyze(w.ts.xn);
            err=abs(dot(reshape(w.synthesize(y),[],1),y(:))-dot(reshape(w.adjointSynthesis(y),[],1),y(:)));
            ok=err<1e-9;
        end
        
    end
    
    %% graphics
    
    methods
        
        
        function [mat, levelshown]=C2graph(obj, lvl)
            % mat=~() convert coefficient vector to format suitable for
            % graphical output
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || lvl<=obj.level,'local', 'smaller than max. #levels');
            
            if nargin < 2 || isempty(lvl)
                lvl=obj.ftrafo.level;
            end
            levelshown=lvl;
            L=obj.ts.size(3);
            img=obj.ftrafo.C2graph(lvl);
            mat=zeros(size(img,1),size(img,2),L);
            obj.ftrafo.C=obj.C{obj.zref};
            mat(:,:,obj.zref)= obj.ftrafo.C2graph(lvl);
            maxval=max(reshape(mat(:,:,obj.zref),[],1));
            for j=1:L
                if j~=obj.zref
                    img=obj.C{j};
                    img(abs(img)>1e-3)=maxval;
                    mat(:,:,j)= img;
                end
            end
            
        end
        
        function show_trafo(obj)
            % show transform C  in open figure window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            suptitle(obj.modified_title(),14);
            
            lvl=obj.level;
            L=obj.ts.size(3);
            sd=factor_subplots(L);
            maxAbsC=zeros(L,1);
            
            for j=1:L
                subplot(sd(1),sd(2),j);
                if j==obj.zref
                    obj.ftrafo.C=obj.C{j};
                    M1=obj.ftrafo.C2graph(lvl);
                    imagesc(obj.repfun(M1));
                    colormap('default');
                    colorbar;
                    if ~isempty(obj.repfun)
                        cblabel_str=func2str(obj.repfun);
                        cblabel(cblabel_str);
                    end
                else
                    M1=reshape(obj.C{j},obj.ts.size(1),obj.ts.size(2));
                    imagesc(obj.repfun(M1));
                    colormap('default');
                    colorbar;
                end
                maxAbsC(j)=max(abs(M1(:)));
                
                
                maxstr=num2str(maxAbsC(j),'%3.1e');
                if j==obj.zref
                    title(['z=',num2str(j),', MRA_{xy}=',obj.ftrafo.algo.name,...
                        ' (',obj.ftrafo.basisname,'), max=',maxstr],...
                        'fontsize',12);
                else
                    title(['z=',num2str(j),', ',obj.algo.name,', max=',...
                        maxstr],'fontsize',12);
                end
            end
            
        end
        
    end
    
    %% invariant
    
    
end


