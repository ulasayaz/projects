classdef MotionVF <MotionField
    % motion vector field operating on frames (matrices)
    % Example
    % ========
    %{
    %--- create vector field
        %-- larger matrix:
        params.M=6; params.N=7;
        %-- or smaller default matrix:
        % params=struct;
        % params.compute_weightField=true;
        vf=MotionVF.test_SimpleIllustration(params);
        vf.set_method_compensate(1);
        vf.show_AllAssociatedFields;
    
        %-- using external method to compensate motion:
        vf.set_method_compensate(0);
        vf.show_AllAssociatedFields;
    
    %--- test adjointness of adjoint moves with random vector fields:
        ok=MotionVF.test_adjointMoves(params);
        disp(['result of adjoint test: ',num2str(ok)]);
    
    %--- motion simulations:
        countRec=1;
        %countRec=10;
        params=struct;  params.sig=0; M=128; N=M; L=3; %noise: params.sig=0.01;
        % optional: params.shiftsize=10; params.ulc=round([M/2,N/2]);
        sf=Signal_Factory();
        [s3,motion]=sf.make_MovingRectangles(countRec,M,N,L,params);
        s3.fig=3;
        %-- show frames in subplots
        s3.show_signal;
        %-- visualize motion by computing pixelwise total variation:
        s3.graph_totalVar;
        %-- visualize motion by playing hte video:
        s3.play_signal;
    
        motion.fig=3;
        %--- show vector field in 1 plot:
        motion.graph_field;
        %--- show effect of motion compensation:
        motion.show_occlusions();
    % motion compensation setting occluded points to 0
        motion.show_motionCompensation();
    % motion compensation using the quotient field for occluded points:
        s3.shiftSignalValues(-1);
        motion.set_weightFieldFromSignal3(s3);
        motion.graph_weightField(1);
        motion.show_motionCompensation();
    
        motion.test_motionCompensation;
    
    %-- motion from videos:
        L=3; %fn='tennis.avi';
        fn='riksch1.avi';
        s3=Signal3D.make_fromVideo(fn,L);
        s3.shiftSignalValues(-1);
        opts=struct;
        %-- choose method (OF is best):
        opts.method_estim='OF'; % method_estim='FSBM';
        vf=MotionVF.make_fromSignal3(s3,opts);
        vf.fig=5;
        vf.show_NormsAnglesDiff(1);
        vf.show_normWithAngle;
        vf.show_occlusions();
        vf.show_crossings();
    % -- treat occlusions by keeping moving pixels
        vf.show_motionCompensation();
        vf.test_motionCompensation(1);
        vf.test_motionCompensation();
    
    % -- treat occlusions by using the quotient field:
        vf.set_weightFieldFromSignal3(vf.signal3);
        vf.show_motionCompensation();
        vf.test_motionCompensation(1);
        vf.test_motionCompensation();
            
    %}
    
    properties (SetAccess=protected)
        weightField       %@<MotionQuotientField>
    end
    
    properties (SetAccess=protected, Hidden)
        adjvf,adjvfPast   %@<MotionVF>  adjoint vector field is accessed by obj.adjField
        trafotype         %@<integer> (0 is original) is vf result of some basic transform?
        domainVstar       %@<MotionQuotientField>
    end
    
    properties
        zref   % reference frame
        scale  %@<double> scaling factor for quiver plot (default 1 automatic, 0 no rescaling)
        LineWidth  %@<int> linewidth of quiver plot
    end
    
    %% commands
    methods
        
        function obj= MotionVF (vf)
            % constructor
            assert(nargin <1 || (length(size(vf))==4 && size(vf,4)==2),'vf dimensions ok');
            if nargin==0
                vf=zeros(1,1,1,2);
            end
            obj = obj@MotionField();
            
            obj.set_field(vf);
            obj.scale=0;
            obj.motionfname='vector field V';
            obj.LineWidth=1;
            obj.meth_comp=1;  % internal method is 1;
            obj.trafotype=0;  % 0 is original
            % fix pointer to weightField so that it can be copied to adjField:
            
            
            obj.ensure(true,'local','check invariant');
        end
        
       
        function set_field(obj,mvf,opts)
            % ~(matvf,opts]) set the vector field to round(mvf);
            % if opts is given and opts.zref is not empty
            % a special rounding scheme based on compensated summation will 
            % be used (see below).
            % drawback: tests have shown that the number of occlusions
            % rises when using compensated summation.
            %
            obj.requireOnly(isempty(mvf) || length(size(mvf))==4 && size(mvf,4)==2,...
                'local','matrix dimensions ok');            
            if nargin <3 || ~isstruct(opts), opts=struct; end
            obj.requireOnly(~isfield(opts,'zref')|| isempty(opts.zref)...
                || (opts.zref>=1 && opts.zref<=size(mvf,3)+1),...
                'local','reference plane admitted');
            
            if ~isfield(opts,'zref')
                opts.zref=[];  
                % temporal size of mvf is one less than temporal size of
                % signal:
                %zref=max(1,ceil((size(mvf,3)+1)/2));  % mid point z-plane
            end
            z0=opts.zref;
            obj.zref=z0;
            L=size(mvf,3);              
            
            % round already here to integers to avoid rounding the weightField:
            if isempty(opts.zref) || isempty(mvf) || L==1
                obj.mat=round(mvf);            
            else % L>1, zref==L+1 possible!                
                % special rounding scheme based on compensated summation;
                % i.e. the cumulative sum of the
                % rounded vector fields stays as close as possible to
                % the cumulative sum of the unrounded vector field:
                              
                obj.mat=zeros(size(mvf));
                mvf1=[];
                if z0>1
                    mvf1=flipdim(mvf(:,:,1:z0-1,:),3);
                end
                mvf2=[];
                if z0<L+1                
                    mvf2=mvf(:,:,z0:end,:);
                end
                
                LH2=size(mvf2,3);
                LH1=size(mvf1,3);
                                
                if LH2<=1
                    mvf2=round(mvf2);
                else
                    mvf2=diff(round(cumsum(mvf2,3)),1,3);                    
                end
                if LH1<=1
                    mvf1=round(mvf1);
                else
                    mvf1=flipdim(diff(round(cumsum(mvf1,3)),1,3),3);                    
                end
                if z0>2
                    obj.mat(:,:,1:z0-2,:)=mvf1;
                end
                if z0<L
                    obj.mat(:,:,z0+1:L,:)=mvf2;
                end
                if z0>1
                    obj.mat(:,:,z0-1,:)=round(mvf(:,:,z0-1,:));
                end
                if z0<L+1
                    obj.mat(:,:,z0,:)=round(mvf(:,:,z0,:));
                end
                % test:
                % ok1=z0>L || all(abs(reshape(cumsum(mvf(:,:,z0:end,:),3)-cumsum(obj.mat(:,:,z0:end,:),3),[],1))<=0.5);
                % ok2=z0<2 || all(abs(reshape(cumsum(flipdim(mvf(:,:,1:z0-1,:),3),3)-cumsum(flipdim(obj.mat(:,:,1:z0-1,:),3),3),[],1))<=0.5);
                % ok3= all(abs(reshape(mvf(:,:,min(z0,L),:)-obj.mat(:,:,min(z0,L),:),[],1))<=0.5);
            end
            obj.reset_adjField();
            if isempty(obj.weightField)
                obj.weightField=MotionQuotientField();
            end
            
            
        end
        
        function set_signal(obj,s3)
            % ~(s3) set 3d signal
            obj.requireOnly(isa(s3,'Signal3D'),'local','3d signal');
            obj.signal3=s3;
            obj.weightField.set_signal(s3);
        end
        
        function set_adjvf(obj,vf)
            obj.adjvf=vf;
        end        
        
        function set_adjvfPast(obj,vf)
            obj.adjvfPast=vf;
        end  
        
        function set_method_compensate(obj,n)
            % ~(n) sets active method used for motion compensation
            % retrieve active method with obj.method_compensate
            obj.requireOnly(obj.isnatural0(n),'local','n is nono-neagative integer');
            obj.meth_comp=n;
        end
        
        function set_weightField(obj,wf)
            % ~(wf) set weight field
            obj.requireOnly(isa(wf,'MotionQuotientField'),...
                'local','data type MotionQuotientField');
            obj.weightField=wf;
        end
        
        function set_weightFieldFromSignal3(obj,s3)
            % ~(s3) compute and set weightField from 3d signal s3
            obj.requireOnly(isa(s3,'Signal3D'),'local','s3 is 3d signal');
            obj.requireOnly(s3.min*s3.max>1e-6,'local',...
                'signal values must be bounded away from 0. Use s3.shiftSignalValues.');
            obj.weightField=MotionQuotientField.make_fromSignal3(s3);
            obj.adjField.set_weightField(obj.weightField);
            obj.weightField.set_domain(obj.occlusions);
        end
        
        function set_weightField2scalar(obj,val)
            % ~(value) set all weightField values to scalar value
            % val=0 sets all occluded points to 0
            % val=1 sets all occluded points to their value before moving
            obj.requireOnly(isscalar(val),'local','needs scalar value');
            obj.weightField.set_field(val*ones([obj.size_space,obj.size_time]));
            obj.adjField.set_weightField(obj.weightField);
            obj.weightField.set_domain(obj.occlusions);
        end
        
        
    end
    
    
    %% transformations
    methods
        
        function vf=inverseField(obj)
            % vf=~() reflected (inverse) vector field for move backward
            vf=obj.clone;
            vf.mat=-obj.mat;
            vf.set_weightField(obj.weightField.inverseField);
            vf.weightField.set_domain([]);
            vf.adjvf=[];
            vf.motionfname='inverse field -V';
            obj.trafotype=2;
        end
        
        function vf=adjField(obj)
            % vf=~() gets or computes adjoint field of present field,
            % i.e.vectors pointing backwards from target to source;
            % result will only be well-defined if each
            % target can be reached by at most one non-zero vector.
            % points where no vector of the original field ends get NaN.
            % adjoining is not yet an involution (only after the call
            % of obj.set_motionVF2Biadjoint).
                        
            if ~isempty(obj.adjvf)  % need not be computed
                vf=obj.adjvf;
            else
                % must be computed
                vf=obj.computeAdjoint(obj); % compute adj of obj
                % after the next line, adjoining isn't yet involution
                % (only after  obj.set_motionVF2Biadjoint has been called).
                obj.adjvf=vf;   % save result
            end
        end
        
        function vf=adjFieldOfInverse(obj)
            % vf=~() computes adjoint (-V)^* of inverse field -V;
            % (this is different from the inverse of the adjoint -(V^*)
            %
            obj.requireOnly(~isempty(obj.mat),'local','needs forward field');
            
            if ~isempty(obj.adjvfPast)  % need not be computed
                vf=obj.adjvfPast;
            else
                % must be computed
                vi=obj.inverseField;
                occ=vi.occlusions;
                vf=obj.computeAdjoint(vi);
                % requires new domain of weightField:
                vf.weightField.set_domain(occ);
                % after the next line, adjoining isn't yet involution
                % (only after  obj.set_motionVF2Biadjoint has been called).
                vf.motionfname='adjoint of inverse (-V)^*';
                obj.adjvfPast=vf;   % save result
                
            end
        end
        
        function reset_motionField(obj)
            % ~() reset motion fields
            obj.mat=[];
            obj.weightField.reset_motionField();
            obj.reset_adjField();
        end
        
        function reset_adjField(obj)
            % ~() set adjoint fields to empty value (allows recalculation)
            obj.adjvf=[];
            obj.adjvfPast=[];
        end
        
    end
    
    %% queries
    methods
        
        function mx=xvals(obj)
            % mx=~() x-values of vector field
            mx=obj.mat(:,:,:,2);
        end
        
        function my=yvals(obj)
            % my=~() y-values of vector field
            my=obj.mat(:,:,:,1);
        end
                
        
        function ok=sizes_match(obj)
            % ok=~() check whether sizes of motion field, weight field and adjoint field match
            ok=isempty(obj.weightField) || (isequal(obj.size_space, obj.weightField.size_space) && ...
                isequal(obj.size_time,obj.weightField.size_time));
            ok=ok && ( isempty(obj.adjvf) || ((isequal(obj.size_space, obj.adjvf.size_space) && ...
                isequal(obj.size_time,obj.adjvf.size_time)) ) );
        end
        
        function occ=occlusions(obj,j)
            % occ=~([j]) occlusions of frame j, i.e. where the adjoint field is undef.
            obj.require(nargin <2 || isempty(j) || (j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if nargin <2
                j=[];
            end
            N=obj.size_space;
            L=obj.size_time;
            occ=false(N(1),N(2),L);
            vfa=obj.adjField;
            if isempty(j)
                for k=1:L
                    occ(:,:,k)=isnan(vfa.mat(:,:,k,2)); % nan in x val
                end
            else
                occ=isnan(vfa.mat(:,:,j,2));
            end
        end
        
        function cp=crossings(obj,j)
            % occ=~([j]) crossings of frame j (adjoint field ambiguous)
            obj.require(nargin <2 || isempty(j) || (j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if nargin <2
                j=[];
            end
            N=obj.size_space;
            L=obj.size_time;
            cp=false(N(1),N(2),L);
            vfa=obj.adjField;
            if isempty(j)
                for k=1:L
                    cp(:,:,k)=isnan(vfa.mat(:,:,k,2)) & ... % nan in x val
                        ~isnan(vfa.mat(:,:,k,1));           % number in y
                end
            else
                cp=isnan(vfa.mat(:,:,j,2))& ~isnan(vfa.mat(:,:,j,1));
            end
        end
        
        function vf=computeAdjoint(obj,vforig)
            % vf=~(vforig computes adjoint field of vforig, i.e. vectors backwards
            % from target to source;
            % calls mex function;
            % result will only be well-defined if each
            % target can be reached by at most one non-zero vector.
            % points where no vector of the original field ends get NaN.
            
            vf=vforig.clone;
            vf.reset_adjField;
            L=vforig.size_time;
            
            for j=1:L
                Dy=vforig.mat(:,:,j,1);
                Dx=vforig.mat(:,:,j,2);
                
                if obj.use_mex
                    % mex function called:
                    [vf.mat(:,:,j,2),vf.mat(:,:,j,1)]=AdjointMotion(Dx,Dy);
                else
                    [vf.mat(:,:,j,2),vf.mat(:,:,j,1)]=obj.AdjointMotion_local(Dx,Dy);
                end
                
            end
            vf.motionfname='adjoint field V^*';
            vf.trafotype=2;
        end
        
        function str= get_trafotype(obj)
            % str=~() transform applied to vf:
            % 0 is original, 1 is adjoint, 2 is inverse, 3 is biadjoint
            switch obj.trafotype
                case 0, str='original';
                case 1, str='adjoint';
                case 2, str='inverse';
                case 3, str='biadjoint';
                otherwise, str= 'unknown';
            end
        end
        
        function str=weightField_type(obj)
            % str=~() type of weightField
            U=unique(reshape(obj.weightField.mat,[],1));
            if length(U)==1
                str=['Q=',num2str(U)];
            else
                str='Q=auto';
            end
        end
        
        function str= method_compensate(obj)
            % str=~() active method used for motion compensation
            if obj.meth_comp==1
                str='pull on V^*';
            elseif obj.meth_comp==2
                str='push V';
            else
                str= 'reconstruct on V^*';
            end
            str=[str,', ',obj.weightField_type];
        end
        
        function vf2=clone(obj)
            % clone object
            obj.require(true,'local','check invariant');
            vf2=clone@MotionField(obj);
            vf2.set_adjvf(obj.adjvf);
            vf2.set_adjvfPast(obj.adjvfPast);
            vf2.set_weightField(obj.weightField);
            vf2.meth_comp=obj.meth_comp;
            vf2.scale=obj.scale;
            vf2.LineWidth=obj.LineWidth;
            obj.ensureOnly(~eq(obj,vf2),'local','new object created');
        end
        
        function vf2=clone_deep(obj)
            % clone object
            obj.require(true,'local','check invariant');
            vf2=obj.clone;
            if ~isempty(obj.signal3)
                vf2.signal3=obj.signal3.clone;
            end
            if ~isempty(obj.adjvf)
                vf2.set_adjvf(obj.adjvf.clone);
            end
            if ~isempty(obj.adjvfPast)
                vf2.set_adjvfPast(obj.adjvfPast.clone);
            end
            if ~isempty(obj.weightField)
                vf2.set_weightField(obj.weightField.clone);
            end     
            obj.ensureOnly(~eq(obj,vf2),'local','new object created');
        end
                
        function vf2=make_like(obj)
            % vf2=~() same as clone but without matrix content
            vf2=make_like@MotionField(obj);
            vf2.adjvf =[];
            vf2.adjvfPast=[];
            vf2.weightField=[];
            vf2.meth_comp=obj.meth_comp;
            vf2.scale=obj.scale;
            vf2.LineWidth=obj.LineWidth;
        end     
        
        function s=frame(obj,j)
            % s=~(j) retrieves frame number j from mat
            obj.requireOnly(min(j)>=1 && max(j)<=obj.size_time,'local',...
                'j admissible'); 
            
            st=obj.size_time;
            s=obj.clone_deep;  
            
            if ~isempty(obj.mat)
                s.mat=obj.mat(:,:,j,:);
            end            
            if ~isempty(obj.signal3)
                j2=[j,j(end)+1]; % signal has one element more
                s.signal3.xn=obj.signal3.xn(:,:,j2);
            end            
            if ~isempty(obj.adjvf)
                s.adjField.set_field(obj.adjField.mat(:,:,j,:));
            end
            if ~isempty(obj.adjvfPast)
                s.adjFieldPast.set_field(obj.adjFieldPast.mat(:,:,j,:));
            end
            if ~isempty(obj.weightField)  % must be last after adjField !
                s.set_weightField(obj.weightField.frame(j));
                s.adjField.set_weightField(s.weightField);
            end
            
            obj.ensureOnly(~isequal(obj,s),'local','new object created');
            obj.ensureOnly(obj.size_time==st && s.size_time==length(j),'local',...
                'frame length ok');            
            obj.ensureOnly(isempty(s.signal3) || s.signal3.size(3)==length(j)+1,...
                'local','new signal length ok');
            obj.ensureOnly(obj.sizes_match,'local','sizes of obj ok');
            obj.ensureOnly(s.sizes_match,'local','sizes of s ok');            
            obj.ensureOnly(isempty(s.adjvf) || isequal(s.adjvf.weightField,s.weightField),...
                'local', 'same weightField');            
%             obj.ensureOnly({s,j},'isempty(obj.mat) || isequaln(local{1}.adjvf.mat,obj.adjvf.mat(:,:,local{2},:))',...
%                              'vf copied');
%             obj.ensureOnly({s,j},'isempty(obj.adjvf) || isequaln(local{1}.adjvf.mat,obj.adjvf.mat(:,:,local{2},:))',...
%                             'adjvf copied');
%             obj.ensureOnly({s,j},'isempty(obj.weightField) || isequaln(local{1}.weightField.mat,obj.weightField.mat(:,:,local{2}))',...
%                             'weightField copied');
            
        end
        
        
        function y=moveForwardAdj(obj,j,x)
            % y=(j,x) adjoint of move forwards== move forward of V^*
            % test adjointness with test_adjointMoves.
            if nargin <3
                x=obj.signal3.xn(:,:,j);
            end
            y=obj.adjField.moveForward(j,x);
        end
        
        function y=moveBackwardAdj(obj,j,x)
            % y=(j,x) adjoint of move backwards== move forward of (-V)^*
            % (one needs the adjoint of the inverse field (-V)^*;
            % cannot use moveBackward, because (-V)^*
            % is different from the inverse of the adjoint -(V^*);
            % test adjointness with test_adjointMoves.
            if nargin <3
                % for backward move use signal index 1 ahead of motion field
                % index:
                x=obj.signal3.xn(:,:,j+1);
            end
            y=obj.adjFieldOfInverse.moveForward(j,x);
        end
        
        function y=moveForward(obj,j,x)
            % s2=~(j) compensate motion of j-th frame to get a
            % postdiction for the (j-1)-th frame.
            obj.requireOnly(j>=1 && j<=obj.size_time+1,'local',' j permitted');
            obj.requireOnly(nargin>=3 || isa(obj.signal3,'Signal3D'),'local','move argument available');
            % use only information from time step j:
            if nargin <3
                x=obj.signal3.xn(:,:,j);
            end
            
            if obj.meth_comp==1
                % internal method: pull operating on adjoint field V^*
                vf=obj.adjField();
                Ax = vf.mat(:,:,j,2);
                Ay = vf.mat(:,:,j,1);
                Q=vf.weightField.FieldOnDomain(j);
                % mex function called:
                y=PullMotion(x,Ax,Ay,Q);
                % occluded pixels possibly marked by NaN in Dxdbquit
                
                % e.g. for adjoint field adjvf
                % y(isnan(Ax))=0;
            elseif obj.meth_comp==2
                % internal method: push on original vector field V
                Ax = round(obj.mat(:,:,j,2));
                Ay = obj.mat(:,:,j,1);
                % mex function called:
                y=PushMotion(x,Ax,Ay);
            else
                % the external method <<reconstruct>>
                % operating on the adjoint vector field V^*:
                vf=obj.adjField();
                Dy=vf.mat(:,:,j,1);
                Dx=vf.mat(:,:,j,2);
                % calling external method:
                y=reconstruct(x,Dx,Dy);
            end
            
            % test values: w should be closer to xn(j-1) than v:
            % test=[obj.signal3.xn(j-1,:)',x(:),y(:),motionV'];
        end
        
        function y=moveBackward(obj,j,x)
            % s2=~(j) compensate motion of j-th frame to get a
            % prediction for the (j+1)-th frame.
            obj.requireOnly(j>=1 && j<=obj.size_time,'local',' j permitted');
            obj.requireOnly(nargin>=3 || isa(obj.signal3,'Signal3D'),'local','move argument available');
            % use only information from time step j:
            if nargin <3
                % for backward move use signal index 1 ahead of motion field
                % index:
                x=obj.signal3.xn(:,:,j+1);
            end
            
            if obj.meth_comp==1
                % internal method: pull on (-V)^*
                vf=obj.adjFieldOfInverse; % -> (-V)^*
                Ax = vf.mat(:,:,j,2);
                Ay = vf.mat(:,:,j,1);
                Q=vf.weightField.FieldOnDomain(j);
                % mex function called:
                y=PullMotion(x,Ax,Ay,Q);
                % occluded pixels possibly marked by NaN in Dx
                % e.g. for adjoint field adjvf
                %y(isnan(Ax))=0;
            elseif obj.meth_comp==2
                % internal method: push operating on -V
                Ax = -obj.mat(:,:,j,2);
                Ay = -obj.mat(:,:,j,1);
                % mex function called:
                y=PushMotion(x,Ax,Ay);
            else
                % the external method <<reconstruct>> operating on the adjoint vector field:
                vf=obj.adjFieldOfInverse; % -> (-V)^*
                Dy=vf.mat(:,:,j,1);
                Dx=vf.mat(:,:,j,2);
                y=reconstruct(x,Dx,Dy);
            end
            
            % test values: w should be closer to xn(j+1) than v:
            % test=[obj.signal3.xn(j+1,:)',x(:),y(:),motionV'];
        end
        
        
        function nf = norm(obj,p)
            % val=~([p]) l2-norm of vector field
            if nargin <2
                p=2;
            end
            ss=obj.size;
            hM=ss(1);
            hN=ss(2);
            hL=ss(3);
            
            vf2d=obj.mat;
            % sum over the vector components (2nd dim):
            nf=reshape(sum(abs(vf2d).^p,4).^(1/p),hM,hN,hL);
        end
        
        function av = angle(obj)
            % val=~() angles of vector field, is NaN for norms<0.5
            ss=obj.size;
            
            av=atan2(obj.yvals,obj.xvals);
            nv=obj.norm;
            cutoff=0.5;
            av(nv<cutoff)=NaN;
        end
        
        function val=meannz(obj, cutoff)
            % val=~() mean length of non-zero vectors
            if nargin <2
                cutoff=1e-9;
            end
            nf=obj.norm(2);
            nf=nf(nf>cutoff);
            val=mean(nf);
        end
        
        function d=diffVF(obj, vf2)
            % d=~(vf2) difference of 2 vector fields obj.mat-vf2.mat
            obj.requireOnly(isa(vf2,'MotionVF') && isequal(obj.size,vf2.size),...
                'local', '2 vector fields of the same size');
            d=obj.clone;
            d.set_field(obj.mat-vf2.mat);
            d.weightField.set_field((obj.weightField.mat./vf2.weightField.mat));
            d.motionfname=['difference ',obj.motionfname,' - ',vf2.motionfname];
        end
        
        function [s3,z0]=eulerIntegration(obj,z0)
            % s3= ~([z0]) reconstruct the video using the euler method;
            % starting from an initial condition on the reference
            % frame z0 the other frames are reconstructed by using
            % the euler method applied to the motion vector field.
            obj.requireOnly(~isempty(obj.signal3),'local','signal needed for IC');
            obj.requireOnly(nargin <2 || isempty(z0) || obj.isnatural(z0),'local',...
                ' z0 is index');
            if nargin <2 || isempty(z0)
                z0=max(1,ceil((obj.size(3)+1)/2));  % mid point z-plane
            end
            L=obj.signal3.size(3);
            z0=max(1,min(z0,L));
            yn=zeros(obj.signal3.size);
            yn(:,:,z0)=obj.signal3.xn(:,:,z0);
            
            for j=z0+1:L
               yn(:,:,j)=obj.moveForward(j-1,yn(:,:,j-1));                 
            end
            for j=z0-1:-1:1
               yn(:,:,j)=obj.moveBackward(j,yn(:,:,j+1));                 
            end
            
            s3=Signal3D(yn);
            s3.signalname=['Euler integration: ',obj.motionfname,', ', obj.signal3.signalname];
            s3.colormap_active=obj.signal3.colormap_active;
            s3.colorlimits=[s3.frame(z0).min-eps,s3.frame(z0).max+eps];
            
        end
        
    end
    
    %% test
    
    methods
        
        
    end
    
    methods (Static)
        
        function vf=test_SimpleIllustration(params)
            % vf=~([params]) test used for illustration purposes
            % if params is given it needs fields A,vy,vx.
            % representing teh first video frame A (possibly a second B) and a
            % vector field with x-component vx and y-component vy.
            if ~exist('params','var') || isempty(params)
                params=struct;
            end
            if ~isfield(params,'implement_crossings')
                params.implement_crossings=false;
            end
            if ~isfield(params,'isrand')
                params.isrand=false;
            end
            if ~isfield(params,'compute_weightField')
                params.compute_weightField=false;
            end
            if ~isfield(params,'weightField')
                params.weightField=[];
            end
            if ~isfield(params,'A')
                params.A=[0.01,   0.2, 0.4;
                    0.6, 0.8, 1];
                if params.isrand
                    c=numel(params.A);
                    params.A= reshape(params.A(randperm(c)),size(params.A));
                end
            end
            sizeA=size(params.A);
            if ~isfield(params,'vy') || ~isequal(sizeA,size(params.vy))
                if isequal(sizeA,[2,3]) && ~params.isrand
                    if ~params.implement_crossings
                        % y values:
                        params.vy=[ 0, 0,  0;
                            -1, 0, -1];
                        % x values:
                        params.vx=[0,  1,  0 ;
                            1,  0, -2 ];
                    else
                        % y values:
                        params.vy=[ 0, 0,  0;
                            -1, -1, -1];
                        % x values:
                        params.vx=[0,   1,  0 ;
                            1,  -1, -1 ];
                    end
                else
                    params.vy=randi(5,sizeA) -3;
                    params.vx=randi(5,sizeA) -3;
                    if params.implement_crossings
                        c=floor(sizeA/2);
                        params.vy(c(1),c(2))=-1;
                        params.vx(c(1),c(2))=1;
                        params.vy(c(1),c(2)+1)=-1;
                        params.vx(c(1),c(2)+1)=0;
                    end
                end
            end
            if ~isfield(params,'M')
                params.M=size(params.A,1);
            end
            if ~isfield(params,'N')
                params.N=size(params.A,2);
            end
            sv=[params.M,params.N];
            Mdiff=params.M-size(params.A,1);
            Ndiff=params.N-size(params.A,2);
            if  Mdiff>0 || Ndiff>0
                params.vy(2,2)=1;
                mp=ceil(Mdiff/2);
                np=ceil(Ndiff/2);
                params.A=wkeep(padarray(params.A,[mp,np]),sv);
                params.vy=wkeep(padarray(params.vy,[mp,np]),sv);
                params.vx=wkeep(padarray(params.vx,[mp,np]),sv);
            end
            
            assert(isequal(size(params.A),size(params.vy)) && ...
                isequal(size(params.vy),size(params.vx)),'same sizes');
            
            [M,N]=size(params.A);
            L=2;
            
            vf=MotionVF();
            vf.signal3=Signal3D();
            vf.signal3.signalname='sim';
            vf.color='k';
            vf.signal3.xn=zeros(M,N,L);
            vfmat=zeros(M,N,L-1,2);
            vfmat(:,:,1,1)=params.vy;
            vfmat(:,:,1,2)=params.vx;
            
            vf.signal3.xn(:,:,1)=params.A;
            vf.weightField.set_motionQF(zeros(M,N,L-1));
            vf.set_field(vfmat);
            vf.scale=0;
            vf.LineWidth=2;
            if isfield(params,'B') && isequal(size(params.A),params.B)
                vf.signal3.xn(:,:,2)=params.B;
            else
                vf.signal3.xn(:,:,2)=vf.moveForward(1);
                try
                    % set the occluded points to random numbers:
                    occ=vf.occlusions();
                    M=vf.signal3.xn(:,:,2);
                    idx=occ(:,:,1);
                    M(idx)=rand(sum(idx(:)),1);
                    vf.signal3.xn(:,:,2)=M;
                catch
                end
            end
            if ~isempty(params.weightField)
                vf.weightField.set_motionQF(params.weightField);
            else
                vf.weightField.set_motionQF(zeros(vf.signal3.size-[0,0,1]));
            end
            if isempty(params.weightField) || params.compute_weightField
                if vf.signal3.min*vf.signal3.max<1e-6
                    vf.signal3.shiftSignalValues(-1);
                end
                vf.set_weightFieldFromSignal3(vf.signal3);
            end
            if isempty(vf.weightField.domain)
                vf.weightField.set_domain(vf.occlusions);
            end
            assert(vf.invariant,'result satisfies invariant');
        end
        
        function [Ax,Ay]=AdjointMotion_local(Dx,Dy)
            % comput adjoint vector field using a matlab function instead
            % of C++.
            % points where no vectors end get nan:
            % hence initialize output matrices to nan:
            [M,N]=size(Dx);
            Ax=nan*zeros(M,N); % default domain is empty, i.e. all nan
            Ay=Ax;
            
            % loop over matrix
            for i=0:M-1
                for j=0:N-1
                    idx=(j*M)+i+1;   % present 1d index
                    i2=i+Dy(idx);  % shift in y direction (rows)
                    j2=j+Dx(idx);  % shift in x direction (columns)
                    if ( (i2>=0) && (i2<M) && (j2>=0) && (j2<N) )
                        Dnorm=abs(Dx(idx))+abs(Dy(idx));
                        idx2=(j2*M)+i2+1; % shifted 1d index
                        Anorm=abs(Ax(idx2))+abs(Ay(idx2));
                        Aisnan=( (isnan(Ax(idx2))) && (isnan(Ay(idx2))) ) ;
                        if ( Aisnan || (Anorm==0) )
                            % nan and 0 can be overwritten
                            Ax(idx2)=-Dx(idx);
                            Ay(idx2)=-Dy(idx);
                        elseif (Dnorm>0)
                            % mark crossings of 2 non-zero motions (adjoint vf is ambiguous)
                            Ax(idx2)=nan;
                            Ay(idx2)=2;
                            % assert: !Aisnan && Anorm!=0
                        end
                    end
                end
            end                                        
        end
        
    end
    
    %% more constructors
    methods (Hidden)
        
         function [signalMovedBack]=make_LikefromSignal3(obj,s3,opts)
            % called from make_fromSignal3
            obj.motionfname=s3.signalname;
            obj.method_estim=opts.method_estim;
            
            motionmatrix=MotionEstimation (s3.xn, obj.method_estim,opts);
            
            % motion estimation computes the vector field
            % according to our definition (NOT the adjoint field!!!)
            M=zeros(size(motionmatrix)-[0,0,0,1]);
            M(:,:,:,1)= motionmatrix(:,:,:,2);  % y
            M(:,:,:,2)= motionmatrix(:,:,:,1);  % x
            % setting the vector field to the output of MotionEstimation!
            opts0.zref=opts.zref;
            obj.set_field(M,opts0);
            if min(s3.size(1:2))>256
                obj.scale=1;  % scale automatically
            end
            
            obj.set_signal(s3);
            obj.motionfname=['Motion Est. ',obj.method_estim];
            default_weight=opts.weight_scalar;
            obj.set_weightField2scalar(default_weight);
            if nargout >0
                signalMovedBack=s3.make_like;
                signalMovedBack.xn=squeeze(motionmatrix(:,:,:,3));
            end
        end
        
    end
    
    methods (Static)
        
        function vf=make()
            % contructor with universal name (callable by parent)
            vf=MotionVF();
        end
        
        function [vf,signalMovedBack]=make_fromSignal3(s3,opts)
            % vf=~(s3,[opts]) make from 3d signal using opts.method_estim
            % as motion estimator;
            % sets weightField to scalar opts.weight_scalar (default 1).
            % opts.method_estim = 'OF' - Optical Flow
            % .method_estim = 'FSBM' - Full Search Block Matching
            requireOnly(isa(s3,'Signal3D'),'local','needs 3d signal');
            requireOnly(nargin <2 || ~isfield(opts,'method_estim')|| isempty(opts.method_estim) || ...
                ismember(opts.method_estim,{'OF','FSBM'}),'local','method_estim admitted');
            if ~exist('opts','var') || isempty(opts)
                opts=struct;
            end
            if ~isfield(opts,'method_estim') || isempty(opts.method_estim)
                opts.method_estim='OF';
            end
            if ~isfield(opts,'weight_scalar') || isempty(opts.weight_scalar)
                opts.weight_scalar=1;
            end
            if ~isfield(opts,'zref') || isempty(opts.zref)
                % opts.zref=s3.midplane;
                opts.zref=[];  % cf. set_field: zref=[] means no compensated summation
            end
            
            vf=MotionVF();
            signalMovedBack=vf.make_LikefromSignal3(s3,opts);
        end
        
        function vf= simulate(hM,hN,hL,pos,shiftsize)
            % vf2d=~() simulates 2d motion vector field
            requireOnly(nargin<4 || length(pos)==2,'local','pos is coord. pair');
            if ~exist('pos','var') || isempty(pos)
                pos=[1,1];
            end
            if ~exist('shiftsize','var') || isempty(shiftsize)
                shiftsize=[1,1];
            end
            if isscalar(shiftsize)
                shiftsize=shiftsize*[1,1];
            end
            
            rectangle=ceil(0.25*[hM, hN]);
            vf2d=zeros(hM,hN,hL,2);
            
            for j=1:hL;
                offs=(j-1)*shiftsize;
                vf2d(offs(1)+pos(1):offs(1)+pos(1)-1+rectangle(1),...
                    offs(2)+pos(2):offs(2)+pos(2)-1+rectangle(2),...
                    j,1)=shiftsize(1)*ones(rectangle);
                
                vf2d(offs(1)+pos(1):offs(1)+pos(1)-1+rectangle(1),...
                    offs(2)+pos(2):offs(2)+pos(2)-1+rectangle(2),...
                    j,2)=shiftsize(2)*ones(rectangle);
            end
            
            vf=MotionVF(vf2d(1:hM,1:hN,:,:));
            vf.method_estim='sim';
        end
        
    end
    
    %% graphics
    methods
        
        function show_NormsAnglesDiff(obj,j)
            % ~(j) show norms, angles and difference for frame pair(j,j+1)
            if ~exist('j','var')
                j=1;
            end
            prepfigure(obj.fig,[],obj.figopt);
            suptitle(obj.motionfname,12);
            subplot(2,2,1);
            obj.graph_norm(j,false);
            subplot(2,2,2);
            obj.graph_angle(j,false);
            subplot(2,2,3);
            s1=obj.diffSignal(j);
            s1.axis_type='off';
            s1.graph_signal(false);
            subplot(2,2,4);
            obj.graph_occlusions(j,false);
            
        end
        
        function graph_fieldOnFrame(obj,j, quantile, open_new,opt)
            % ~([j,quantile,open_new,opt]) show vectors on top of frame j (1)
            obj.require(nargin <2 || isempty(j) || (j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if ~exist('j','var')
                j=1;
            end
            if ~exist('quantile','var') || isempty(quantile)
                if numel(obj.mat)>128
                    quantile=0.9;  % 10% largest values
                else
                    quantile=0;  % show all
                end
            end
            if ~exist('opt','var') || isempty(opt)
                opt=struct;
            end
            if ~isfield(opt,'on_full_grid') || isempty(opt.on_full_grid)
                opt.on_full_grid= max(obj.size_space)<=64;
            end
            if ~isfield(opt,'label_on')
                opt.label_on=-1;
            end
            
            if open_new
                prepfigure(obj.fig,[],obj.figopt);
            end
            
            s2=obj.signal3.frame(j);
            s2.colorlimits=[obj.signal3.min,obj.signal3.max];
            s2.colormap_freeze=false;
            s2.ison_colorbar=false;
            s2.graph_signal(false);
            hold on;
            obj.graph_field(j,quantile,false,opt);
            hold off;
            xlim([0,obj.signal3.size(2)+1]);
            ylim([0,obj.signal3.size(1)+1]);
            caxis([obj.signal3.min, obj.signal3.max]);
            if opt.label_on==1 || (s2.numel<=9 && opt.label_on~=0)
                % show cell labels
                [x,y]=meshgrid(1:s2.size(2),1:s2.size(1));
                opt.data=1:s2.numel; opt.format='%d'; opt.fontsize=14;
                textvec(x,y,opt);
            end
            
        end
        
        function graph_field(obj,j, quantile, open_new, opt )
            % ~([j,full_grid,quantile]) show vector field of all frame pairs or of (j,j+1);
            % j .... frame index (all if missing)
            % opt.on_full_grid ... show vector field on all grid positions?
            % quantile  ... shows only quantile of longest vectors
            %
            obj.require(nargin <2 || isempty(j) || (j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if ~exist('j','var')
                j=[];
            end
            if ~exist('quantile','var') || isempty(quantile)
                if numel(obj.mat)>100
                    quantile=0.9;  % 10% largest values
                else
                    quantile=0;  % show all
                end
            end
            if ~exist('opt','var') || isempty(opt)
                opt=struct;
            end
            if ~isfield(opt,'on_full_grid') || isempty(opt.on_full_grid)
                opt.on_full_grid= max(obj.size_space)<=64;
            end
            
            if open_new
                prepfigure(obj.fig,[],obj.figopt);
            end
            maxcount=50;
            ss=obj.size;
            if opt.on_full_grid
                triage=[1,1];
            else
                triage=ceil(ss(1:2)/maxcount);
            end
            selstr=[];
            if ~opt.on_full_grid
                selstr='(on subgrid) ';
            end
            
            vf2=obj.mat;
            % show only subgrid to make sure that arrows are
            % distinguishable:
            
            [X,Y]=meshgrid(1:triage(2):ss(2),1:triage(1):ss(1));
            
            do_all_frame_pairs =isempty(j);
            
            
            if ~do_all_frame_pairs
                zstr=['z=',num2str(j)];
                vfX=vf2(1:triage(1):end,1:triage(2):end,j,2);
                vfY=vf2(1:triage(1):end,1:triage(2):end,j,1);
                nf= vfX.^2+vfY.^2;
                nfsorted=sort(nf(:),'ascend');
                L=length(nfsorted);
                cutoff_idx=min(L,max(1,ceil(quantile*L)));
                
                if L<25
                    marker='o';
                    cutoff=nfsorted(cutoff_idx);
                else
                    cutoff=max(eps,nfsorted(cutoff_idx));
                    marker='none';
                end
                nzpos=nf>=cutoff;  % show only largest values
                if ~isempty(obj.color)
                    quiver(X(nzpos),Y(nzpos),vfX(nzpos),vfY(nzpos),obj.scale,...
                        'LineWidth',obj.LineWidth,'color',obj.color,'Marker',marker);
                else
                    quiver(X(nzpos),Y(nzpos),vfX(nzpos),vfY(nzpos),obj.scale,...
                        'LineWidth',obj.LineWidth,'Marker',marker);
                end
            else
                legstr=cell(1,obj.size_time);
                L=obj.size_time;
                zstr=['z=1:',num2str(L)];
                for j=1:L
                    vfX=vf2(1:triage(1):end,1:triage(2):end,j,2);
                    vfY=vf2(1:triage(1):end,1:triage(2):end,j,1);
                    nf= vfX.^2+vfY.^2;
                    nfsorted=sort(nf(:),'ascend');
                    L=length(nfsorted);
                    cutoff_idx=min(L,max(1,ceil(quantile*L)));
                    cutoff=max(eps,nfsorted(cutoff_idx));
                    nzpos=nf>=cutoff;  % show only largest values
                    if ~isempty(obj.color)
                        quiver(X(nzpos),Y(nzpos),vfX(nzpos),vfY(nzpos),obj.scale,...
                            'LineWidth',obj.LineWidth,'color',obj.color);
                    else
                        quiver(X(nzpos),Y(nzpos),vfX(nzpos),vfY(nzpos),obj.scale,...
                            'LineWidth',obj.LineWidth);
                    end
                    legstr{j}=['z=',num2str(j)];
                    hold all;
                end
                hold off;
            end
            axis ij; % image convention
            axis image;
            xlim([1,obj.size(2)]);
            ylim([1,obj.size(1)]);
            
            if do_all_frame_pairs
                legend(legstr,'location','best');
            end
            tit=[obj.motionfname,' ',selstr,'at ',zstr];
            if ~strcmp(obj.method_estim,'sim')
                tit=[tit,' (method ',obj.method_estim];
            else
                tit=[tit,' (',obj.method_estim];
            end
            if quantile>0
                tit=[tit,', largest ',num2str(round((1-quantile)*100)),'%)'];
            else
                tit=[tit,')'];
            end
            
            title(tit,'fontsize',12);
        end
        
        
        function show_field(obj)
            % show vector field of all frame pairs
            prepfigure(obj.fig,[],obj.figopt);
            
            ss=obj.size;
            vf2=obj.mat;
            [X,Y]=meshgrid(1:ss(2),1:ss(1));
            nzcutoff=1e-9; % defines non-zero condition
            
            hL=ss(3);
            suptitle(obj.motionfname,14);
            sd=factor_subplots(hL);
            for j=1:hL
                subplot(sd(1),sd(2),j);
                vfX=vf2(:,:,j,2);
                vfY=vf2(:,:,j,1);
                nzpos=vfX.^2+vfY.^2>nzcutoff; % show only non-zero values
                if ~isempty(obj.color)
                    quiver(X(nzpos),Y(nzpos),vfX(nzpos),vfY(nzpos),obj.scale,...
                        'LineWidth',obj.LineWidth,'color','k');
                else
                    quiver(X(nzpos),Y(nzpos),vfX(nzpos),vfY(nzpos),obj.scale,...
                        'LineWidth',obj.LineWidth);
                end
                axis ij; % image convention
                xlim([1,obj.size(2)]);
                ylim([1,obj.size(1)]);
                title(['z=',num2str(j)],'fontsize',12);
            end
        end
        
        
        function graph_occlusions(obj,j,open_new)
            % ~(j) show occlusions of motion from frame j to frame j+1
            obj.require(nargin <2 || isempty(j) || (j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig,[],obj.figopt);
            end
            if nargin <2 || isempty(j)
                j=1;
            end
            occ=Signal2D(obj.occlusions(j));
            cc=sum(occ.xn(:))/numel(occ.xn);
            occ.signalname=obj.motionfname;
            occ.signalname=[num2tex(cc*100,'%3.1e','none'),'% occlusions, z=',...
                num2str(j),': ', obj.signal3.signalname];
            occ.colormap_active='gray';
            occ.axis_type='on';  % not ot confuse white border
            occ.graph_signal(false);
            
        end
        
        function graph_weightField(obj,j,open_new)
            % ~(j) show weight field of motion from frame j to frame j+1
            obj.require((j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig,[],obj.figopt);
            end
            if nargin <2 || isempty(j)
                j=1;
            end
            wf=Signal2D(obj.weightField.FieldOnDomain(j));
            wf.signalname=obj.motionfname;
            %domain_size=numel(obj.weightField.domain);
            if isempty(obj.weightField.domain)
                warning('no domain set');
            end
            wf.colormap_freeze=false;
            wf.axis_type='off';
            wf.signalname={['weight field outside occlusions, z=',num2str(j)],obj.signal3.signalname};
            wf.graph_signal(false);
        end
        
        function graph_crossings(obj,j, open_new)
            % ~() show occlusions of motion from frame j to frame j+1
            obj.require(nargin <2 || isempty(j) || (j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig,[],obj.figopt);
            end
            if nargin <2 || isempty(j)
                j=1;
            end
            cp=Signal2D(obj.crossings(j));
            cp.signalname=obj.motionfname;
            cp.fig=obj.fig;
            cc=sum(cp.xn(:))/numel(cp.xn);
            cp.signalname={[num2tex(cc*100,'%3.1e','none'),'% crossings: z=',...
                num2str(j)],obj.signal3.signalname};
            cp.axis_type='off';
            cp.graph_signal(false);
        end
        
        function graph_EulerIntegrationError(obj,open_new)
            % ~() PSNR between signal and Euler integration
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig,[],obj.figopt);
            end
            [s3,z0]=obj.eulerIntegration;
            rs=0.1;            
            L=obj.signal3.size(3);
            yn=zeros(1,L);
            N=prod(s3.size(1)*s3.size(2));
            for j=1:L
                dev=obj.signal3.frame(j).diff(s3.frame(j));
                yn(j)=dev.sparsityDefect(rs)/N;
            end
            d=TimeSignal(yn);
            %d=TimeSignal(obj.signal3.PSNR(s3));
            
            d.signalname=['error as rel. sparsity defect(',s3.signalname,') (z_0=',num2str(z0),')'];
            d.xn(z0)=NaN;
            %d.amp_quantity='PSNR';
            %d.amp_unit='dB';
            d.amp_quantity='\sigma_s(x-x_{Euler})_p';
            d.marker='.';
            d.graph_signal(false);
            xlabel('frame index');
            
        end
        
        function show_occlusions(obj)
            % ~() show occlusions of motion
            occ=Signal3D(obj.occlusions);
            occ.signalname=obj.motionfname;
            occ.fig=obj.fig;
            occ.ison_colorbar=false;
            occ.colormap_active='gray';
            occ.signalname=['occlusions: ', obj.signal3.signalname];
            occ.show_signal;
        end
        
        function show_crossings(obj)
            % ~() show occlusions of motion
            cp=Signal3D(obj.crossings);
            cp.signalname=obj.motionfname;
            cp.ison_colorbar=false;
            cp.colormap_active='gray';
            cp.fig=obj.fig;
            cp.signalname=['crossings: ', obj.signal3.signalname];
            cp.show_signal;
        end
        
        function show_normWithAngle(obj)
            % ~() show norms and angles of vector field of all frame pairs
            prepfigure(obj.fig,[],obj.figopt);
            
            p=2;
            nf=obj.norm(p);
            af=obj.angle;
            nfsorted=sort(nf(~isnan(nf)),'ascend');
            ne=numel(nfsorted);
            clims=[nfsorted(max(1,floor(0.01*ne)))-eps,nfsorted(min(ne,ceil(0.99*ne)))+eps];
            
            hL=obj.size_time;
            suptitle(['motion vector field ',obj.motionfname],14);
            sd=factor_subplots(2*hL);
            colormap('default');
            for j=1:hL
                n=2*j-1;
                subplot(sd(1),sd(2),n);
                imagescnan(nf(:,:,j),clims);
                colorbar;
                axis off; % axis ij; % image convention
                title(['norms at z=',num2str(j)],'fontsize',12);
                
                subplot(sd(1),sd(2),n+1);
                imagescnan(af(:,:,j),[-pi,pi]);
                colorbar;
                axis off; % image convention
                title(['angles at z=',num2str(j)],'fontsize',12);
                
            end
        end
        
        function graph_norm(obj,j,open_new)
            % ~(j) show norms of vector field of frame pair (j,j+1)
            obj.require(nargin <2 || isempty(j) || (j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig,[],obj.figopt);
            end
            if nargin <2 || isempty(j)
                j=1;
            end
            
            p=2;
            nf=obj.norm(p);
            colormap('default');
            imagesc(nf(:,:,j));
            colorbar;
            freezeColors;
            cbfreeze;
            axis off;
            title(['norms at z=',num2str(j)],'fontsize',12);
            
        end
        
        function graph_angle(obj,j,open_new)
            % ~(j) show angles of vector field of frame pair (j,j+1)
            obj.require(nargin <2 || isempty(j) || (j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig,[],obj.figopt);
            end
            if nargin <2 || isempty(j)
                j=1;
            end
            
            p=2;
            nf=obj.angle;
            colormap('default');
            imagescnan(nf(:,:,j));
            colorbar;
            freezeColors;
            cbfreeze;
            axis off;
            title(['angles at z=',num2str(j)],'fontsize',12);
            
        end
        
        function show_norm(obj)
            % show norms of vector field of all frame pairs
            prepfigure(obj.fig,[],obj.figopt);
            p=2;
            nf=obj.norm(p);
            nfsorted=sort(nf(~isnan(nf)),'ascend');            
            ne=numel(nfsorted);
            clims=[nfsorted(max(1,floor(0.01*ne)))-eps,nfsorted(min(ne,ceil(0.99*ne)))+eps];
            
            hL=obj.size_time;
            suptitle('norm of motion vector field',14);
            sd=factor_subplots(hL);
            colormap('default');
            for j=1:hL
                subplot(sd(1),sd(2),j);
                imagesc(nf(:,:,j),clims);
                colorbar;
                axis off; % image convention
                title(['z=',num2str(j)],'fontsize',12);
            end
        end
        
        function show_angle(obj)
            % ~(j) show probability distribution of norm and angle
            prepfigure(obj.fig,[],obj.figopt);
            
            af=obj.angle;
            
            hL=obj.size_time;
            suptitle('angle of motion vector field',14);
            sd=factor_subplots(hL);
            colormap('default');
            for j=1:hL
                subplot(sd(1),sd(2),j);
                imagescnan(af(:,:,j),[-pi,pi]);
                colorbar;
                axis off; %axis ij; % image convention
                title(['z=',num2str(j)],'fontsize',12);
            end
        end
        
        
        function show_AllAssociatedFields(obj,show_all)
            % ~([show_all]) show field and [all] associated fields
            if nargin <2
                show_all=true;
            end
            
            prepfigure(obj.fig,[],obj.figopt);
            countplots=6;
            if ~show_all
                countplots=4;
            end
            sd=factor_subplots(countplots);
            
            subplot(sd(1),sd(2),1);
            oldcolor=obj.color;
            obj.color='k'; % black arrows
            obj.graph_fieldOnFrame(1,[],false);
            
            subplot(sd(1),sd(2),3);
            obj.adjField.color='k';
            obj.adjField.graph_fieldOnFrame(1,[],false);
            
            subplot(sd(1),sd(2),4);
            obj.inverseField.graph_fieldOnFrame(1,[],false);
            
            if countplots>4
                
                subplot(sd(1),sd(2),5);
                vf=obj.adjField.inverseField;
                vf.motionfname='inverse of adjoint : -(V^*)';
                vf.graph_fieldOnFrame(1,[],false);
                
                subplot(sd(1),sd(2),6);
                vf=obj.inverseField.adjField;
                vf.motionfname='adjoint of inverse: (-V)^*';
                vf.graph_fieldOnFrame(1,[],false);
            end
            
            xlims=xlim;
            ylims=ylim;
            
            subplot(sd(1),sd(2),2);
            k=1;
            s= obj.moveForwardSignal(k);
            s.colormap_freeze=false;
            s.ison_colorbar=false;
            s.graph_signal(false);
            xlim(xlims); ylim(ylims);
            
            obj.color=oldcolor;
        end
        
        
    end
    
    
end

