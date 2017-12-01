classdef MotionField <DC
    % abstract class for motion fields operating on frames (matrices)
    
    properties (SetAccess=protected)
        mat      %@<matrix>   matrix or tensor representing vector field
    end
    
    properties (Access=protected)
        meth_comp         %@<integer>method used for motion compensation
    end
    
    properties
        method_estim  %@<string> method of motion estimation
        
        motionfname     %@<string> name of field
        color
        fig
        figopt
        signal3       %<Signal3D>  3d signal
        use_mex=true  %<logical> use mex functions if there is also a matlab function?
    end
    
    %% commands
    methods
        
        function obj= MotionField (cf)
            % constructor
            if nargin==0
                cf=[];
            end
            obj.set_field(cf);
            obj.fig=1;
            obj.figopt=struct;
            obj.figopt.pixsizeX=1000;
            obj.motionfname='sim';
            
            obj.method_estim='sim';
        end
        
        function set_field(obj,cf)
            % ~(cf,[M]) set the field
            obj.mat=cf;
        end
                
        function set_defaultfield(obj,msize)
            % ~(cf,[M]) set 0s as default field
            obj.mat=zeros([msize,2]);
        end
        
        function set_signal(obj,s3)
            % ~(s3) set 3d signal
            obj.requireOnly(isa(s3,'Signal3D'),'local','3d signal');
            obj.signal3=s3;
        end
        
        function reset_motionField(obj)
            % ~() reset motion fields
            obj.mat=[];
        end
        
    end
    
    %% abstract methods
    
    methods (Abstract)
        
        y=moveForward(obj,j,x)  % move-forward using field
        y=moveBackward(obj,j,x) % move-backward
        y=moveForwardAdj(obj,j,x) % adjoint of move-forward
        y=moveBackwardAdj(obj,j,x) % adjoint of move-backward
        qf=inverseField(obj)    % inverse field for move-backward
        
    end
    
    methods (Abstract, Static)
        qf=make()               % contructor with universal name
        qf=make_fromSignal3(s3,method_estim,opts)  % constructutor
        vf=test_SimpleIllustration(params)    % constructs example
    end
    
    %% transformations
    methods
        
        function s=moveBackwardSignal(obj,j)
            % s=~(j) move backward frame j using present vector field
            s=Signal2D(obj.moveBackward(j,obj.signal3.xn(:,:,j)));
            % to make colors comparable:
            s.colorlimits=[obj.signal3.min,obj.signal3.max];
            s.signalname='image moved backward along (V,Q)';
        end
        
        function s=moveForwardSignal(obj,j)
            % s=~(j) move backward frame j using present vector field
            s=Signal2D(obj.moveForward(j,obj.signal3.xn(:,:,j)));
            % to make colors comparable:
            s.colorlimits=[obj.signal3.min,obj.signal3.max];
            s.signalname={'image moved forward along (V,Q)',obj.method_compensate};
        end
        
        function M=matrix(obj,moveOp,j)
            % M= ~(moveOp,j) matrix representation of function moveOp
            % columns of M are the result of applying moveOp to the
            % natural basis vectors e_i.
            % Should be only called with small image sizes.
            % ex: cf=MotionVF.test_SimpleIllustration();
            % MF=cf.matrix(@cf.moveForward,1);
            obj.require(isa(moveOp,'function_handle'),'local','moveOp is function');
            obj.requireOnly(j>=1 && j<=obj.size_time,'local','index admitted');
            oldsignal=obj.signal3;
            SP=Signal2D();
            
            hN=obj.size_space;
            hnum=prod(hN);
            M=zeros(hnum,hnum);
            wait_handle = waitbar(0,...
                ['computing matrix representation (dim=',num2str(hnum),') ...']);
            for i=1:hnum
                waitbar(i / hnum);
                ei=SP.make_delta(i,hN,0); % next vector of natural basis
                obj.signal3.xn(:,:,j)=ei.xn;
                obj.signal3.xn(:,:,j+1)=ei.xn;  % for backward moves
                s=moveOp(j);
                % set row vector of matrix M:
                M(:,i)=s(:);
            end
            close(wait_handle );
            
            obj.set_signal(oldsignal);
            
        end
        
    end
    
    %% queries
    methods
        
        function vf2=clone(obj)
            % clone object
            obj.require(true,'local','check invariant');
            vf2=obj.make_like();  % make same class
            vf2.mat=obj.mat;  
            obj.ensureOnly(~eq(obj,vf2),'local','new object created');
        end
        
        function vf2=clone_deep(obj)
            % deep-clone object
            vf2=obj.clone;
            obj.ensureOnly(~eq(obj,vf2),'local','new object created');
        end
        
        function vf2=make_like(obj)
            % vf2=~() same as clone but without matrix content
            vf2=obj.make();  % make same class
            vf2.mat=[];
            vf2.color=obj.color;
            vf2.motionfname =obj.motionfname;
            vf2.fig=obj.fig;
            vf2.figopt=obj.figopt;
            vf2.signal3 =obj.signal3;
            vf2.use_mex=obj.use_mex;
            obj.ensureOnly(~eq(vf2,obj),'local','new object created');
        end
        
        function s=frame(obj,j)
            % s=~(j) retrieves frame number j from mat
            obj.requireOnly(min(j)>=1 && max(j)<=obj.size_time,'local',...
                'j admissible');
            s=obj.clone_deep;
            s.mat=s.mat(:,:,j);
        end
        
        function ok=sizes_match(obj)
            % ok=~() chekc whether singla size and motion field size match
            ok=isequal(obj.signal3.size(1:2), obj.size_space) && ...
                obj.size_time==obj.signal3.size(3);
        end
        
        function str= method_compensate(obj)
            % str=~() active method used for motion compensation
            str=obj.meth_comp;
        end
        
        function ss=size(obj,d)
            % s=~([d]) size along dimension d
            %ss=[obj.M, size(obj.mat,1)/obj.M,size(obj.mat,2)];
            ss= size(obj.mat);
            if nargin>1
                ss=ss(d);
            end
        end
        
        function L=size_time(obj)
            L=size(obj.mat,3);
        end
        
        function N=size_space(obj)
            N=[size(obj.mat,1), size(obj.mat,2)];
        end
        
        function d=delta(obj,j)
            % x=~(j) difference MATRIX: frame(j+1)-frame(j) USING motion
            % compensation
            obj.requireOnly(j>=1 && j<=obj.size_time,'local',' j permitted');
            obj.requireOnly(isa(obj.signal3,'Signal3D'),'local','2d signal is set');
            % which one is best ? Does that depend on the motion estimator?
            % For MotionEstimation with 'OF' movebackward seems to be
            % better!
            d=obj.signal3.xn(:,:,j+1)-obj.moveForward(j);
            %d=obj.moveBackward(j,obj.signal3.xn(:,:,j+1))-obj.signal3.xn(:,:,j);
        end
        
        function s2=deltaSignal(obj,j)
            % s=~(j) difference SIGNAL: frame(j+1)-frame(j) USING motion
            % compensation.
            obj.requireOnly(j>=1 && j<=obj.size_time,'local',' j permitted');
            obj.requireOnly(isa(obj.signal3,'Signal3D'),'local','3d signal is set');
            d=obj.delta(j);
            s2=Signal2D(d);
            s2.signalname=['motion compensated \Delta_{\tau}(',num2str(j+1),',',num2str(j),') of ',...
                obj.signal3.signalname];
        end
        
        function sd=sparsityDefect_delta(obj,rs)
            % d=~(rs) rel. sparsity defect of motion compensated diff.
            % for relative sparsity rs in [0,1].
            % a.k.a. lp/error of best s/term approximation of delta
            if nargin <2
                rs=0.1; % take sparsity s=rs*lengthOfSignal
            end
            L=obj.size_time;
            p=1;
            sd=zeros(1,L);
            cc=numel(obj.signal3.xn(:,:,1));
            % absolute sparsity:
            s=max(1,floor(rs*cc));
            for j=1:L
                xstar=sort(reshape(abs(obj.delta(j)),[],1),'descend');                                              
                sd(j)=norm(xstar(s:end),p)/cc;                
            end
        end
        
        function [sd,p]=sparsityDefect_diff(obj,rs)
            % d=~(rs) rel. sparsity defect of uncompensated diff.
            % for relative sparsity rs in [0,1].
            % a.k.a. lp/error of best s/term approximation of difference
            if nargin <2
                rs=0.1; % take sparsity s=rs*lengthOfSignal
            end
            L=obj.size_time;
            p=1;
            sd=zeros(1,L);
            cc=numel(obj.signal3.xn(:,:,1));
            % absolute sparsity:
            s=max(1,floor(rs*cc));
            for j=1:L
                xstar=sort(reshape(abs(obj.diff(j)),[],1),'descend');                             
                sd(j)=norm(xstar(s:end),p)/cc;                
            end
        end
        
        function v=norm_delta(obj,p)
            % v=~([p]) compute rel. lp-norm of motion compensated difference
            % images
            obj.requireOnly(isa(obj.signal3,'Signal3D'),'local','3d signal is set');
            if nargin <2
                p=1;
            end
            L=obj.size_time;
            v=zeros(1,L);
            for j=1:L
                v(j)=norm(reshape(obj.delta(j),[],1),p)/...
                    norm(reshape(obj.signal3.frame(j),[],1),p);
            end
        end
        
        function [v,l1rel]=totalVarDelta(obj)
            % v=~(dim) 1-d pixelwise  total variation of motion compensated delta
            % optional output: l1-norm
            obj.require(isa(obj.signal3,'Signal3D'),'local','3d signal is set');
            v=0;
            L=obj.size_time;
            for j=1:L
                v= v+abs(obj.delta(j));
            end
            if nargout>1
                L=obj.signal3.size(3);
                n1=obj.signal3.norm(1)*(L-1)/L;
                l1rel=sum(v(:))/n1;
            end
        end
        
        function [v,l1rel]=totalVarDiff(obj)
            % v=~() 1-d pixelwise total variation of diff
            % optional output: l1-norm
            obj.require(isa(obj.signal3,'Signal3D'),'local','3d signal is set');
            if nargin<2
                dim=3;
            end
            v= sum(abs(diff(obj.signal3.xn,1,dim)),dim);
            if nargout>1
                L=obj.signal3.size(3);
                n1=obj.signal3.norm(1)*(L-1)/L;
                l1rel=sum(v(:))/n1;
            end
        end
        
        function d=diff(obj,j)
            % x=~(j) difference MATRIX: frame(j+1)-frame(j) uncompensated
            obj.requireOnly(j>=1 && j<=obj.size_time,'local',' j permitted');
            obj.requireOnly(isa(obj.signal3,'Signal3D'),'local','2d signal is set');
            
            d=obj.signal3.xn(:,:,j+1)-obj.signal3.xn(:,:,j);
        end
        
        function s2=diffSignal(obj,j)
            % s=~(j) difference signal: frame(j+1)-frame(j) WITHOUT using motion
            obj.requireOnly(isa(obj.signal3,'Signal3D'),'local','3d signal is set');
            v=obj.signal3.xn(:,:,j+1);
            w=obj.signal3.xn(:,:,j);
            % subtract signal values
            s2=Signal2D(v-w);
            s2.signalname=['D(',num2str(j+1),',',num2str(j),') of ',...
                obj.signal3.signalname];
        end
        
        function v=norm_diff(obj,p)
            % v=~([p]) compute rel. lp-norm of uncompensated difference
            % images
            obj.requireOnly(isa(obj.signal3,'Signal3D'),'local','3d signal is set');
            if nargin <2
                p=1;
            end
            L=obj.size_time;
            v=zeros(1,L);
            for j=1:L
                v(j)=norm(reshape(obj.diff(j),[],1),p)/...
                    norm(reshape(obj.signal3.frame(j),[],1),p);
            end
        end
        
        function s2=diffWithIC(obj, zref)
            % s2=~(zref) form differences with signal at j=zref as initial condition
            obj.requireOnly(isa(obj.signal3,'Signal3D'),'local','2d signal is set');
            s2=obj.signal3.clone;
            L=s2.size(3);
            for j=zref+1:L
                s2.xn(:,:,j)=obj.signal3.xn(:,:,j)-obj.signal3.xn(:,:,j-1);
            end
            for j=zref-1:-1:1
                s2.xn(:,:,j)=obj.signal3.xn(:,:,j+1)-obj.signal3.xn(:,:,j);
            end
            s2.signalname=['extended diff_{IC} ',obj.signal3.signalname];
            
        end
        
        function s2= integrateDiff(obj,diffSig,zref)
            % s2=~(delta,zref) integrate difference signal with initial
            % condition at zref
            obj.requireOnly(isa(diffSig,'Signal3D'),'local','needs delta signal');
            % save original signal
            sorig=obj.signal3.clone;
            % replace signal by delta-signal:
            obj.signal3=diffSig.clone;
            % integration of deltaSig operates on obj.signal3:
            L=obj.signal3.size(3);
            % we leave obj.signal.xn(:,:,objzref) as it is.
            for j=zref+1:L
                obj.signal3.xn(:,:,j)=obj.signal3.xn(:,:,j)+obj.signal3.xn(:,:,j-1);
            end
            for j=zref-1:-1:1
                obj.signal3.xn(:,:,j)=obj.signal3.xn(:,:,j+1)-obj.signal3.xn(:,:,j);
            end
            obj.signal3.signalname=['Integrated ',obj.signal3.signalname];
            s2=obj.signal3;
            % restore original signal
            obj.signal3=sorig;
        end
        
        
        function s2=deltaWithIC(obj, zref)
            % s2=~(zref) form motion compensated differences with signal at j=zref as initial condition
            obj.requireOnly(isa(obj.signal3,'Signal3D'),'local','2d signal is set');
            s2=obj.signal3.clone;
            L=s2.size(3);
            for j=zref+1:L
                s2.xn(:,:,j)=obj.signal3.xn(:,:,j)-obj.moveForward(j-1);
            end
            for j=zref-1:-1:1
                s2.xn(:,:,j)=obj.moveBackward(j+1)-obj.signal3.xn(:,:,j);
            end
            s2.signalname=['extended \Delta_{IC}: ',obj.signal3.signalname];
        end
        
        function s2= integrateDelta(obj,deltaSig,zref)
            % s2=~(delta,zref) integrate delta signal with initial
            % condition at zref
            obj.requireOnly(isa(deltaSig,'Signal3D'),'local','needs delta signal');
            % save original signal
            sorig=obj.signal3.clone;
            % replace signal by delta-signal:
            obj.signal3=deltaSig.clone;
            % integration of deltaSig operates on obj.signal3:
            L=obj.signal3.size(3);
            % we leave obj.signal.xn(:,:,objzref) as it is.
            for j=zref+1:L
                obj.signal3.xn(:,:,j)=obj.signal3.xn(:,:,j)+obj.moveForward(j-1);
            end
            for j=zref-1:-1:1
                obj.signal3.xn(:,:,j)=obj.moveBackward(j+1)-obj.signal3.xn(:,:,j);
            end
            obj.signal3.signalname=['Integrated ',obj.signal3.signalname];
            s2=obj.signal3;
            % restore original signal
            obj.signal3=sorig;
        end
        
        
    end
    
    %% test
    
    methods
        
        function test_motionCompensation(obj,j)
            % ~([j]) test motion compensation of frame j (default all)
            obj.require(isa(obj.signal3,'Signal3D'),'local','2d signal is set');
            obj.requireOnly(nargin<2 || (j>=1 && j<=obj.size_time),'local',' j permitted');
            if nargin <2
                j=[];
            end
            compute_all=isempty(j);
            if compute_all
                L=obj.size_time;
                n1=zeros(L,1);
                n2=n1;
                for j=1:L
                    s1=obj.diffSignal(j);
                    s2=obj.deltaSignal(j);
                    fjnorm=obj.signal3.frame(j).norm(1);
                    n1(j)=s1.norm(1)/fjnorm;
                    n2(j)=s2.norm(1)/fjnorm;
                end
                mean1=mean(n1);
                mean2=mean(n2);
            else
                d1=obj.diffSignal(j);
                d2=obj.deltaSignal(j);
                fjnorm=obj.signal3.frame(j).norm(1);
                l1rel_d1=d1.norm(1)/fjnorm;
                l1rel_d2=d2.norm(1)/fjnorm;
            end
            prepfigure(obj.fig,[],obj.figopt);
            
            if compute_all
                plot(1:L,n1,'o-',1:L,n2,'x-');
                xlabel('time');
                ylabel('||:||_{1,rel}');
                legend(['diff, \mu=',num2str(mean1,'%3.1e')],...
                    ['delta, \mu=',num2str(mean2,'%3.1e'),...
                    sprintf('\n motion compensated')],'location','best');
                tit=['compensation method=',obj.method_compensate,', ',obj.signal3.get_signalname];
                title(tit,'fontsize',12);
            else
                
                suptitle(['compensation method=',obj.method_compensate],14);
                
                subplot(2,2,1);
                d1.colormap_freeze=false;
                d1.graph_signal(false);
                h=get(gca,'Title');
                tit1=get(h,'String');
                tit={tit1, ['||.||_{1,rel}=',num2str(l1rel_d1,'%3.1e')]};
                title(tit,'fontsize',12);
                
                subplot(2,2,3);
                d2.colormap_freeze=false;
                d2.graph_signal(false);
                h=get(gca,'Title');
                tit1=get(h,'String');
                tit={'motion compensated', ['||.||_{1,rel}=',num2str(l1rel_d2,'%3.1e')]};
                title(tit,'fontsize',12);
                
                subplot(2,2,[2,4]);
                d1.graph_kde(false);
                hold all;
                d2.graph_kde(false);
                hold off;
                title('kernel density estimator');
                legend(sprintf('w/o motion \ncompensation'),'with ~','location','NorthWest');
            end
            
        end
        
    end
    
    methods (Static)
        
        function [ok,MF,MFA,MB,MBA]=test_adjointMoves(params)
            % ok=~([isrand]) test if adjoint moves are actually adjoint to moves
            % using random vector fields if isrand=true;
            % e.g. [ok,MF,MFA,MB,MBA]=MotionVF.test_adjointMoves(false);
            %
            if nargin <1 || isempty(params)
                params=struct;
            end
            if ~isfield(params,'isrand')
                params.isrand=true;  % random vector fields
            end
            if ~isfield(params,'M') || isempty(params.M)
                params.M=2; params.N=3;
            end
            
            vf=MotionVF.test_SimpleIllustration(params);
            % matrix representation of moveForward and moveForwardAdj:
            MF=vf.matrix(@vf.moveForward,1);
            MFA=vf.matrix(@vf.moveForwardAdj,1);
            % is moveForwardAdj adjoint to moveForward?
            ok1=max(abs(reshape(MF'-MFA,[],1)))<1e-9;
            
            % matrix representation of moveBackward and moveBackwardAdj:
            MB=vf.matrix(@vf.moveBackward,1);
            MBA=vf.matrix(@vf.moveBackwardAdj,1);
            % is moveBackwardAdj adjoint to moveBackward?
            ok2=max(abs(reshape(MB'-MBA,[],1)))<1e-9;
            
            ok=ok1 & ok2;
        end
        
        function [x,xforw,xbackw,vf]=test_moves()
            % ~() test moveForward and moveBackward with simple example
            params.rand=true;
            vf=test_SimpleIllustration(params);
            xforw=vf.moveForward(1,x);
            
            xbackw=vf.moveBackward(2,y);
            
        end
        
    end
    
    %% more constructors
    methods (Static)
        
        
        
    end
    
    %% transforms
    methods
        
        
    end
    
    %% graphics
    methods
        
        function show_field_TV(obj,j,quantile)
            % ~([j,quantile]) show vector field (j or all) and total variation
            if ~exist('quantile','var')
                quantile=[];
            end
            if ~exist('j','var')
                j=[];
            end
            prepfigure(obj.fig,[],obj.figopt);
            subplot(1,2,1);
            obj.graph_field(j,quantile,false);
            subplot(1,2,2);
            obj.graph_totalVarDiff(false);
            
        end
        
        function graph_totalVarDiff(obj, open_new, clims)
            % ~([open_new, clims]) show total variation of diff operator (not motion compensated)
            obj.require(isa(obj.signal3,'Signal3D'),'local','3d signal is set');
            if ~exist('open_new','var') || isempty(open_new)
                prepfigure(obj.fig,[],obj.figopt);
            end
            if ~exist('clims','var')
                clims=[];
            end
            [img,l1]=obj.totalVarDiff();
            linf=max(img(:));
            s2=Signal2D(img);
            s2.signalname=['TV(D) of ',obj.signal3.signalname,', ',...
                vec2str(obj.signal3.size)];
            s2.colorlimits=clims;
            s2.colormap_freeze=false;
            s2.graph_signal(false);
            h=get(gca,'Title');
            present_titstr{1}=get(h,'String');
            present_titstr{2}=['|||.||_{1,rel}=', num2str(l1,'%3.1e'),', ||.||_\infty=',...
                num2str(linf,'%3.1e')];
            title(present_titstr,'fontsize',12);
            
        end
        
        function graph_totalVarDelta(obj, open_new, clims)
            % ~([open_new,cax]) show total variation of diff operator (not motion compensated)
            obj.require(isa(obj.signal3,'Signal3D'),'local','3d signal is set');
            if ~exist('open_new','var') || isempty(open_new)
                prepfigure(obj.fig,[],obj.figopt);
            end
            if ~exist('clims','var')
                clims=[];
            end
            [img,l1]=obj.totalVarDelta();
            linf=max(img(:));
            s2=Signal2D(img);
            s2.signalname=['motion compensated TV(\Delta_\tau) of ',obj.signal3.signalname];
            s2.colorlimits=clims;
            s2.colormap_freeze=false;
            s2.graph_signal(false);
            h=get(gca,'Title');
            present_titstr{1}=get(h,'String');
            present_titstr{2}=['|||.||_{1,rel}=', num2str(l1,'%3.1e'),', ||.||_\infty=',...
                num2str(linf,'%3.1e')];
            title(present_titstr,'fontsize',12);
            
        end
        
        
        function show_motionCompensation(obj)
            % ~() show effect of motion compensation on differences and
            % total variation
            obj.requireOnly(isa(obj.signal3,'Signal3D'),'local','3d signal is set');
            sdiff=obj.diffSignal(1);
            sdiff.colormap_freeze=false;
            sdelta=obj.deltaSignal(1);
            sdelta.colormap_freeze=false;
            
            prepfigure(obj.fig,[],obj.figopt);
            suptitle(['compensation method=',obj.method_compensate],14);
            
            subplot(2,2,1);
            sdiff.graph_signal(false);
            cax= caxis;
            subplot(2,2,2);
            sdelta.colorlimits=cax;
            sdelta.graph_signal(false);
            subplot(2,2,3);
            obj.graph_totalVarDiff(false);
            cax= caxis;
            subplot(2,2,4);
            obj.graph_totalVarDelta(false,cax);
            
        end
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal size match';
            ok= isempty(obj.signal3) || isequal(obj.signal3.size(1:2),obj.size_space);
        end
    end
    
end


