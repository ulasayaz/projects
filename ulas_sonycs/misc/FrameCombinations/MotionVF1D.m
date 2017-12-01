classdef MotionVF1D <DC
    % 1d-motion vector field operating on vectors
    % Example
    % ========
    % --- motion simulations:
    %     countI=1;
    %     N=100;L=20;
    %     params=struct;
    %     params.w=3; params.shiftsize=7; params.sig=0;
    %     params.w=7; params.shiftsize=3; params.sig=0;
    %     [s2,motion]=Signal2D.make_MovingIntervals(countI,N,L,params);
    %     s2.fig=3;
    %     -- show frames in subplots
    %     s2.graph_signal;
    %
    %     motion.fig=3;
    %     --- show vector field in 1 plot:
    %     motion.graph_field;
    %     motion.test_motionCompensation;
    %
    
    
    properties (SetAccess=protected)
        mat   %@<matrix>   matrix representing 1-d vector field
        matBack % backward vector field
    end
    
    properties
        fig
        figopt
        signal2  %<Signal2D>  2d signal
    end
    
    %% commands
    methods
        
        function obj= MotionVF1D (vf)
            % constructor
            assert(nargin <1 || ismatrix(vf),'vf dimensions ok');
            if nargin==0
                hM=10;
                vf=zeros(hM,hL);
            end
            
            obj.set_motionVF(vf);
            obj.fig=1;
            obj.figopt=struct;
            obj.figopt.pixsizeX=1000;
            
            obj.ensure(true,'local','check invariant');
        end
        
        function set_motionVF(obj,vf)
            % ~(vf,[M]) set the vector field
            obj.requireOnly(ismatrix(vf),'local','vf dimensions ok');
            obj.mat=vf;
            obj.matBack=[];                   
        end
        
        function compute_backwardField(obj)
            % ~() computes field of vectors backwards from target to source
            % result will only be well-defined if each
            % the targets of the forward field are unique.
            obj.requireOnly(~isempty(ibj.mat),'local','needs forward field');
            N=obj.size_space;
            L=obj.size_time;
            obj.matBack=zeros(size(obj.mat));  
            for k=1:L
                for j=1:N
                    target=j+obj.mat(k,j);
                    if target>=1 && target <=N
                        obj.matBack(k,target)=-obj.mat(k,j);
                    end
                end
            end
        end
        
        function set_signal(obj,s2)
            % ~(s2) set 2d signal
            obj.requireOnly(isa(s2,'Signal2D'),'local','2d signal');
            obj.signal2=s2;
            
        end
        
        
    end
    
    %% queries
    methods
        
        
        function ss=size(obj,d)
            % s=~([d]) size along dimension d
            ss=size(obj.mat);
            if nargin>1
                ss=ss(d);
            end
        end
        
        function L=size_time(obj)
            L=size(obj.mat,1);
        end
        
        function M=size_space(obj)
            M=size(obj.mat,2);
        end
        
        function flow=motion2flow(obj)
            % flow=~() motion represented as indices
            
            hM=obj.size_space;
            hL=obj.size_time;
            
            flow=obj.mat;
            plane1=1:hM;
            for j=1:hL
                flow(:,j)=plane1+flow(:,j);
            end
        end
        
        function y=moveBackward(obj,j)
            % s2=~(j) compensate motion of j-th frame to get a
            % postdiction for the (j-1)-th frame.
            obj.requireOnly(j>=2 && j<=obj.size_time,'local',' j permitted');
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            % use only information from time step j:
            x=obj.signal2.xn(j,:);
            % use backward field on the basis of time j !!! not on the basis of
            % time j-1 !!!
            motionV=obj.mat(j,:); 
                      
            N=obj.size_space;
            y=x; % prior is nothing moved
            idx=find(motionV);  % limit to non zero motion 
            for k=idx
                v=motionV(k);
                target=k-v;
                if target>=1 && target<=N
                    % decide how to replace position k
                    % which will become vacant due to motion:
                    backtarget=min(N,max(1,k+v));
                    yt=x(backtarget);
                    % move value of k to target
                    y(target)=x(k);
                    % replace value of vacant position k
                    y(k)=yt;
                end
            end
%             motionV=obj.mat(j,:);
%             N=obj.signal2.N;
%             y=zeros(size(x));
%             % exchange moved locations against
%             idxMotion=find(motionV);
%             idx=(1:N(2))-motionV;  % indices moved backward
%             idxBack=min(N(2),max(1,idxMotion-motionV(idxMotion)));
%             idx(idxBack)=idxMotion; % exchange of indices
%             % assign permuted values:
%             idx=min(N(2),max(1,idx));
%             w(idx)=v;
            % test values: w should be closer to xn(j-1) than v:
            % test=[obj.signal2.xn(j-1,:)',x(:),y(:),motionV'];
        end
        
        function y=moveForward(obj,j)
            % s2=~(j) compensate motion of j-th frame to get a
            % prediction for the (j+1)-th frame.
            obj.requireOnly(j>=1 && j<=obj.size_time,'local',' j permitted');
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            % use only information from time step j:
            x=obj.signal2.xn(j,:);
            motionV=obj.mat(j,:);    % use forward field on the basis of knowing time j
            N=obj.size_space;
            y=x; % prior is nothing moved
            idx=find(motionV);  % limit to non zero motion 
            for k=idx
                v=motionV(k);
                target=k+v;
                if target>=1 && target<=N
                    % decide how to replace position k
                    % which will become vacant due to motion:
                    backtarget=min(N,max(1,k-v));
                    yt=x(backtarget);
                    % move value of k to target
                    y(target)=x(k);
                    % replace value of vacant position k
                    y(k)=yt;
                end
            end
            
%             N=obj.signal2.N;
%             y=zeros(size(x));
%             idxMotion=find(motionV);
%             idx=(1:N(2))+motionV; % indices moved forward
%             idxForward=min(N(2),max(1,idxMotion+motionV(idxMotion)));
%             idx(idxForward)=idxMotion; % exchange of indices
%             % assign permuted values:
%             idx=min(N(2),max(1,idx));
%             y(idx)=x;
            % test values: w should be closer to xn(j+1) than v:
            % test=[obj.signal2.xn(j+1,:)',x(:),y(:),motionV'];
            
        end
        
        function d=delta(obj,j)
            % x=~(j) difference MATRIX: frame(j+1)-frame(j) USING motion
            % occluded areas get values copied from foreground instead of NaN.
            obj.requireOnly(j>=1 && j<=obj.size_time,'local',' j permitted');
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            d=obj.signal2.xn(j+1,:)-obj.moveForward(j);
        end
        
        function s2=deltaSignal(obj,j)
            % s=~(j) difference SIGNAL: frame(j+1)-frame(j) USING motion
            % occluded areas get values of foreground instead of NaN.
            obj.requireOnly(j>=1 && j<=obj.size_time,'local',' j permitted');
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            d=obj.delta(j);
            s2=TimeSignal(d);
            s2.signalname=['motion compensated \Delta_{\tau}(',num2str(j+1),',',num2str(j),') of ',...
                obj.signal2.signalname];
        end
        
        function [v,l1]=totalVarDelta(obj)
            % v=~(dim) 1-d pixelwise  total variation of motion compensated delta
            % optional output: l1-norm
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            v=0;
            L=obj.size_time;
            for j=1:L
                v= v+abs(obj.delta(j));
            end
            if nargout>1
                l1=sum(v(:));
            end
        end
        
        function [v,l1]=totalVarDiff(obj)
            % v=~() 1-d pixelwise total variation of diff
            % optional output: l1-norm
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            if nargin<2
                dim=2;
            end
            v= sum(abs(diff(obj.signal2.xn,1,dim)),dim);
            if nargout>1
                l1=sum(v(:));
            end
        end
        
        function s2=diffSignal(obj,j)
            % s=~(j) difference signal: frame(j+1)-frame(j) WITHOUT using motion
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            obj.requireOnly(j>=1 && j<=obj.size_time,'local',' j permitted');
            v=obj.signal2.xn(j+1,:);
            w=obj.signal2.xn(j,:);
            % subtract signal values
            s2=TimeSignal(v-w);
            s2.signalname=['\Delta(',num2str(j+1),',',num2str(j),') of ',...
                obj.signal2.signalname];
        end
        
        function s2=diffWithIC(obj, zref)
            % s2=~(zref) form differences with signal at j=zref as initial condition
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            s2=obj.signal2.clone;
            L=s2.size(1);            
            for j=zref+1:L
               s2.xn(j,:)=obj.signal2.xn(j,:)-obj.signal2.xn(j-1,:);
            end            
            for j=zref-1:-1:1
               s2.xn(j,:)=obj.signal2.xn(j+1,:)-obj.signal2.xn(j,:); 
            end
            s2.signalname=['diff-signal: ',obj.signal2.signalname];
            
        end
        
        function s2= integrateDiff(obj,diffSig2,zref)
             % s2=~(delta,zref) integrate difference signal with initial
            % condition at zref
            obj.requireOnly(isa(diffSig2,'Signal2D'),'local','needs delta signal');
            % save original signal
            sorig=obj.signal2.clone;
            % replace signal by delta-signal:
            obj.signal2=diffSig2.clone;
            % integration of deltaSig operates on obj.signal2:
            L=obj.signal2.size(1);            
            for j=zref+1:L
               obj.signal2.xn(j,:)=obj.signal2.xn(j,:)+obj.signal2.xn(j-1,:);              
            end            
            for j=zref-1:-1:1
               obj.signal2.xn(j,:)=obj.signal2.xn(j+1,:)-obj.signal2.xn(j,:);                 
            end            
            obj.signal2.signalname=['Integrated ',obj.signal2.signalname];
            s2=obj.signal2;
            % restore original signal
            obj.signal2=sorig;
        end
        
        
        function s2=deltaWithIC(obj, zref)
            % s2=~(zref) form motion compensated differences with signal at j=zref as initial condition
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            s2=obj.signal2.clone;
            L=s2.size(1);            
            for j=zref+1:L
               s2.xn(j,:)=obj.signal2.xn(j,:)-obj.moveForward(j-1);              
            end            
            for j=zref-1:-1:1
               s2.xn(j,:)=obj.moveBackward(j+1)-obj.signal2.xn(j,:);                 
            end            
            s2.signalname=['\Delta-signal: ',obj.signal2.signalname];
        end
        
        function s2= integrateDelta(obj,deltaSig2,zref)
            % s2=~(delta,zref) integrate delta signal with initial
            % condition at zref
            obj.requireOnly(isa(deltaSig2,'Signal2D'),'local','needs delta signal');
            % save original signal
            sorig=obj.signal2.clone;
            % replace signal by delta-signal:
            obj.signal2=deltaSig2.clone;
            % integration of deltaSig operates on obj.signal2:
            L=obj.signal2.size(1);           
            for j=zref+1:L
               obj.signal2.xn(j,:)=obj.signal2.xn(j,:)+obj.moveForward(j-1);              
            end            
            for j=zref-1:-1:1
               obj.signal2.xn(j,:)=obj.moveBackward(j+1)-obj.signal2.xn(j,:);                 
            end            
            obj.signal2.signalname=['Integrated ',obj.signal2.signalname];
            s2=obj.signal2;
            % restore original signal
            obj.signal2=sorig;
        end
        
        function nf = norm(obj)
            % val=~() in 1d all lp norms are the same for p<infty
            nf=abs(obj.mat);
        end
        
        function val=meannz(obj, cutoff)
            % val=~() mean length of non-zero vectors
            if nargin <2
                cutoff=1e-9;
            end
            nf=obj.norm();
            nf=nf(nf>cutoff);
            val=mean(nf);
        end
        
    end
    
    %% test
    
    methods 
        
        function test_motionCompensation(obj,j)
            % ~(j) test motion compensation
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
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
                    n1(j)=s1.norm(1);
                    n2(j)=s2.norm(1);
                end
            else
                d1=obj.diffSignal(j);
                d2=obj.deltaSignal(j);
            end
            prepfigure(obj.fig,[],obj.figopt);
            
            if compute_all                
                plot(1:L,n1,'o-',1:L,n2,'x-');
                xlabel('time');
                ylabel('l_1-norm');
                legend('diff','delta motion compensated');
                tit=obj.signal2.get_signalname;
            else
                
                subplot(1,2,1);
                d1.graph_signal(false);
                h=get(gca,'Title');
                tit1=get(h,'String');
                tit={tit1, ['||.||_1=',num2str(d1.norm(1),'%3.1e')]};
                title(tit,'fontsize',12);
                
                subplot(1,2,2);
                d2.graph_signal(false);
                h=get(gca,'Title');
                tit1=get(h,'String');
                tit={tit1, ['||.||_1=',num2str(d2.norm(1),'%3.1e')]};                
            end
            
            title(tit,'fontsize',12);
        end
        
        
    end
    
    %% transforms
    methods
        
        
       
        
    end
    
    %% graphics
    methods
        
        function graph_field(obj,j, on_full_grid, open_new)
            % ~([j]) show vector field of all frame pairs or of (j,j+1);
            obj.require(nargin <2 || isempty(j) || (j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if ~exist('j','var')
                j=[];
            end
            if ~exist('on_full_grid','var')
                on_full_grid=true;
            end
            
            if open_new
                prepfigure(obj.fig,[],obj.figopt);
            end
            maxcount=50;
            ss=obj.size;
            if on_full_grid
                triage=1;
            else
                triage=ceil(ss(1)/maxcount);
            end
            selstr=[];
            if ~on_full_grid
                selstr='(on subgrid) ';
            end
            
            vf2=obj.mat;
            nzcutoff=1e-9; % defines non-zero condition
            do_all_frame_pairs =isempty(j);
            
            if ~do_all_frame_pairs
                X=1:triage:ss(2);
                zstr=['z=',num2str(j)];
                vfX=vf2(j,1:triage(1):end);
                nzpos=abs(vfX)>nzcutoff;  % show only non-zero values
                xX=X(nzpos);
                quiver(X(nzpos),zeros(size(xX)),vfX(nzpos),zeros(size(xX)));
            else
                L=obj.size_time;
                % show only subgrid to make sure that arrows are
                % distinguishable:
                [X,Y]=meshgrid(1:triage:ss(2),1:L);
                zstr=['z=1:',num2str(L)];
                vfX=vf2(:,1:triage(1):end);
                nzpos=abs(vfX)>nzcutoff;
                vfY=zeros(size(vfX));
                quiver(X(nzpos),Y(nzpos), vfX(nzpos),vfY(nzpos));
            end
            axis ij; % image convention
            xlim([1,obj.size(2)]);
            ylim([1,obj.size(1)]);
            xlabel('z');
            
            title(['vector field ',selstr,'at ',zstr],'fontsize',12);
        end
        
        function graph_totalVarDiff(obj, open_new)
            % ~() show total variation of diff operator (not motion compensated)
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            if ~exist('open_new','var') || isempty(open_new)
                prepfigure(obj.fig,[],obj.figopt);
            end
            [img,l1]=obj.totalVarDiff();
            s2=TimeSignal(img);
            s2.signalname=['TV(\Delta) of ',obj.signal2.signalname];
            s2.graph_signal(false);
            h=get(gca,'Title');
            present_titstr{1}=get(h,'String');
            present_titstr{2}=['|||.||_1=', num2str(l1,'%3.1e')];
            title(present_titstr,'fontsize',12);
            
        end
        
        function graph_totalVarDelta(obj, open_new)
            % ~() show total variation of diff operator (not motion compensated)
            obj.requireOnly(isa(obj.signal2,'Signal2D'),'local','2d signal is set');
            if ~exist('open_new','var') || isempty(open_new)
                prepfigure(obj.fig,[],obj.figopt);
            end
            [img,l1]=obj.totalVarDelta();
            s2=TimeSignal(img);
            s2.signalname=['motion compensated TV(\Delta_\tau) of ',obj.signal2.signalname];
            s2.graph_signal(false);
            h=get(gca,'Title');
            present_titstr{1}=get(h,'String');
            present_titstr{2}=['|||.||_1=', num2str(l1,'%3.1e')];
            title(present_titstr,'fontsize',12);
            
        end
        
        
    end
    
end



