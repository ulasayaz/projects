classdef MotionQuotientField <MotionField
    % motion quotient field operating on frames (matrices)
    % Example
    % ========
    %{
    % --- create quotient field
    %     -- larger matrix:             
             params.M=6; params.N=6; 
    %     -- or smaller default matrix: 
             % params=struct;
        qf=MotionQuotientField.test_SimpleIllustration(params);        
        qf.show_AllAssociatedFields;
    %     -- using external method to compesnate motion:
        qf.set_method_compensate(0);
        qf.show_AllAssociatedFields;
    %
    % -- motion from videos:
        L=3; fn='tennis.avi'; %fn='riksch1.avi';
        s3=Signal3D.make_fromVideo(fn,L);  
        s3.shiftSignalValues(-1);
        qf=MotionQuotientField.make_fromSignal3(s3);
        qf.fig=5;
        qf.show_field_TV();       
        qf.show_motionCompensation();
        qf.test_motionCompensation(1);
    	qf.test_motionCompensation();
    %}
    
    
    properties (SetAccess=protected)
        domain         
    end
    
    properties (Access=protected)
         
    end
    
    properties
        
    end
    
    %% commands
    methods
        
        function obj= MotionQuotientField (qf)
            % constructor
            assert(nargin <1 || ismember(length(size(qf)),[2,3]),'qf dimensions ok');
            if nargin==0
                qf=zeros(2,2,2);
            end
            obj = obj@MotionField(qf);         
            obj.motionfname='quotient field V';     
            obj.meth_comp='intensitity';
            
            obj.ensure(true,'local','check invariant');
        end
        
        function set_motionQF(obj,qf)
            % ~(qf,[M]) set the quotient field
            obj.requireOnly(isempty(obj.signal3) || size(qf,3)==obj.signal3.size(3)-1,...
                'local','qf dimensions ok');
            obj.mat=qf;            
        end
        
        function set_domain(obj,dom)
           % ~(dom) set domain of field 
           obj.requireOnly(isempty(dom) || isequal(obj.size(1),size(dom,1)),'local',...
                'matches data size');
           obj.domain=dom; 
           obj.ensureOnly(obj.sizes_match,'local','sizes ok');
        end
                      
       
    end
    
   
    %% transformations
    methods
        
        function qf=inverseField(obj)
            % qf=~() inverse quotient field for move-backward
            qf=obj.clone;
            qf.mat=1./obj.mat;            
            qf.motionfname='inverse field 1/V';                                   
        end
                
       
        
    end
    
    %% queries
    methods
                                                                 
        function v= FieldOnDomain(obj,j)
            % v=~(j) field values on domain (rest set to 0);
            v=obj.mat(:,:,j);
            if ~isempty(obj.domain)
                v(~obj.domain(:,:,j))=0;
            end
        end
        
        function ss=size(obj,d)
            % s=~([d]) size along dimension d
            %ss=[obj.M, size(obj.mat,1)/obj.M,size(obj.mat,2)];
            ss= size(obj.mat);
            if length(ss)<3
                ss(3)=1;
            end
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
                 
        function ok=sizes_match(obj)
            % ok=~() check whether sizes of weightField and domain match
            ok=isempty(obj.domain) || (isequal(obj.size_space, [size(obj.domain,1),size(obj.domain,2)]) && ...
                isequal(obj.size_time,size(obj.domain,3)));           
        end
        
        function s=frame(obj,j)
            % s=~(j) retrieves frame number j from mat
            obj.requireOnly(min(j)>=1 && max(j)<=obj.size_time,'local',...
                'j admissible');
            st=obj.size_time;
            s=obj.clone_deep;
            s.mat=obj.mat(:,:,j);
            s.domain=obj.domain(:,:,j);
            obj.ensureOnly(~eq(obj,s),'local','new object created');
            obj.ensureOnly(obj.size_time==st && s.size_time==length(j),'local',...
                'frame length ok');
            obj.ensureOnly(obj.sizes_match,'local','sizes ok');
        end
        
        function vf2=make_like(obj)
            % vf2=~() same as clone but without matrix content
            vf2=make_like@MotionField(obj);
            vf2.domain=[];
        end     
        
        function vf2=clone(obj)
            % clone object
            obj.require(true,'local','check invariant');
            vf2=clone@MotionField(obj);
            vf2.domain=obj.domain;        
            obj.ensureOnly(~eq(obj,vf2),'local','new object created');
        end
        
        function y=moveForward(obj,j,x)
            % s2=~(j) compensate motion of j-th frame to get a
            % postdiction for the (j-1)-th frame.
            obj.requireOnly(j>=1 && j<=obj.size_time,'local',' j permitted');
            obj.requireOnly(nargin>=3 || isa(obj.signal3,'Signal3D'),'local','move argument available');
            % use only information from time step j:
            if nargin <3
                x=obj.signal3.xn(:,:,j);
            end    
            y=x.*obj.mat(:,:,j);
           
            % test values: w should be closer to xn(j-1) than v:
            % test=[obj.signal3.xn(j-1,:)',x(:),y(:)];
        end
        
        function y=moveBackward(obj,j,x)
            % s2=~(j) compensate motion of j-th frame to get a
            % prediction for the (j+1)-th frame.
            obj.requireOnly(j>1 && j<=obj.size_time+1,'local',' j permitted');
            obj.requireOnly(nargin>=3 || isa(obj.signal3,'Signal3D'),'local','move argument available');
            % use only information from time step j:
            if nargin <3
                % for backward move use signal index 1 ahead of motion field
                % index:
                x=obj.signal3.xn(:,:,j+1);
            end
            
            y=x./obj.mat(:,:,j);
            
            % test values: w should be closer to xn(j+1) than v:
            % test=[obj.signal3.xn(j+1,:)',x(:),y(:),motionV'];
        end     
        
        function y=moveForwardAdj(obj,j,x)
            % y=(j,x) adjoint of move forwards, here sell-adjoint!
            if nargin <3
                x=obj.signal3.xn(:,:,j);
            end
            y=obj.moveForward(j,x);
        end
        
        function y=moveBackwardAdj(obj,j,x)
            % y=(j,x) adjoint of move backwards, here self-adjoint!
            if nargin <3
                % for backward move use signal index 1 ahead of motion field
                % index:
                x=obj.signal3.xn(:,:,j+1);
            end
            y=obj.moveBackward(j,x);
        end
        
    end
    
    %% test
    
    methods
        
       
    end
    
    methods (Static)
                
        function qf=make()
            % contructor with universal name (callable by parent)
            qf=MotionQuotientField();
        end
        
       function qf=test_SimpleIllustration(params)
            % qf=~([params]) test used for illustration purposes
            % if params is given it needs fields A,qf,v2.
            % representing teh first video frame A (possibly a second B) and a 
            % vector field with x-component v2 and y-component qf.
            if ~exist('params','var') || isempty(params)
                params=struct;
            end
            if ~isfield(params,'rand')
                params.rand=false;
            end
            if ~isfield(params,'A')
                params.A=[0,   0.2, 0.4;
                    0.6, 0.8, 1];
                if params.rand
                    c=numel(params.A);
                    params.A= reshape(params.A(randperm(c)),size(params.A));
                end
            end
            sizeA=size(params.A);
            if ~isfield(params,'qf') || ~isequal(sizeA,size(params.qf))
                if isequal(sizeA,[2,3]) && ~params.rand                   
                    params.qf=[0.1,  0.1,  0.1;
                        2, 0.1, 1];                   
                else
                    params.qf=2*abs(randn(sizeA));                    
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
                params.qf(2,2)=2;
                mp=ceil(Mdiff/2);
                np=ceil(Ndiff/2);
                params.A=wkeep(padarray(params.A,[mp,np],1),sv);                
                params.qf=wkeep(padarray(params.qf,[mp,np],1),sv);                
            end
           
            assert(isequal(size(params.A),size(params.qf)),'same sizes');
            
            [M,N]=size(params.A);
            L=2;
            
            qf=MotionQuotientField();
            qf.signal3=Signal3D();
            qf.signal3.xn=zeros(M,N,L);
            
            qf.signal3.xn(:,:,1)=params.A;                
            qf.set_field(params.qf);            
            if isfield(params,'B') && isequal(size(params.A),params.B)
                qf.signal3.xn(:,:,2)=params.B;
            else    
                try
                    qf.signal3.xn(:,:,2)=qf.moveForward(1);
                catch
                end
            end
            assert(qf.invariant,'result satisfies invariant');
        end
        
       
    end
    
    %% more constructors
    methods (Static)
        
        function qf=make_fromSignal3(s3,opts)
            % qf=~(s3,[method_estim]) make from 3d signal using method_estim as motion estimator
            % method_estim = 'OF' - Optical Flow
            % method_estim = 'FSBM' - Full Search Block Matching
            requireOnly(isa(s3,'Signal3D'),'local','needs 3d signal');
            requireOnly(nargin <2 || isempty(method_estim) || ...
                ismember(method_estim,{'OF','FSBM'}),'local','method_estim admitted');
            
            if ~exist('opts','var') || isempty(opts)
                opts=struct;
            end
            
            qf=MotionQuotientField();
            M=s3.quotientSignal.xn; 
            qf.set_signal(s3);
            qf.set_motionQF(M);  
                        
            
            
        end               
        
        
    end
    
    %% transforms
    methods
        
       
    end
    
    %% graphics
    methods                
        
        function graph_field(obj,j, quantile, open_new)
            % ~([j]) show quotient field of all frame pairs or of (j,j+1);
            % j .... frame index (default 1)
            %
            obj.require(nargin <2 || isempty(j) || (j>=1 && j<=obj.size_time),...
                'local', 'valid index');
            if ~exist('open_new','var') || isempty(open_new)
                open_new=true;
            end
            if ~exist('j','var') || isempty(j)
                j=1;
            end
            if ~exist('quantile','var') || isempty(quantile)
                quantile=0;  % show all               
            end
            
            if open_new
                prepfigure(obj.fig,[],obj.figopt);
            end
            
            zstr=['z=',num2str(j)];
            imagesc(obj.mat(:,:,j)); colorbar;
            
            N=obj.size_space;
            if max(N)<10
                set(gca,'XTick',1:N(2));
                set(gca,'XTickLabel',1:N(2));
                set(gca,'YTick',1:N(1));
                set(gca,'YTickLabel',1:N(1));
            end
            
            axis ij; % image convention
            axis image;     
            caxis([min(obj.mat(:)), max(obj.mat(:))]);
                        
            tit=[obj.motionfname,' at ',zstr,...
                ' (method ',obj.method_estim,')'];            
            title(tit,'fontsize',12);
        end
                       
        
        function show_field(obj)
            % show quotient field of all frame pairs
            prepfigure(obj.fig,[],obj.figopt);
            
            ss=obj.size;
            cx=[min(obj.mat(:)), max(obj.mat(:))];          
            hL=ss(3);
            suptitle(obj.motionfname,14);
            sd=factor_subplots(hL);
            for j=1:hL
                subplot(sd(1),sd(2),j);
                imagesc(obj.mat(:,:,j)); colorbar;
                caxis(cx);
                axis ij; % image convention
                xlim([1,ss(2)]);
                ylim([1,ss(1)]);
                title(['z=',num2str(j)],'fontsize',12);
            end
        end
                 
        function show_AllAssociatedFields(obj)
            % ~() show field and all associated fields
             prepfigure(obj.fig,[],obj.figopt);
                    
             N=obj.size_space;
             cx=[obj.signal3.min, obj.signal3.max];
            
             subplot(2,2,1);             
             obj.graph_field(1,[],false);
             
             subplot(2,2,3);
             obj.inverseField.graph_field(1,[],false);
             xlims=xlim;
             ylims=ylim;
             
             subplot(2,2,4);
             k=1;
             s= obj.moveForwardSignal(k);   
             s.colormap_freeze=false;
             s.graph_signal(false);
             xlim(xlims); ylim(ylims);
             caxis(cx);
             if max(N)<10
                set(gca,'XTick',1:N(2));
                set(gca,'XTickLabel',1:N(2));
                set(gca,'YTick',1:N(1));
                set(gca,'YTickLabel',1:N(1));
             end
             
             subplot(2,2,2);
             k=1;
             s= obj.signal3.frame(k);   
             s.colormap_freeze=false;
             s.graph_signal(false);
             xlim(xlims); ylim(ylims);
             caxis(cx);
             if max(N)<10
                set(gca,'XTick',1:N(2));
                set(gca,'XTickLabel',1:N(2));
                set(gca,'YTick',1:N(1));
                set(gca,'YTickLabel',1:N(1));
             end
            
        end
        
        
    end
    
    
    
end



