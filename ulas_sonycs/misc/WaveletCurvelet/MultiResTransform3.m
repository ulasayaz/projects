classdef MultiResTransform3 < MultiResTransform
    % abstract 3-dim. multi-resolution transformations
    %
    
    
    properties
        
    end
    
    %% commands
    methods
        
        function obj=MultiResTransform3(signal)
            % constructor
            if nargin==0 
                signal=Signal3D();
            end
            assert(isa(signal,'Signal3D'),'3d signal');
            obj = obj@MultiResTransform(signal);
            
            % do not yet check invariant for abstract class
        end
        
        
    end
    
    
    %% queries
    
    methods
        
        function L=deepestlevAnisotropic(obj)
            % deepest level in the deepest direction
            L=obj.wmaxlev-obj.wcoarsestlev;
        end
        
    end
    
    methods (Hidden)
        
        function t=detcoef_rangemax(obj,level)
            % t=~() max value of parameter type of the detail coefficients
            t=7;
        end
        
    end
    
    %% filter
    
    methods
        
    end
    
    %% transforms
    methods
                
        
    end
    
    %% statics
    methods (Static)
        
        
    end
    
    %% graphics
    
    methods
        
        
        function graph_trafo(obj, open_new, lev, z)
            % ~([open_new, lev, z) show projection to z-plane of level lev of transform C 
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin <3 || (isnumeric(lev) && lev>=1 && lev <=obj.level),...
                'local','level computed');
            obj.requireOnly(nargin <4 || obj.isnatural(z),...
                'local','z is permissible z coordinate');
            if ~exist('lev', 'var') || isempty(lev)
                lev=1;
            end
            if ~exist('z', 'var') || isempty(z)
                z=1;
            end
            if ~exist('open_new', 'var') || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            L=1+obj.detcoef_rangemax();                       
            
            P=cell(1,8);
            P{1}=obj.appcoef();
            P{1}=P{1}(:,:,z);            
            strend=size(P{1});
            
            smin=[Inf,Inf];
            for j=1:L-1
               P{j+1}= obj.detcoef(lev,j);
               P{j+1}=P{j+1}(:,:,z);
               smin=min(smin,size(P{j+1}));
            end                                               
            
            xn=zeros(smin(1),smin(2),L);
            xn(1:strend(1),1:strend(2),1)=P{1};
            for j=2:L
                xn(:,:,j)=P{j}(1:smin(1),1:smin(2));
            end
            s3=Signal3D(xn);
            
            s3.graph_signal(false,L);
            tittext=[obj.algo.name,'(',obj.basisname,'), lvl=',num2str(lev),...
                ', z=',num2str(z),obj.add_signalname];
            title(tittext,'fontsize',12);
          
        end
        
        function graph_trend(obj,open_new)
            % ~([n]) show trend signal n in open figure window
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            lvl=obj.level;
            
            A= Signal3D(obj.appcoef());
            A.signalname=[obj.ts.signalname,':  Trend'];
            A.graph_signal(false);
            
            trend_str=['Trend ',num2str(lvl),' (LLL)'];
            % norm of trend signal (relative to norm of full coefficient vector):
            Anorm=A.norm/norm(obj.C2vec);
            tittext=['||',trend_str,'||_2= ',num2tex(Anorm,'%3.1e','none')];
            title(tittext,'fontsize',obj.fontsize);
            
        end
        
        function graph_detail(obj,lvl,type,open_new)
            % ~(orient,[n])  show fluctuation at level n in open figure window           
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(lvl>=1 && lvl<=obj.level,'local',...
                'level is below computed level');
            obj.requireOnly(isnumeric(type) && type>=1 && type <=7,'local',...
                'detail type is number in 1..7');
            if nargin <4 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            D=Signal3D(obj.detcoef(lvl,type));
            
            if ~ischar(type)
                type=obj.detail_idx2str(type);
            end
            detail_str=['Detail ',num2str(lvl),' (',type,')'];
            D.signalname=[obj.ts.signalname,':  Detail ',vec2str([lvl,type])];
            
            czmax=3000;
            cz=D.nnz;
            %filter_quantile=max(0.99,(cz-czmax)/cz);
            D.graph_signal(false,4);
            % norm of detail signal (relative to norm of full coefficient vector):
            Dnorm=D.norm/norm(obj.C2vec);
            
            tittext=['||',detail_str,'||_2= ',num2tex(Dnorm,'%3.1e','none')];
            title(tittext,'fontsize',obj.fontsize);
            
        end
        
        function show_trafo(obj,lev)
            % show transform C  in open figure window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || (isnumeric(lev) && lev>0),'local', 'integer');
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            
            suptitle([obj.ts.signalname,': ', obj.transformName],14);
            
            subplot(2,4,1);
            obj.graph_trend(false);
            
            for j=1:7
                subplot(2,4,j+1);
                obj.graph_detail(lev,j,false);
            end
            
        end
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 3d';
            ok= isa(obj.ts,'Signal3D');
        end
    end
    
    
end

