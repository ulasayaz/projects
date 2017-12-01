classdef ConvexOptim <DC
    % solving convex optimization problems using the primal-dual algorithm
    %
    % example:
    % -----------
    % cp=ConvexOptim();
    % cp.testdataWithOutliers();
    % --- change parameters:
    %   -- start value 0:            cp.x0=zeros(size(cp.A,2),1);
    %   -- start value LSQ solution: cp.x0= cp.A \cp.y;
    %   -- number of iterations:     cp.paramsRecon.maxIters=2000;
    %   -- call solver         :     cp.SolvePrimalDualAlgo();
    %   -- plot solution:            cp.graph_linearRegression;
    %
    % restrictions:
    % as yet no stopping criterion, algo proceeds always to paramsRecon.maxIters
    % TO DO: find a good stopping criterion
    %
    
    properties (SetAccess=protected)
        A      %@<matrix> matrix A in the optimization problem min F(Ax)+G(x)
        PFstar %@<function_handle> proximal map of F* in th OP: min F(Ax)+G(x)
        PG     %@<function_handle> proximal map of G in th OP: min F(Ax)+G(x)
        y      %@<vector> y in the OP min ||Ax-y||_1
    end
    
    properties
        xdag      %@<vector> solution of primal problem min F(Ax)+G(x)
        xidag     %@<vector> solution of dual problem max -F^*(xi)-G^*(-A^*(xi))
        x0        %@<vector> start value for primal problem
        algo         %@<struct> finger print of algorithm
        paramsRecon   %@<struct> parameter list
        solutionProps %@<struct> solution properties
        optimProblem %@(string)
        
        fig       %<integer) graphical output
        figopt    %<integer> (default empty struct) figure options
        fontsize  %<integer)
    end
    
    
    %% commands
    methods
        
        function obj=ConvexOptim(hA)
            % constructor obj=~([hA]) setting linear map A in OP min F(Ax)+G(x)
            obj.requireOnly(nargin <1 || isempty(hA) || ismatrix(hA), 'local',...
                'A is matrix');
            
            if nargin <1 || isempty(hA)
                obj.A=[];
            else
                obj.A=hA;
            end
            
            obj.solutionProps=struct;
            obj.paramsRecon=struct;
            
            obj.reset_paramsRecon;
            obj.set_algo;
            
            obj.fig=1;
            obj.fontsize=12;
            obj.figopt=struct;
            obj.figopt.pixsizeX=1000;
            
        end
        
        function set_algo(obj)
            % ~() set identification data of class and general algorithm
            obj.algo.version=0.9;
            obj.algo.versiondate='30.6.2014';
            obj.algo.name='primal-dual';
            obj.algo.toolbox=[];
            obj.algo.mfilename=class(obj);
            obj.algo.author='G. Troll';
        end
        
        function reset_paramsRecon(obj)
            % ~() set default parameter values for reconstruction methods
            
            obj.paramsRecon.verbose = true;
            obj.paramsRecon.fullPath = false; % 1 returns entire solution path, 0 only final
            obj.paramsRecon.maxIters = 5000;
            
            obj.paramsRecon.theta = 1;
            obj.paramsRecon.sigma=0.01; % depends on noise sig, e.g. eps= 0.25*sig*sqrt(signal.numel);
            obj.paramsRecon.tau=0.01;  % FISTA
            
        end
        
        function set_A(obj,hA)
            % ~(A) set linear map A in the OP: min F(Ax)+G(x)
            obj.requireOnly(ismatrix(hA),'local','A is matrix');
            obj.A=hA;
            obj.initialize;
        end
        
        function set_y(obj,hy)
            % ~(y) set y in the OP min ||Ax-y||_1
            obj.requireOnly(isvector(hy),'local','y isvector');
            obj.y=hy;
            obj.initialize;
        end
        
        function set_PFstar(obj,hP)
            % ~(PFstar) set proximal map of the convex conjugate F* in the OP: min F(Ax)+G(x)
            obj.requireOnly(isa(hP,'function_handle'),'local','PFstar is function handle');
            obj.PFstar=hP;
            obj.initialize;
        end
        
        function initialize(obj)
            obj.xdag=[];
            obj.xidag=[];
            obj.x0=[];
        end
        
        function set_PG(obj,hP)
            % ~(PG) set proximal map of G in th OP: min F(Ax)+G(x)
            obj.requireOnly(isa(hP,'function_handle'),'local','PG is function handle');
            obj.PG=hP;
            obj.initialize;
        end
        
        
        function set_linearRegression(obj,hA,hy)
            % ~(A,y) set optimization problem: min ||Ax-y||_1
            obj.set_A(hA);
            obj.set_y(hy);
            obj.set_PFstar(@obj.PFstar_linearRegression);
            obj.set_PG(@obj.PG_linearRegression);
            obj.optimProblem='min_x ||Ax-y||_1';
        end
        
        
    end
    
    %% queries
    methods
        
        
        function ok=isTrafodone(obj)
            % ok=~(): is frame decomposition of signal available?
            ok=~isempty(obj.xdag);
        end
        
        function s= PFstar_linearRegression(obj,z)
            % ~(z) proximal map for F in the linear problem min ||Ax-y||_1, with
            % i.e. min (F(Ax)+G(x) with  F(x)=||x-y||_1, G(x)=0
            sig=obj.paramsRecon.sigma;
            z0=z/sig-obj.y;
            s=z - sig*(obj.y+(abs(z0) > 1/sig).*sign(z0).*(abs(z0)-1/sig));
            % case 1: z0>1/sig: s=1
            % case 2: abs(z0)<1/sig: s=z-sig*obj.y
            % case 3: z0<-1/sig: s=-1
        end
        
        function s= PG_linearRegression(obj,z)
            % ~(z) proximal map for G in the linear problem min ||Ax-y||_1,
            % i.e. min (F(Ax)+G(x) witt  F(x)=||x-y||_1, G(x)=0
            s=z;
        end
        
        function [ok,maxpossible]=check_params(obj)
            % ok=~() one needs tau*sigma*||A||_2^2<1
            % it is faster to compute the frobenius norm
            % for which one knows: ||A||_2^2 <= ||A||_F
            AF=norm(obj.A,'fro');
            ok= obj.paramsRecon.sigma*obj.paramsRecon.tau*AF^2<1;
            maxpossible=1/AF;
        end
        
        function adapt_parameters(obj)
            % ~() adapt parameters sigma and tau to A
            AF=norm(obj.A,'fro');
            ok= obj.paramsRecon.sigma*obj.paramsRecon.tau*AF^2<1;
            if ~ok
                obj.paramsRecon.sigma=1/(2*AF);
                obj.paramsRecon.tau=obj.paramsRecon.sigma;
            end
        end
        
    end
    
    methods (Static)
        
        function str=msgDoTrafo()
            str='result available (call SolvePrimalDualAlg).';
        end
        
    end
    
    
    %% solvers
    methods
        
        function SolvePrimalDualAlgo(obj)
            % ~() solve optimization problem  min F(Ax)+G(x)
            obj.requireOnly(~isempty(obj.A) && ~isempty(obj.PFstar) &&...
                ~isempty(obj.PG),'local',' all variables set');
            
            N=size(obj.A);
            if isempty(obj.x0) || length(obj.x0) ~=N(2)
                x=zeros(N(2),1); % x=obj.A'*obj.y;
            else
                x=obj.x0;
            end
            xi=zeros(N(1),1);
            xdash=x;
            theta=obj.paramsRecon.theta;
            sig=obj.paramsRecon.sigma;
            tau=obj.paramsRecon.tau;
            maxiters=obj.paramsRecon.maxIters;
            
            iter=1;
            
            if obj.paramsRecon.verbose
                wait_handle = waitbar(0,['DouglasRachford ',...
                    '(itmax=',num2str(maxIters),') ...']);
            end
            
            while  (iter <= maxiters)
                
                if obj.paramsRecon.verbose
                    waitbar(iter /maxIters);
                end
                
                xi=obj.PFstar(xi+sig*obj.A*xdash);
                xn=x;
                x=obj.PG(x- tau*obj.A'*xi);
                xdash=x+ theta*(x-xn);
                
                iter=iter+1;
                
            end
            
            if obj.paramsRecon.verbose
                close(wait_handle );
            end
            
            obj.xdag=x;
            obj.xidag=xi;
            obj.solutionProps.it=iter-1;
            
        end
        
    end
    
    %% tests
    methods
        
        function [a,b,sigy,x]=testdata_scalar(obj,a,b,sigy,x)
            % ~(a,b,sigy) testdata scattered around straight line
            if ~exist('a','var') || isempty(a)
                a=0.3;
            end
            if ~exist('b','var') || isempty(b)
                b=2;
            end
            if ~exist('sigy','var') || isempty(sigy)
                sigy=0.3;
            end
            if ~exist('x','var') || isempty(x)
                x=(1:5)';
            end
            x=x(:);
            hy=a+b*x+sigy*randn(size(x));
            hA=[ones(length(x),1),x];
            obj.set_linearRegression(hA,hy);
            
            obj.adapt_parameters();
            obj.x0=zeros(size(obj.A,2),1);
            obj.paramsRecon.maxIters=5000;
            
            obj.SolvePrimalDualAlgo();
            
            obj.graph_linearRegression(false);
            
            
        end
        
        function [a,b,sigy,x]=testdataWithOutliers(obj,a,b,sigy,x)
            % ~(a,b,sigy)  scattered around straight line with outlier
            if ~exist('a','var')
                a=0.3;
            end
            if ~exist('b','var')
                b=2;
            end
            if ~exist('sigy','var')
                sigy=0.3;
            end
            if ~exist('x','var') || isempty(x)
                x=(1:10)';
            end
            x=x(:);
            hy=a+b*x+sigy*randn(size(x));
            % create an outlier:
            hy(end)=a+b*x(end)+2*(x(end)-x(1));
            hA=[ones(length(x),1),x];
            
            obj.set_linearRegression(hA,hy);
            obj.adapt_parameters();
            obj.x0=obj.A \obj.y; % initialize with LSQ-solution
            obj.paramsRecon.maxIters=1000;
            
            obj.SolvePrimalDualAlgo();
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            
            obj.graph_linearRegression(false);
            
        end
        
    end
    
    %% graphics
    methods
        
        function graph_linearRegression(obj,open_new)
            % plot solution of linear regression
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(length(obj.xdag)==2,'local','needs 2d problem');
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            a1=obj.xdag(1);
            b1=obj.xdag(2);
            % compare l2_min
            xl2= obj.A \obj.y;
            a2=xl2(1);
            b2=xl2(2);
            
            x=obj.A(:,2);
            
            plot(x,obj.y,'o', x, a1+b1*x, x,a2+b2*x,':');
            legend('data','l_1-fit','l_2-fit','location','best');
            title([obj.algo.name,'(',obj.optimProblem,')',...
                ', it=',num2str(obj.solutionProps.it)],'fontsize',obj.fontsize);
            
        end
    end
    
end

