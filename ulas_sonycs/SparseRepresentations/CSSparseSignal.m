classdef CSSparseSignal <DC
    % CS reconstruction of a pre-sparsified signal from an underdetermined y=Ax.
    % The basis for reconstruction in this class is the property xsparse
    % which will be constructed from xorig either by sparsification in the
    % natural basis or via some sparsifying transform like the wavelet transform.
    %
    % --- Jan. 2014: numerical tests indicate that method 1 (check with get_methodActive)
    %     (Douglas-Rachford splitting) converges faster than method 2 (FISTA) and
    %     number of iterations depends only weakly on image size.
    %
    % Example
    % 1.) define a signal, e.g.
    %     signal=Signal2D.make_fromImage('Mondrian.tif');
    %     signal=Signal2D.make_fromImage('cameraman.bmp');
    %     signal=Signal2D.make_zernike(48);
    %
    % 2.) create an object of this class:
    %     cs= CSSparseSignal();
    %     N=32; cs.set_original(signal,N);
    %     c=2; cs.set_compression(c);
    %     --- check 99% level of sparsifying compression
    %     ct=cs.s1w_Donoho
    %
    % 3a) TEST1:
    %     --- sparsify given signal in the natural basis
    %         with sparsifying compression on the 99%-level curve given by s1w_Donoho:
    %     Phi=IdTransform(signal);
    %     cs.sparsifyPhiFrame(signal,Phi,N);
    %     --- choose, if final step is a LSQ restricted to the max.
    %         reconstructible active set:
    %     --- optional (final lsqr on support): cs.params.lssolution=true;
    %     --- optional (remove influence of sparsity defects): cs.removeSparsityDefects=true;
    %     --- set encoder (random distribution for the measurement matrix):
    %     SNR=[]; cs.test_Gaussian(SNR);
    %             cs.test_Bernoulli(SNR);
    %     -- increase sparsity (decrease sparsity compression) beyond
    %        critical level:
    %     ct=floor(cs.s1w_Donoho)-2; N=32; cs.sparsifyPhiFrame(signal,Phi,N,ct);
    %     SNR=[];  cs.test_Gaussian(SNR);
    %
    % 3b) TEST2:
    %    -- sparsify given signal using some sparsifying transform Phi:
    %    -- matlab toolbox: Phi=Wavelet2D_mlab(signal);
    %    -- wavelab toolbox: WavePath; Phi=Wavelet2D_wlab(signal);
    %    -- curvelab toolbox: Phi=Curvelet2_clab(signal);
    %    -- Fourier transform: Phi=Fourier2(signal);
    %    cs.removeSparsityDefects=true;
    %    N=32; cs.sparsifyPhiFrame(signal,Phi,N);
    %    --- check sparsifying capability of Phi:
    %    cs.graph_sparsityDefect();
    %    --- reconstruct with some random distribution for the measurement matrix:
    %    cs.params.lssolution=true;
    %    SNR=[]; cs.test_Gaussian(SNR);
    %    --- compare without final LSQ-step (stability vs. sparsity defects)
    %    cs.params.lssolution=false;
    %    cs.test_Gaussian(SNR);
    %
    
    properties %(SetAccess=protected)
        
        A       %@<FrameTrafo> encoder: measurement transform
        Phi     %@<FrameTrafo> decoder: sparsifying transform
        
        y       %@<vector> measured signal
        xorig   %@<SignalClass> original signal
        xsparse %@<SignalClass> sparsified signal used in this class as basis for reconstruction
        xopt    %@<SignalClass> reconstruction at quality bound (Donoho-Tanner phase transition)
        Atype   %@<string> type of measurement transform A
        c       %@<real> measurement compression rate (>=1)
        ctrafo  %@<real> sparsify compression rate (>=1)
        
        % computed properties
        x      %@<SignalClass> reconstructed signal in standard basis
        activeSet %@<vector> support of solution
        solutionProps  %@<struct> .numIters (#iterations)
        solver_paramstr %@<char> parameter output string for solver
    end
    
    properties
        % parameters of algorithm
        removeSparsityDefects %@<logical> measurement vector with/out s-defects
        params   %@<struct> parameter list
        methodsRecon %@<struct> implemented reconstruction methods
        methodActive %@<integer> active reconstruction method
        algo         %@<struct> finger print of algorithm
        
        % output parameters
        keep_reconsparse=false %@<logical> keep reconstruction in Phi-frame
        fig           %<integer) graphical output
        figopt        %<integer> (default empty struct) figure options
        fontsize      %<integer)
        
    end
    
    
    
    %% constructors and commands
    methods
        
        function obj=CSSparseSignal(hA,hPhi)
            % constructor obj=~(hA, [hPhi]) setting encoder and decoder
            % transforms
            obj.requireOnly(nargin <1 || isempty(hA) || isa(hA,'FrameTrafo'), 'local',...
                'decoder class is FrameTrafo');
            obj.requireOnly(nargin <2 || isempty(hPhi) || isa(hPhi,'FrameTrafo'), 'local',...
                'encoder class is FrameTrafo');
            
            if nargin <1 || isempty(hA)
                obj.A=[];
                obj.Atype=[];
            else
                obj.A=hA;
                obj.Atype=obj.A.frame.descr;
            end
            
            if nargin <2
                hPhi=[];
            end
            obj.Phi=hPhi;
            obj.xorig=SignalClass();
            obj.c=2;   % default measurement compression rate e.g. 2
            
            obj.reset_paramsRecon;
            obj.removeSparsityDefects=true;
            
            obj.fig=2;
            obj.fontsize=12;
            obj.figopt=struct;
            obj.figopt.pixsizeX=1000;
            
            obj.set_methodsRecon;
            % numerical tests showed that method 1
            % (Douglas-Rachford splitting) converges fastest and
            % number of iterations depends only weakly on image size
            obj.set_methodActive(1);
        end
        
        function set_algo(obj)
            % ~() set identification data of class and general algorithm
            obj.algo.version=0.9;
            obj.algo.versiondate='1.2.2014';
            obj.algo.name=[];
            obj.algo.toolbox=[];
            obj.algo.mfilename=class(obj);
            obj.algo.author='G. Troll';
        end
        
        function set_original(obj, hxorig, Nnew)
            % ~(hxorig,[Nnew])set signal to hxorig and rescale if Nnew is
            % given
            obj.requireOnly(isa(hxorig,'SignalClass'),'local',' arg is a Signal');
            if nargin <3
                Nnew=Inf;
            end
            if length(Nnew)==1
                Nnew=Nnew*ones(1,2);
            end
            obj.xorig=hxorig;
            obj.xopt=[];
            obj.xsparse=[];
            N=obj.xorig.N;
            
            if min(N)>min(Nnew)
                obj.xorig.resize(Nnew);
            end
            % do not normalize in case we need origional range of
            % intensities!
            %obj.xorig.normalize;
            % use set_signal in client classes to be sure of a reset
            if ~isempty(obj.Phi)
                obj.Phi.set_signal(obj.xorig);
            end
            if ~isempty(obj.A)
                obj.A.set_signal(obj.xorig);
            end
            obj.y= [];
        end
        
        
        function set_y(obj,hy)
            % ~(hA) set encoder (measurement vector)
            obj.requireOnly(~isempty(obj.A) ,'local', 'needs encoder A');
            obj.requireOnly(isvector(hy) && length(hy)==obj.A.frame_dim,'local',...
                'size matches with encoder A');
            obj.y=hy;
        end
        
        function set_Phi(obj,hPhi)
            % ~(hPhi) set decoder (sparsifying transform) to hPhi
            obj.requireOnly(isa(hPhi,'FrameTrafo'), 'local',...
                'decoder is a FrameTrafo');
            obj.Phi=hPhi;
            if ~isequal(obj.xorig,obj.Phi.ts)
                obj.Phi.set_signal(obj.xorig); % use set_signal to be sure of a reset
            end
            obj.x=[];  % reset reconstruction;
            obj.solutionProps=struct;
            obj.ensureOnly(~obj.isReconDone,'local', 'reconstruction reset');
        end
        
        function set_APhi(obj,hA,hPhi)
            % ~(ha,hPhi) set encoder A and decoder Phi
            obj.requireOnly(isa(hA,'FrameTrafo'),'local', 'encoder A is a FrameTrafo');
            obj.requireOnly(isa(hPhi,'FrameTrafo'),'local', 'decoder Phi is a FrameTrafo');
            
            obj.A=hA;
            if ~isequal(obj.xorig,hA.ts) && ~obj.xorig.isemptydata
                obj.A.set_signal(obj.xorig);
                obj.xopt=[];
            end
            obj.Phi=hPhi;
            if ~isequal(obj.xorig,hPhi.ts) && ~obj.xorig.isemptydata
                obj.Phi.set_signal(obj.xorig); % use set_signal to be sure of a reset
            end
            cold=obj.c;
            obj.c=obj.A.frame_redundancy;  % actual measurement compression
            if cold~=obj.c
                obj.xopt=[];
            end
            if obj.c>1
                obj.ctrafo=obj.s1w_Donoho();
            else
                obj.ctrafo=NaN;
            end
            obj.x=[];  % reset reconstruction;
            obj.solutionProps=struct;
            
            obj.ensureOnly(~isempty(obj.A) && ~isempty(obj.Phi),'local',...
                'enoder and decoder have been set');
            obj.ensureOnly(isempty(obj.A.frame) || ...
                isequal(obj.A.frame.length,obj.Phi.frame_dim),'local', 'frame dimensions match');
            obj.ensureOnly(~obj.isReconDone,'local', 'reconstruction reset');
        end
        
        function reset_A(obj)
            % ~(hA) reset encoder (encoder should only be set together with
            % decoder Phi, cf. setAPhi)
            obj.A=[];
            obj.y=[];  % reset measured signal
            obj.x=[];  % reset reconstruction;
            obj.solutionProps=struct;
            obj.ensureOnly(~obj.isReconDone,'local', 'reconstruction reset');
        end
        
        function reset_xsparse(obj)
            obj.xsparse=[];
            obj.xopt=[];
        end
        
        function set_compression(obj,c)
            % ~(c) set measurement compression c
            obj.requireOnly(c>=1,'local','compression rate >=1');
            if c~=obj.c
                obj.c=c;
                obj.ctrafo=[]; % reset sparsify compression
                obj.A=[];      % reset measurement matrix
                obj.xopt=[];   % new c implies new xopt
                obj.xsparse=[]; % new c implies new sparsified signal
            end
        end
        
        function set_Bernoulli(obj,c,SNR)
            % ~(c,SNR) simulate a measurement vector y by applying
            % a Bernoulli matrix
            % with measurement compression c and noise with a signal to
            % noise ration of SNR.
            obj.requireOnly(c>=1,'local','compression rate >=1');
            obj.requireOnly(~isempty(obj.baseSignal),'local',...
                'needs baseSignal');
            if nargin <4
                SNR=[];
            end
            obj.c=c;
            n=obj.baseSignal.numel;
            m=round(n/c);
            obj.Atype='Bernoulli';
            if isempty(obj.A) || isempty(obj.A.frame) ||  isempty(strfind(obj.A.frame.descr,'Bernou')) || ...
                    ~isequal(size(obj.A.frame),[m,n])
                frame=Frame.make_Bernoulli(m,n);
                obj.set_FrameOfA(frame,SNR);
            end
        end
        
        function set_Gaussian(obj,c,SNR)
            % ~(c,SNR) simulate a measurement vector y by applying a Gaussian
            % matrix with measurement compression c and noise with a signal to
            % noise ration of SNR.
            obj.requireOnly(c>=1,'local','compression rate >=1');
            obj.requireOnly(~isempty(obj.baseSignal),'local',...
                'needs baseSignal');
            if nargin <4
                SNR=[];
            end
            obj.c=c;
            n=obj.baseSignal.numel;
            m=round(n/c);
            obj.Atype='Gaussian';
            if isempty(obj.A) || isempty(obj.A.frame) ||  isempty(strfind(obj.A.frame.descr,'Gauss')) || ...
                    ~isequal(size(obj.A.frame),[m,n])
                frame=Frame.make_Gaussian(m,n);
                obj.set_FrameOfA(frame,SNR);
            end
        end
        
        
        function set_methodActive(obj,n)
            % ~(n) set active reconstruction method to number n
            % cf. query get_methodActive
            obj.requireOnly(n>=1 && n<= length(obj.methodsRecon),'local',...
                'method number n exists');
            obj.methodActive=n;
            obj.algo.name=obj.methodsRecon(obj.methodActive).algo_name;
            obj.x=[];  % reset reconstruction;
            obj.solutionProps=struct;
            obj.ensureOnly(~obj.isReconDone,'local', 'reconstruction reset');
        end
        
        function reset_paramsRecon(obj)
            % ~() set default parameter values for reconstruction methods
            obj.params.positivity = false; % enforcing positivity of solution
            obj.params.verbose = false;
            obj.params.fullPath = false; % 1 returns entire solution path, 0 only final
            obj.params.lssolution = false; % return the least-squares estimate once the support is obtained
            obj.params.DR.maxIters = 100;
            obj.params.PD.maxIters = 300;
            obj.params.Greedy.maxIters=100; % should equal expected sparsity
            obj.params.FISTA.maxIters = 500;
            obj.params.lambdaStop = 0; % coeffs <= lambdaStop removed for the LS solution
            obj.params.mu = 1;
            obj.params.tightFrame =0; % constant of the tight frame if known
            obj.params.DR.gamma = 1; % Relaxation parameter for Douglas-Rachford
            obj.params.OptTol=1e-2; % depends on signal
            % quadratical constraint depends on noise sig, e.g. eps= 0.25*sig*sqrt(signal.numel);
            obj.params.epsilon=0;
            
            obj.params.L=[]; % Lipschitz constant of encoder A
            obj.params.FISTA.tau=0.01;  % FISTA: weight of ||x||_1
            obj.params.PD.lambda = 0.1;  % primal-dual: weight of ||x||_1
            obj.params.PD.sigma = 0.1;
            obj.params.PD.tau = 0.1;
            obj.params.PD.theta = 1;
            
        end
        
        function set_methodsRecon(obj)
            % ~() set data of all reconstruction methods implemented
            field1 = 'function';
            value1 = {@obj.SolveBPDN_DouglasRachford ...
                ,@obj.SolveBPDN_FISTA...
                ,@obj.SolveGreedy_OMP ...
                ,@obj.SolveBPDNl1_PD ...
                ,@obj.SolveBPDN_PD ...
                ,@obj.SolveBPDNl1_PD_Cons ...
                ,@obj.SolveBPDN_PD_Cons
                };
            
            field2 = 'algo_descr';
            value2 = {'Douglas-Rachford splitting',...
                'fast iterative shrinkage-thresholding algorithm',...
                'orthogonal matching pursuit','primal-dual','primal dual',...
                'primal dual','primal dual'};
            
            field3 = 'algo_name';
            value3 = {'DRS','FISTA','OMP','PD','PD','PD','PD'};
            
            field4 = 'problemName';
            value4 = {'quadr. constrained basis pursuit',...
                'basis pursuit l2-denoising','greedy method',...
                'basis pursuit l1-denoising', 'basis pursuit l2-denoising',...
                'l1-constrained basis pursuit','quadr. constrained basis pursuit'};
            
            field5= 'problemName_short';
            value5 = {'BP_\epsilon','BPDN','greedy','BPDN_{l1}','BPDN',...
                'BPDN_{l1}_constrained','BPDN_{l2}_constrained'};
            
            field6 = 'problem';
            value6 = {'min_x ||x||_1 s.t. ||y -Ax||_2 <= epsilon',...
                'min_x \tau ||x||_1 + 1/2 |Ax-y||_2^2',...
                'max_j |(A^*(y-Ax))_j|',...
                'min_x \lambda ||x||_1+||Ax-y||_1',...
                'min_x \lambda ||x||_1 + 1/2 ||Ax-y||_2^2',...
                'min_x \lambda ||x||_1 s.t ||Ax-y||_1 <=epsilon',...
                'min_x \lambda ||x||_1 s.t. ||Ax-y||_2 <= epsilon'};
            
            field7 = 'best_usage';
            value7 = {'high dimensions, large sparsity',...
                'high dimensions, large sparsity',...
                'small sparsity','l_1 constraint',...
                'high dimensions, large sparsity',...
                'high dimensions, large sparsity',...
                'high dimensions, large sparsity'};
            
            obj.methodsRecon=struct(field1,value1,field2,value2,...
                field3,value3,field4,value4,field5,value5,field6,value6,...
                field7, value7);
        end
        
        function fPhi=compute_QualityBound(obj)
            % ~() computes obj.xopt at compression ctrafo, optional output
            % is a copy of the decoder with signal xopt and filtered
            % coefficients.
            % e.g. cs=CompressiveSensing(); cs.set_original(signal);
            %      cs.set_APhiNodes(Phi,B,nodes); fPhi=cs.compute_QualityBound;
            obj.requireOnly(~obj.xorig.isemptydata,'local','needs signal');
            if isempty(obj.ctrafo)
                obj.ctrafo=obj.s1w_Donoho();
            end
            obj.Phi.set_signal(obj.xorig);
            obj.Phi.dec;
            if nargout==0
                obj.xopt=obj.Phi.sparseApprox(obj.ctrafo);
            else
                [obj.xopt,~,~,~,fCoeff]=obj.Phi.sparseApprox(obj.ctrafo);
                fPhi=obj.Phi.clone;
                fPhi.set_signal(obj.xopt);
                fPhi.C=fCoeff;
            end
            obj.xopt.signalname='Quality Bound';
            obj.ensureOnly(isequal(obj.xorig.N,obj.xopt.N),'local',...
                'size matches that of xorig');
        end
        
        function compute_reconQuality(obj)
            % ~() computes PSNR, SSIM of reconstruction
            obj.require(obj.isReconDone,'local',obj.msgDoTrafo);
            obj.requireOnly(~obj.xorig.isemptydata,'local','needs original signal');
            obj.solutionProps.quality.PSNR=obj.xorig.PSNR(obj.x);
            obj.solutionProps.quality.SSIM=obj.xorig.SSIM(obj.x);
        end
        
        function rec(obj, heps)
            % ~([eps]) reconstruct signal using method obj.methodActive
            obj.requireOnly(nargin<2 || heps>=0,'local','epsilon >=0');
            if nargin <2
                heps=[];
            end
            
            % ---- parameters for Solvers
            redundancy=obj.Phi.frame_redundancy;
            obj.params.mu = 0.99/redundancy;
            [~,smax]=s1w_Donoho(obj);
            obj.params.Greedy.maxIters=round(1.5*smax);
            
            if ~isempty(heps)
                % reset parameter for BP_epsilon
                obj.params.epsilon=heps;
            end
            
            % parameter for BPDN
            % L must be smaller than the square of the largest singular
            % value of obj.A.fram.data;
            % LipConst=(svds(obj.A.frame.data,1))^2;  % slow
            % but we know: ||A||_F^2/n<= ||A||_{2->2}^2 <= ||A||_F^2;
            LipConst=1;
            if ~isempty(obj.A.frame)
                LipConst=norm(obj.A.frame.data,'fro')^2;
            end
            obj.params.L=LipConst;
            
            % execute active resonstruction method
            recfun=obj.methodsRecon(obj.methodActive).function;
            recfun();
            
            %obj.x.normalize;  % pv 1, mean 0 like xorig
            
            % compute reconstruction quality
            if ~obj.xorig.isemptydata
                obj.compute_reconQuality();
            end
            
            obj.ensure(obj.isReconDone,'local','reconstruction available');
        end
        
        
        function frame2Grass(obj)
            % ~() transforms present frame to one with minimal
            % t-averaged coherence, which is a better indicator than the
            % coherence for the CS-reconstruction capability of th frame.
            obj.require(~isempty(obj.A) && ~isempty(obj.A.frame),'local',...
                'needs frame');
            obj.A.frame.param_coh.itmax=30;
            obj.A.frame.param_coh.t=0.8;
            obj.A.frame.param_coh.shrink=0.7;
            obj.A.frame.frame2Grass;
            
        end
        
        
    end
    
    %% solvers
    
    methods (Access=protected)
        
        function s= PFstar_l2(obj,z,sig)
            % ~(z) proximal map for l2-denoising;
            % used in Solve_PD.
            s=(z-sig*obj.y)/(1+sig);
        end
        
        function s= PFstar_l1(obj,z,sig)
            % ~(z) proximal map for l1-denoising;
            % used in Solve_PD
            z0=z/sig-obj.y;
            s=z - sig*(obj.y+(abs(z0) > 1/sig).*sign(z0).*(abs(z0)-1/sig));
            % case 1: z0>1/sig: s=1
            % case 2: abs(z0)<1/sig: s=z-sig*obj.y
            % case 3: z0<-1/sig: s=-1
        end
        
        function s= PFstar_l1_Cons(obj,z,sig)
            % ~(z) proximal map for l1-constrained;
            % used in Solve_PD.
            s=(z-sig*obj.y)/(1+sig);
        end
        
        function s= PFstar_l2_Cons(obj,z,sig)
            % ~(z) proximal map for l2-constrained;
            % used in Solve_PD.
            s=(z-sig*obj.y)/(1+sig);
        end
        
        function Solve_PD(obj, PFstar)
            % ~(PFStar) fast iterative primal-dual algo applied to basis pursuit
            % plus additional regularisation term FoA given by a proximal
            % map PFStar.
            % min_x F(Ax)+ lambda *||x||_1
            % (original x must be a sparse signal!)
            
            function s= PG(z)
                % ~(z) proximal map for G is the soft threshold map:
                r=lambda*tau;
                s=(abs(z) > r).*sign(z).*(abs(z)-r);
            end
            
            obj.require(~isempty(obj.y),'local', 'measurement vector set');
            obj.require(isa(PFstar,'function_handle'),'local', 'proximal map is function');
            
            OptTol=obj.params.OptTol;
            obj.solutionProps=struct;
            
            theta=obj.params.PD.theta;
            sig=obj.params.PD.sigma;
            tau=obj.params.PD.tau;
            lambda=obj.params.PD.lambda;
            maxIters=obj.params.PD.maxIters;
            
            
            p1=obj.A.frame_length;
            p2=obj.A.frame_dim;
            x1 = zeros(p1,1);
            xdash=x1;
            xi=zeros(p2,1);
            
            normy=norm(obj.y,1);
            iter=1;
            
            while  (iter <= maxIters)
                
                if isempty(obj.A.frame)
                    % fast op
                    xdual=obj.A.recAnalysis(xdash);
                    xprimal=obj.A.recSynthesis(xi);
                else
                    % matrix op
                    xdual=obj.A.frame.data*xdash;
                    xprimal=obj.A.frame.data'*xi;
                end
                
                xi=PFstar(xi+sig*xdual,sig);
                xn=x1;
                x1=PG(x1- tau*xprimal);
                xdash=x1+ theta*(x1-xn);
                
                iter=iter+1;
                
            end
            
            % final x1 is best approximation to solution
            obj.x=obj.xorig.make_like();
            obj.x.replace_signal(x1);
            
            %obj.solutionProps.normresrel=norm(obj.A.recAnalysis(x1)-obj.y,1)/normy;
            obj.solutionProps.numIters = iter-1;
            obj.solver_paramstr=['[\sigma,\tau]=',vec2str([sig,tau],'%3.1g'),', '...
                '\lambda=', num2str(lambda,'%3.1g')];
            
            obj.postprocess();
            
        end %SolveBPDN_PD
        
    end
    
    methods
        
        function SolveBPDN_FISTA(obj)
            % ~() fast iterative shrinkage-thresholding algorithm:
            % accelerated proximal gradient method applied to basis pursuit
            % l2-denoising:
            % min_x 1/2*||Ax-y||_2^2 + tau *||x||_1
            % (original x must be a sparse signal!)
            
            obj.require(~isempty(obj.y),'local', 'measurement vector set');
            L=obj.params.L;
            gamma=obj.params.FISTA.tau/obj.params.L;
            OptTol=obj.params.OptTol;
            obj.solutionProps=struct;
            
            t1=1;
            
            if ~obj.Phi.isTrafodone
                obj.Phi.dec;
            end
            p=obj.Phi.frame_length;
            z = zeros(p,1);
            
            x1=z;
            normy=norm(obj.y);
            iter=1;
            
            while  (iter <= obj.params.FISTA.maxIters)
                
                if isempty(obj.A.frame)
                    % fast op
                    zdiff=obj.A.recAnalysis(z)-obj.y;
                    xdiff=z-obj.A.recSynthesis(zdiff/L);
                else
                    % matrix op
                    zdiff=obj.A.frame.data*z-obj.y;
                    xdiff = z-obj.A.frame.data'*zdiff/L;
                end
                
                x0=x1;
                x1  = (abs(xdiff) > gamma) .* ...
                    (xdiff - gamma.*sign(xdiff));
                
                t0=t1;
                t1=(1+sqrt(4*t0^2+1))/2;
                
                lambda=1+(t0-1)/t1;
                
                z=x0+lambda*(x1-x0);
                
                iter=iter+1;
                
            end
            
            
            % final x1 is best approximation to solution
            obj.x=obj.xorig.make_like();
            obj.x.replace_signal(x1);
            
            %obj.solutionProps.normresrel=norm(zdiff)/normy;
            obj.solutionProps.numIters = iter-1;
            obj.solver_paramstr=['\gamma=',num2str(gamma,'%3.1g'),', '...
                '\tau=', num2str(obj.params.FISTA.tau,'%3.1g')];
            
            obj.postprocess();
            
        end %SolveBPDN_FISTA
        
        function SolveBPDNl1_PD(obj)
            % ~() fast iterative primal-dual algo applied to basis pursuit
            % l_1-denoising:
            % min_x ||Ax-y||_1 + lambda *||x||_1
            % (original x must be a sparse signal!)
            
            obj.Solve_PD(@obj.PFstar_l1);
            
        end %SolveBPDNl1_PD
        
        function SolveBPDN_PD(obj)
            % ~() fast iterative primal-dual algo applied to basis pursuit
            % l_1-denoising:
            % min_x ||Ax-y||_1 + lambda *||x||_1
            % (original x must be a sparse signal!)
            
            obj.Solve_PD(@obj.PFstar_l2);
            
        end %SolveBPDNl1_PD
        
        function SolveBPDNl1_PD_Cons(obj)
            % ~() fast iterative primal-dual algo applied to basis pursuit
            % l_1-constrained:
            % min_x lambda *||x||_1 s.t. ||Ax-y||_1 <= epsilon
            % (original x must be a sparse signal!)
            
            obj.Solve_PD(@obj.PFstar_l1_Cons);
            
        end %SolveBPDNl1_PD_Cons
        
        function SolveBPDN_PD_Cons(obj)
            % ~() fast iterative primal-dual algo applied to basis pursuit
            % l_2-constrained:
            % min_x lambda *||x||_1 s.t. ||Ax-y||_2 <= epsilon
            % (original x must be a sparse signal!)
            
            obj.Solve_PD(@obj.PFstar_l2_Cons);
            
        end %SolveBPDNl1_PD_Cons
        
        
        function SolveBPDN_DouglasRachford(obj)
            % ~() Iterative Douglas-Rachford Splitting Proximal iteration
            % to solve the BPDN problem (basis pursuit quadratically constrained):
            % (original x must be a sparse signal!)
            %		min || x ||_1 s.t. || y -Ax ||_2 <= epsilon
            %
            %	The solution is given by the iteration :
            %		x^(t+1/2) = prox_iCoA(x^(t)) = x^(t) + pinv(A)*(P_C - Id)A x^(t)
            %		x^(t+1)   = x^(t) + mu_t (ST[2 prox_iCoA(x^(t)) - x^(t)] - x^(t+1/2)),
            %	where
            %		P_C: is the projector onto the l_2-ball centered at y and of radius epsilon,
            %       such that:
            %			(P_C - Id)(z) = (1 - epsilon/||y - z||_2)* (y - z)
            %
            %		0 < mu_t < 2, fixed to mu_t=1 in this code (optimal convergence speed).
            %		ST[x] is soft-thresholding with threshold gamma for gamma strictly positive.
            %
            %	The solution to the minimization problem is given by prox_iCoA(x^(t)) at convergence.
            %
            % Parameters:
            %	    epsilon	    desired residual error
            %   	gamma	    a strictly positive real, default 1
            %   	lambdaStop  If specified (and > 0), the algorithm removes all coeffs <= lambdaStop in magnitude for the LS solution
            %       maxIters
            % Outputs
            %	 obj.x       solution of BPDN problem
            %    obj.solutionProps  solution/iteration properties
            %
            obj.require(~isempty(obj.y),'local', 'measurement vector set');
            
            
            if isempty(obj.params.epsilon)
                epsilon=0;
            else
                epsilon=obj.params.epsilon;
            end
            gamma=obj.params.DR.gamma;
            
            normy=norm(obj.y);
            
            iter = 1;
            obj.solutionProps=struct;
            normresrel=Inf;
            normresrel0=Inf;
            OptTol=obj.params.OptTol;
            ok=true;
            
            
            if (obj.params.verbose)
                l1norm=zeros(obj.params.DR.maxIters,1);
                rss=l1norm;
            end
            
            % one needs pseudo inverse, backslash (QR) does not work!
            % computational expensive for large frame.data !!!
            if ~isempty(obj.A.frame)
                Ainv=pinv(obj.A.frame.data);
            end
            
            if ~obj.Phi.isTrafodone
                obj.Phi.dec;
            end
            p=obj.Phi.frame_length;
            x0 = zeros(p,1);
            
            while  (iter <= obj.params.DR.maxIters)  && (normresrel > OptTol)
                
                
                %while  (iter <= obj.params.DR.maxIters)
                
                % start with a sparse signal x0
                
                % residuum:
                if isempty(obj.A.frame)
                    % fast op
                    res=obj.y-obj.A.recAnalysis(x0);
                else
                    % matrix op
                    res=obj.y-obj.A.frame.data*x0;
                end
                normres=norm(res);
                % normresrel0=normresrel;
                normresrel=normres/normy;
                
                if (obj.params.verbose) && mod(iter,20)==1
                    disp(normresrel);
                end
                
                rproj = max(1 -epsilon/normres,0).*res;
                
                % transform residuum backwards
                if isempty(obj.A.frame)
                    % fast op
                    corr=obj.A.recSynthesis(rproj);
                else
                    % matrix op
                    % ---> corr=obj.A.frame.data \ rproj; DOES NOT WORK,
                    % one needs MP-pseudoinverse:
                    corr=Ainv* rproj;
                end
                % x1 is proximal mapping P_F of first operator F
                % x1=P_F(x0)=x0+corr:
                x1 = x0 + corr(:);
                
                % compute proximal mapping p_G of 1-norm operator
                % xdiff is its argument
                xdiff = 2*x1 - x0;
                % apply P_G (soft thresholding ) minus corr
                % x0 = P_G(2*x1 - x)-corr:
                if isreal(xdiff(1))
                    x0  = (abs(xdiff) > gamma) .* ...
                        (xdiff - gamma.*sign(xdiff)) - corr(:);
                else
                    x0  = (abs(xdiff) > gamma) .* ...
                        (xdiff - gamma.*angle(xdiff)) - corr(:); % Soft-thresholding also valid for the complex case.
                end
                
                if obj.params.positivity
                    x1(x1 < 0) = 0;
                end
                
                if (obj.params.verbose)
                    l1norm(iter) = norm(x1,1);
                    activeSet = find(abs(x1) > eps);
                    fprintf('Iteration %d: |I| = %d, ||x||_1 = %g\n', iter, length(activeSet), l1norm(iter));
                    rss(iter) = normres;
                    plot(rss);drawnow
                end
                
                if obj.params.fullPath
                    obj.solutionProps.sols = [obj.solutionProps.sols x1];
                end
                
                iter = iter+1;
            end
            
            % final x1 is best approximation to solution
            obj.x=obj.xorig.make_like();
            obj.x.replace_signal(x1);
            
            obj.solutionProps.numIters = iter-1;
            %obj.solutionProps.normresrel=norm(res)/normy;
            obj.solver_paramstr=['\gamma=',num2str(gamma,'%3.1g'),...
                ', \epsilon=',num2str(obj.params.epsilon,'%3.1g')];
            
            obj.postprocess();
        end
        
        function SolveGreedy_OMP(obj)
            % greedy method orthogonal matching pursuit
            % suitable for small sparsities s because
            % runtime is proportional to s.
            
            obj.require(~isempty(obj.y),'local', 'measurement vector set');
            
            normy=norm(obj.y);
            
            iter = 1;
            obj.solutionProps=struct;
            normresrel=Inf;
            
            OptTol=obj.params.OptTol/normy;
            
            
            tol=1e-9;
            maxit=100;
            if ~obj.Phi.isTrafodone
                obj.Phi.dec;
            end
            p=obj.Phi.frame_length;
            x0 = zeros(p,1);
            obj.activeSet=[];
            
            while  (iter <= obj.params.Greedy.maxIters)  && (normresrel > OptTol)
                
                if isempty(obj.A.frame)
                    % fast op
                    zdiff=obj.y-obj.A.recAnalysis(x0);
                    xdiff=obj.A.recSynthesis(zdiff);
                else
                    % matrix op
                    zdiff=obj.y-obj.A.frame.data*x0;
                    xdiff = obj.A.frame.data'*zdiff;
                end
                
                % maximum norm:
                [resmax,jmax]=max(abs(xdiff));
                % extend active set:
                obj.activeSet=sort([obj.activeSet,jmax]);
                
                % next iteration of x0using new active set:
                if isempty(obj.A.frame)
                    [z0,lsqr_flag,lsqr_relres,lsqr_iter] = ...
                        lsqr(@obj.AOp,obj.y,tol,maxit,[],[],x0(obj.activeSet));
                else
                    z0=obj.A.frame.data(:,obj.activeSet) \obj.y;
                end
                x0=zeros(p,1);
                x0(obj.activeSet)=z0;
                
                normresrel=resmax/normy;
                iter=iter+1;
                
            end
            
            % final x0 is best approximation to solution
            obj.x=obj.xorig.make_like();
            obj.x.replace_signal(x0);
            
            obj.solutionProps.numIters = iter-1;
            %obj.solutionProps.normresrel=normresrel;
            obj.solver_paramstr=['tol=',num2str(OptTol,'%3.1g')];
            
            obj.postprocess();
            
        end
        
        
        function lssolution(obj)
            % ~() least-squares estimate of obj.x once its support has been obtained;
            % select smax largest components of the best approximation
            % of the solver where smax is given by Donoho's and Tanner's weak phase
            % transition curve.
            obj.require(~isempty(obj.x),'local','solution exists');
            ct=max(1,obj.ctrafo-1);
            % ct=round(obj.s1w_Donoho())-1;
            smax=ceil(numel(obj.x.xn)/ct); % max. reconstructible sparsity
            x0=sort(abs(obj.x.xn),'descend');
            thresh=max([obj.params.lambdaStop, eps, x0(smax)]);
            activeSet = find(abs(obj.x.xn) > thresh);
            if ~isempty(obj.A.frame)
                AI = obj.A.frame.data(:,activeSet);
                x0=zeros(size(obj.x.xn));
                x0(activeSet) = AI \ obj.y(:);
                obj.x.xn=x0;
            else
                warning('lssolution for fast op. not yet implemented!');
            end
            
        end
        
        
        function postprocess(obj)
            % ~() post prcessing of reconstruction
            
            if obj.params.lssolution
                obj.lssolution();
            end
            
            % shorten signalname if possible
            sn=[];
            if ~obj.xorig.isemptydata && ~isempty(obj.xorig.signalname)
                sn=regexp(obj.xorig.signalname,'[a-zA-Z0-9]*','once','match');
            end
            obj.x.signalname=[obj.methodsRecon(obj.methodActive).problemName_short,...
                ' (',obj.algo_name,') : ',...
                sn];
            
            % reconstructed signal in sparsifying frame:
            % obj.x.graph_signal();
            
            if obj.keep_reconsparse
                obj.xsparse=obj.xorig.make_like();
                obj.xsparse.replace_signal(obj.x.xn);
                obj.xsparse.signalname='reconstruction in \Phi-frame';
            end
            
            % transform back to standard base
            obj.x.xn= obj.Phi.synthesize(obj.x.xn);
            obj.solutionProps.normresrel=norm(obj.A.recAnalysis(obj.x.xn)-obj.y)/norm(obj.y);
            
            if ~obj.xorig.isemptydata
                obj.x.xn=reshape(obj.x.xn,obj.xorig.N);
            end
            
        end
        
        function x2C(obj)
            % ~() apply decoder Phi to reconstructed signal x
            % can be shown as obj.Phi.show_trafo
            obj.Phi.set_signal(obj.x);
            obj.Phi.dec;
        end
        
        function y1=AOp(obj,xtilde,transp_flag)
            % AOp(xhat,'notransp')==A*Phi*x, AOp(xhat,'transp')==(A*Phi)'*x
            % where Phi is rectricetd to obj.activeSet.
            % (function with this signature needed for lsqr)
            if strcmp(transp_flag,'transp')
                % y = A'*x: Range(A)-> activeSet
                assert(length(xtilde)==obj.A.frame_dim,'cond1');
                y1=reshape(obj.A.analyze(xtilde),[],1);
                y1=y1(obj.activeSet);
            else  % y = A*x: activeSet -> Range(A)
                assert(length(xtilde)==numel(obj.activeSet),'cond2');
                x1=zeros(obj.A.frame_length,1);
                x1(obj.activeSet)=xtilde;
                y1=reshape(obj.A.synthesize(x1(:)),[],1);
            end
        end
        
    end
    
    %% queries
    methods
        
        function mr= get_methodActive(obj)
            % mr= ~() get struct describing active reconstruction method
            % cf. command set_methodActive
            mr=obj.methodsRecon(obj.methodActive);
        end
        
        function get_methodsAllAvailable(obj)
            % mr= ~() get struct describing all reconstruction method
            % cf. command set_methodActive
            for j=1:length(obj.methodsRecon)
                mr{j}=([num2str(j),': ',obj.methodsRecon(j).algo_name,...
                    ' (',obj.methodsRecon(j).algo_descr,...
                    '), ', obj.methodsRecon(j).problemName]);
                disp(mr{j});
            end
        end
        
        function str= sparsifyingFrameName(obj)
            % str=~() frame used for sparsification obj.Phi
            str=obj.Phi.transformName();
        end
        
        function str= runMode(obj)
            % str=~() test mode
            if obj.removeSparsityDefects
                str='testmode: presparsified, sparsity defects removed';
            else
                str='testmode: presparsified, sparsity defects kept';
            end
        end
        
        function s= baseSignal(obj)
            % s=~() base signal to be reconstructed by solvers
            s=obj.xsparse;
        end
        
        function ok= isReconDone(obj)
            % ok= ~() has reconstruction been done?
            ok=~isempty(obj.x);
        end
        
        function sd=sparsitydefect(obj,x,p)
            % sd=~(x,s,p) s-sparsity defect of x in p-norm
            % inf{||x-z||; z \in C^N s-sparse}=(\sum_{i>s} |x_i^*|^^p)^{1/p}
            obj.requireOnly(nargin <2 || isnumeric(x),'local','x is array or matrix');
            obj.requireOnly(nargin <3 || isscalar(p),'local','p is scalar');
            if nargin <2 || isempty(x)
                x=obj.xorig.xn;
            end
            if nargin <3
                p=2;
            end
            % apply sparsifying transform
            x=obj.Phi.analyze(x);
            x=sort(abs(x(:)),'ascend').^p;
            sd=flipud(cumsum(x).^(1/p));
        end
        
        
        function [ct,smax]=s1w_Donoho(obj)
            % [ct,smax]=~() max. sparsity for m measurements
            % with measurement compression obj.c (ct=smax/m, c=N/m)
            % using polynomial approx. of  a selection of image sizes Nvals
            obj.requireOnly(obj.c>=1,'local','compression rate >=1');
            obj.requireOnly(~obj.xorig.isemptydata||~isempty(obj.A),'local',...
                'needs dimensions of original signal');
            
            Nvals=[200,1000,5000,1e6];    % actually Inf instead of 1e6
            pcoeff=zeros(10,4);
            pcoeff(:,1)= [48.191991880107516,-1.823956417716720e+02,2.679606511524920e+02,...
                -1.727421918514478e+02,13.033410741844145,51.860683881826920,...
                -34.706421025995270,10.479392456699042,-1.319063895299251,0.063416449794653]';
            pcoeff(:,2)=[74.116329065705240,-2.733104687814432e+02,3.991792422889943e+02,...
                -2.810750689914758e+02,81.769797544656020,9.377110766259726,...
                -12.503667036541442,2.835226671930660,0.277397543891433,-0.006743101205354];
            pcoeff(:,3)=[1.655113261968571e+02,-6.944487339221886e+02,1.231895346815504e+03,...
                -1.205380874217993e+03,7.147525023691115e+02,-2.675212589582132e+02,...
                64.992325825612440,-10.659937056676181,1.666651237682664,-0.009342748260871]';
            pcoeff(:,4)=[2.872347691605911e+02,-1.241998859413596e+03,2.273306408354441e+03,...
                -2.293831878572575e+03,1.394582507509451e+03,-5.255212581657510e+02,...
                1.225199100299579e+02,-17.365584292201174,1.920355292125700,0.088252206152749]';
            
            delta=1/obj.c;
            if ~obj.xorig.isemptydata
                M=obj.xorig.numel;
            else
                M=obj.A.frame_dim;
            end
            
            m=delta*M;
            
            if M>Nvals(end)
                rho=polyval(pcoeff(:,end),delta);
            else % interpolate logarithmically
                s1w=@(i) polyval(pcoeff(:,i),delta);
                rhovals= arrayfun(s1w,1:4);
                rho=interp1(log(Nvals),rhovals,log(M));
            end
            smax=rho*m;
            ct=M/smax;
        end
        
        function [ct, smax]=s1w_Donoho_Asymptotics(obj)
            % [ct,smax]=~() asymptotics of max. sparsity for m measurements
            % with measurement compression c. (ct=smax/m, c=N/m)
            % weak phase transition for Gaussians according to Donoho et al.
            % in the case of M-> infty, m/M -> 0, i.e. c:=M/m -> infty
            % recovery with prob. p close to 1 if s < smax.
            % failure with high prob. if s>smax
            % c_m=M/m; c_p=M/sp;
            % source: [MICS] notes p.305
            obj.requireOnly(~isempty(obj.A),'local', 'frame set');
            M=obj.A.frame_length;
            m=obj.A.frame_dim; % # measurements
            smax=min(m,ceil(0.5*m./log(M./m)));
            ct=M/smax;
        end
        
        function [ct,smax]=s1w_RauhutEq924(obj)
            % [ct,smax]=~() max. sparsity for m measurements with
            % measurement compression c. (ct=smax/m, c=N/m)
            % in the case of large M: m=2*s_m*ln(e*M/s_m) and nonuniform
            % recovery with prob. p close to 1
            % c_m=M/m; c_p=M/sp;
            % source: [MICS] eq. 9.24
            obj.requireOnly(~isempty(obj.A),'local', 'frame set');
            M=obj.A.frame_length;
            m=obj.A.frame.dim; % # measurements
            e=exp(1);
            fun=@(s) m-2*s.*log(e*M./s);
            % fun assumes min at M; sparsity limit smax is in [1, argmin]
            argmin=M;
            fmin=fun(argmin);
            % [argmin,fmin] = fminbnd(fun,1,M);
            if fmin<0
                smax=floor(fzero(fun,[1,argmin]));
            else % no restriction as to sparsity: all vectors can be reconstructed
                smax=M;
            end
            ct=M/smax;
        end
        
        function an=algo_name(obj)
            % an=~() name of active algorithm
            an=obj.methodsRecon(obj.methodActive).algo_name;
            if obj.params.lssolution
                an=[an,'-','LSQ'];
            end
        end
        
    end
    
    
    %% tests
    methods
        
        function sparsifyPhiFrame(obj,signal,hPhi,Nnew, ct)
            % ~(signal,[Phi, Nnew,ct]) sparsify <signal> using transform Phi
            % by setting all but the Phi-transformed coefficients largest in modulus to
            % zero.
            % (c =N/m measurement compression, ct=N/smax sparsifying
            % compression)
            obj.requireOnly(isa(hPhi,'FrameTrafo'),'local',...
                'if given, must be FrameTrafo');
            obj.requireOnly(obj.c>=1,'local','compression rate >=1');
            obj.requireOnly(nargin <5 || ct>=1,'local','sparsify compression rate >=1');
            
            N=signal.N;
            if N(1)~=N(2)
                signal.crop(min(N));
            end
            if ~exist('Nnew','var')
                Nnew=min(N,256);
            end
            if ~exist('hPhi','var') || isempty(hPhi);
                if license('test','Wavelet_Toolbox')
                    hPhi=Wavelet2D_mlab();
                else
                    hPhi=Wavelet2D_wlab();
                end
                hPhi.set_basisname('db1');
            end
            
            if nargin <5 || isempty(ct)
                % set sparsify compression rate to critical one:
                ct=obj.s1w_Donoho();
            end
            
            obj.ctrafo=ct;
            
            if length(Nnew)==1
                Nnew=Nnew*ones(1,2);
            end
            
            obj.xorig=signal.clone;
            
            if min(N)>min(Nnew)
                obj.xorig.resize(Nnew);
            end
            
            if abs(1-obj.xorig.pv)> 1e-6 || abs(obj.xorig.mean)>1e-6
                obj.xorig.normalize;
            end
            
            % create sparse signal
            obj.Phi=hPhi;
            obj.Phi.set_signal(obj.xorig);
            obj.Phi.dec;
            
            % using representation of xorig in Phi basis as sparse signal:
            obj.xsparse=obj.xorig.make_like();
            obj.xsparse.replace_signal(reshape(obj.Phi.C2vec,obj.xorig.N));
            obj.xsparse.signalname=obj.sparsifyingFrameName;
            % figure, imagesc(obj.Phi.C2mat()); colorbar;
            % truncate sparse signal at critical compression rate:
            
            if obj.removeSparsityDefects
                % remove sparsity defects from xsparse
                % observe that the measurement sim will be based on xsparse!
                obj.xsparse.sparsify(obj.ctrafo);
            end
            
            % Phi2=obj.Phi.clone; Phi2.C=obj.xsparse.xn; Phi2.S=obj.Phi.S;
            % figure, imagesc(Phi2.C2mat()-obj.Phi.C2mat()); colorbar;
            
            
        end %sparsifyPhiFrame
        
        
        function test_Bernoulli(obj,SNR,force_newframe)
            % ~([SNR]) simulate a measurement vector y by applying
            % a Bernoulli matrix
            % with measurement compression obj.c and noise with a signal to
            % noise ration of SNR.
            obj.requireOnly(obj.c>=1,'local','compression rate >=1');
            obj.requireOnly(~isempty(obj.baseSignal),'local',...
                'needs baseSignal ');
            if nargin <4
                SNR=[];
            end
            if nargin <5 || isempty(force_newframe)
                force_newframe=false;
            end
            n=obj.baseSignal.numel;
            m=round(n/obj.c);
            obj.Atype='Bernoulli';
            if isempty(obj.A) || isempty(obj.A.frame) ||  isempty(strfind(obj.A.frame.descr,'Bernou')) || ...
                    ~isequal(size(obj.A.frame),[m,n]) || force_newframe
                obj.set_Bernoulli(obj.c,SNR);
            end
            % simulate measurement vector using xsparse
            obj.A.ts=obj.baseSignal; % make sure of that before sim_measurement
            
            obj.Phi.ts=obj.baseSignal;
            
            obj.y=obj.A.sim_measurement;
            obj.test_Frame();
        end
        
        function test_Gaussian(obj,SNR,force_newframe)
            % ~([SNR]) simulate a measurement vector y by applying a Gaussian
            % matrix with measurement compression obj.c and noise with a signal to
            % noise ration of SNR.
            obj.requireOnly(obj.c>=1,'local','compression rate >=1');
            obj.requireOnly(~isempty(obj.baseSignal),'local',...
                'needs baseSignal');
            
            if nargin <4
                SNR=[];
            end
            if nargin <5 || isempty(force_newframe)
                force_newframe=false;
            end
            n=obj.baseSignal.numel;
            m=round(n/obj.c);
            if isempty(obj.A) || isempty(obj.A.frame) ||  isempty(strfind(obj.A.frame.descr,'Gauss')) || ...
                    ~isequal(size(obj.A.frame),[m,n]) || force_newframe
                obj.set_Gaussian(obj.c,SNR);
            end
            % simulate measurement vector using baseSignal
            obj.A.ts=obj.baseSignal;  % make sure of that before sim_measurement
            
            obj.Phi.ts=obj.baseSignal;
            
            obj.y=obj.A.sim_measurement;
            obj.test_Frame();
        end
        
        
        function set_FrameOfA(obj,frame,SNR)
            %~(frame,[SNR] set measurement matrix to <frame>
            % and setting the SNR of the frame to <SNR>
            obj.requireOnly(~isempty(obj.baseSignal),'local',...
                'needs baseSignal');
            obj.requireOnly(isa(frame,'Frame'),'local','needs a frame');
            
            if nargin <3
                SNR=[];
            end
            
            % create encoder A
            obj.reset_A();
            obj.A=FrameTrafo();
            
            obj.A.frame=frame;
            obj.A.SNR=SNR;
            obj.A.ts=obj.baseSignal;
            
        end
        
        function test_Frame(obj)
            % ~() reconstruct signal from a simulated measurement y
            % and show result graphically.
            obj.require(true,'local','check invariant');
            
            m=obj.A.frame.dim;
            epsilon = sqrt(m)*obj.A.sigma;
            
            % reconstruct signal from measurement matrix
            obj.rec(epsilon);
            obj.compute_QualityBound;
            
            obj.show_recon;
            
        end
        
    end
    
    %% graphics
    methods
        
        function str=msgDoTrafo(obj)
            % str=~() message string
            str='reconstruction available (call rec).';
        end
        
        function graph_sparsityDefect(obj,xtest,open_new)
            % ~([xtest]) plot sparsity defect of xtest(=obj.xorig by default)
            if nargin <3 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            if nargin <2 || isempty(xtest)
                xtest=obj.xorig;
            end
            if isa(xtest,'SignalClass')
                xtest=xtest.xn;
            end
            p=2;
            sd1=obj.sparsitydefect(xtest,p);
            sd2=flipud(cumsum(sort(abs(xtest(:)),'ascend')).^(1/p));
            plot(1:numel(sd1),sd1,...
                1:numel(xtest),sd2);
            sparsified_str='sparsified';
            
            sparsified_str=[sparsified_str,' (',obj.sparsifyingFrameName,')'];
            
            legend(sparsified_str,'original','location','best');
            xlabel('sparsity s');
            ylabel(['\sigma_s(x)_{',num2str(p),'}']);
            title(['l_{',num2str(p),'}-error of best s-term approx.'],...
                'fontsize',12);
            
        end
        
        function graph_recon(obj, open_new)
            % ~() show reconstructed image (in natural basis)
            obj.requireOnly(obj.isReconDone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
                        
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            param_str=[' (it=',num2str(obj.solutionProps.numIters)];
            param_str=[param_str, ', ', obj.solver_paramstr,')'];
            
            sparsity_control=obj.runMode();
            
            result_str=[];
            opts.statfunc1=@max; opts.statfunc2=@(x) max(x)-min(x);
            if ~obj.xorig.isemptydata
                if ~isfield(obj.solutionProps,'quality')
                    obj.compute_reconQuality();
                end
                SSIM_str=[];
                if ~isempty(obj.solutionProps.quality.SSIM)
                    opts.format1='%3.2f'; opts.format2='%3.0e'; 
                    SSIM_str=[', SSIM=',stat2str(obj.solutionProps.quality.SSIM,opts)];
                end
                opts.format1='%3.1f'; opts.format2='%3.1f'; 
                result_str=['PSNR=',stat2str(obj.solutionProps.quality.PSNR,opts),...
                    SSIM_str];
            end
            
            obj.x.colormap_active=obj.xorig.colormap_active;
            
            obj.x.graph_signal(false);
            
            if length(sparsity_control)<=result_str
                result_str=[sparsity_control,'-',result_str];
                sparsity_control=[];
            end
            titstr={[obj.x.signalname, param_str],sparsity_control,result_str};
            titstr=titstr(~cellfun(@isempty,titstr))  ;
            title(titstr,'fontsize',obj.fontsize);
            
        end
        
        function graph_opt(obj, open_new)
            % ~() show quality bound (Donoho and Tanner)
            obj.requireOnly(~isempty(obj.xopt),'local',...
                'quality bound xopt computed (use compute_QualityBound)');
            obj.requireOnly(~obj.xorig.isemptydata,'local','original image exists');
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            SSIM_opt=NaN; PSNR_opt=NaN;
            if isequal(obj.xopt.size,obj.xorig.size)
                SSIM_opt=obj.xorig.SSIM(obj.xopt);
                PSNR_opt=obj.xorig.PSNR(obj.xopt);
            end
            
            Phi_str=['\Phi=',obj.sparsifyingFrameName];
            SSIM_str=[];
            opts.statfunc1=@max; opts.statfunc2=@(x) max(x)-min(x);
            if ~isempty(SSIM_opt)
                SSIM_str=[', SSIM=',num2str(mean(SSIM_opt),'%3.2f')];
            end
            commentstr=[];
            if obj.xopt.numel<1000
                commentstr=' (unreliable small sample)';
            end
            opts.format1='%3.1f'; opts.format2='%3.1f'; 
            resultOpt_str=[Phi_str,', PSNR=',stat2str(PSNR_opt,opts),...
                SSIM_str,commentstr];
            cstr=' at {\bf uniform} compression ';
            %             if abs(obj.ctrafo-obj.s1w_Donoho())<eps
            %                 cstr=' at 99%-DT-compression ';
            %             else
            %                 cstr=' at compression ';
            %             end
            param2=[cstr,vec2str([obj.c,obj.ctrafo],'%3.1f')];
            obj.xopt.colormap_active='gray';
            obj.xopt.colormap_freeze=true;
            
            obj.xopt.graph_signal(false);
            title({[obj.xopt.signalname,param2],resultOpt_str},'fontsize',obj.fontsize);
            
        end
        
        function RS=graph_reconInFrame(obj, open_new)
            % ~() show reconstructed signal coefficients in the Phi-frame
            obj.requireOnly(obj.isReconDone,'local',...
                'sparse signal computed');
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            RS=obj.xorig.make();
            obj.Phi.ts=obj.x;
            obj.Phi.dec;
            RS.set_signal(obj.Phi.C2graph());
            RS.repfun=obj.Phi.repfun;
            RS.signalname=[obj.xorig.signalname,' reconstructed signal in frame ',obj.sparsifyingFrameName];
            RS.graph_signal(false);
            
        end
        
        function RS=graph_optInFrame(obj, open_new)
            % ~() show quality-bound signal coefficients in the Phi-frame
            obj.requireOnly(obj.isReconDone,'local',...
                'sparse signal computed');
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2
                open_new=true;
            end
            if open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            RS=obj.xorig.make();
            obj.Phi.ts=obj.xopt;
            obj.Phi.dec;
            RS.set_signal(obj.Phi.C2graph());
            RS.repfun=obj.Phi.repfun;
            RS.signalname=[obj.xorig.signalname,' quality bound signal coeffs. in frame ',obj.sparsifyingFrameName];
            RS.graph_signal(false);
            
        end
        
        function graph_diff(obj, basesignal, open_new)
            % ~([basesignal]) show diffence image between reconstruction and basesignal
            obj.requireOnly(nargin <2 || isa(basesignal,'SignalClass')...
                || (isempty(basesignal) && (~isempty(obj.xopt)||~obj.xorig.isemptydata)),...
                'local','needs basesignal or original or quality bound');
            obj.requireOnly(nargin<3 || islogical(open_new),'local', 'boolean');
            if ~exist('basesignal','var') || isempty(basesignal)
                if  ~isempty(obj.xopt)
                    basesignal=obj.xopt;
                else
                    basesignal=obj.xorig;
                end
            end
            
            xdiff=basesignal.diff(obj.x);
            titstr=[basesignal.signalname, ' - reconstruction'];
            
            if nargin <3 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            
            xdiff.colormap_active='default';
            xdiff.colormap_freeze=false;  % false to have access to data
            xdiff.graph_signal(false);
            title({titstr,xdiff.get_signalname},'fontsize',12);
        end
        
        function graph_encoder(obj, open_new)
            % ~() show encoder A (measurement matrix)
            obj.require(nargin<2 || islogical(open_new),'local', 'boolean');
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            obj.A.graph_distribution(false,obj.y);
            h=get(gca,'Title');
            present_titstr=get(h,'String');
            if iscellstr(present_titstr)
                present_titstr{1}=['Encoder: ',present_titstr{1}];
            else
                present_titstr=['Encoder: ',present_titstr];
            end
            title(present_titstr,'fontsize',12);
        end
        
        function graph_encoder_diff(obj, open_new)
            % ~() difference of actual measurements to those simulated with reconstructed signal
            obj.require(nargin<2 || islogical(open_new),'local', 'boolean');
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            s1=obj.A.measurement2graph(obj.y);
            sig_encoder=obj.A.ts;
            % set signal to reconstructed signal
            obj.A.set_signal(obj.x);
            yn=obj.A.sim_measurement;
            ssim=obj.A.measurement2graph(yn);
            % reset signal
            obj.A.set_signal(sig_encoder);
            
            if isa(s1,'TimeSignal')
                s1.marker='.';
                s1.graph_signal(false);
                hold all;
                obj.x.marker='';
                obj.x.graph_signal(false);
                hold off;
                legend('measurements','reconstruction','location','best');
                title('Error graph','fontsize',12);
            else
                diff_meas=s1.diff(ssim);
                diff_meas.graph_signal(false);
            end
            
            
        end
        
        function graph_decoder(obj, open_new)
            % ~() show decoder Phi
            obj.require(nargin<2 || islogical(open_new),'local', 'boolean');
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            obj.Phi.graph_trafo(false);
            h=get(gca,'Title');
            present_titstr=get(h,'String');
            if iscellstr(present_titstr)
                present_titstr{1}=['Decoder: ',present_titstr{1}];
            else
                present_titstr=['Decoder: ',present_titstr];
            end
            title(present_titstr,'fontsize',12);
        end
        
        function graph_original(obj, open_new)
            % ~() show encoder (measurement matrix)
            obj.require(nargin<2 || islogical(open_new),'local', 'boolean');
            obj.requireOnly(~obj.xorig.isemptydata,'local','original image exists');
            
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            obj.xorig.graph_signal(false);
            title(['Original Signal:',obj.xorig.signalname],'fontsize',12);
        end
        
        function show_recon(obj)
            % ~() show reconstruction and compare with original
            obj.require(obj.isReconDone,'local',obj.msgDoTrafo);
            
            % 5 optional subplot graphics
            % ============================
            % reconstruction, encoder (distribution of measurement matrix),
            % quality bound , original, difference signal
            %
            if obj.fig==0
                return;
            end
            spcount=3;
            if ~isempty(obj.xopt)
                spcount=spcount+1;
            end
            if ~obj.xorig.isemptydata
                spcount=spcount+1;
            end
            if ~obj.xorig.isemptydata || ~isempty(obj.xopt)
                spcount=min(4,spcount+1);
            end
            
            sd=factor_subplots(spcount);
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            
            subplot(sd(1),sd(2),1);
            obj.graph_recon(false);
            
            subplot(sd(1),sd(2),2);
            obj.graph_encoder(false);
            
            variant_IdTransform=isa(obj.Phi,'IdTransform');
            
            if ~obj.xorig.isemptydata && isempty(obj.xopt)
                obj.compute_QualityBound;
            end
            
            if ~isempty(obj.xopt) && obj.xopt.numel>1000
                subplot(sd(1),sd(2),3);
                obj.graph_opt(false);
            elseif ~obj.xorig.isemptydata
                % decoder
                subplot(sd(1),sd(2),3);
                obj.graph_decoder(false);
            else
                subplot(sd(1),sd(2),3);
                obj.graph_encoder_diff(false);
            end
            
            if ~obj.xorig.isemptydata
                subplot(sd(1),sd(2),4);
                if variant_IdTransform
                    obj.graph_diff([],false);
                else
                    obj.graph_original(false);
                end
            end
            
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='encoder frame is set ';
            ok= isa(obj.A,'FrameTrafo');
        end
    end
    
end

