classdef CompressiveSensing_U <CSSparseSignal_U
    % CS-methods to reconstruct a sparsifiable signal from an underdetermined measurement y= A*x.
    % In contrast to class CSSparseSignal_U the sparsifying
    % transform Phi is applied implicitly within the solvers.
    %
    % Phi is the sparsifying transform: implemented as fast transform or frame
    %
    % Example:
    % 1.) define a 2-dim. signal, e.g.
    %     signal=Signal2D.make_fromImage('Mondrian.tif');
    %     signal=Signal2D.make_fromImage('cameraman.bmp');
    %     signal=Signal2D.make_zernike(48);
    %     signal=Signal2D.make_fromImage('parrots.tiff'); signal.crop();L=130;
    % 2.) define some sparsifying transform:
    %     -- matlab wavelet toolbox:       Phi=Wavelet2D_mlab(signal);
    %     -- OR wavelab toolbox: WavePath; Phi=Wavelet2D_wlab(signal);
    %     Phi.set_basisname('db1');
    %     -- OR curvelab toolbox:          Phi=Curvelet2_clab(signal);
    %
    % CASE 1: RANDOM MATRICES as encoder:
    %-------------------------------------
    % 3.) create an object of this class using a signal resized to N:
    %     cs=CompressiveSensing_U();
    %     N=32; cs.set_original(signal,N);
    %     c=2; cs.set_compression(c);
    %     --- check 99% level of sparsifying compression
    %     ct=cs.s1w_Donoho
    % 4.) set sparsifying transform (decoder):
    %     cs.set_Phi(Phi);
    %
    % 5.) set random matrix (encoder) and start simulation:
    %     SNR=[]; cs.test_Gaussian(SNR);
    %
    % CASE 2: SAMPLING ONS as encoder:
    %-------------------------------------
    % 3.) use sampling transform:
    %     a) select measurement frame:
    %     --- natural basis for sampling: B=IdTransform(signal);
    %     --- or Hadamard Transform:      B=WalshHadamard2(signal);
    %     --- or Fourier Transform:       B=Fourier2(signal);
    %     --- or Hartley Transform:       B=HartleyTrafo2(signal);
    %     --- or Wavelet level 1:         B=Wavelet2D_mlab(signal);
    %                                     B.set_deepestlev(1);
    %     b) select sampling nodes
    %     --- measurement compression rate: c=8; params=struct; params.w1=0.999; 0.95;
    %     --- random nodes:                 nodes=B.nodesPDF(c,params);
    %     --- or radial sampling:   RL=50;  nodes=B.nodesOnRadialLines(RL,c);
    % 4.) create object:                           cs=CompressiveSensing_U();
    %     set original signal:                     cs.set_original(signal);
    %     simulate measurement and reconstruction: cs.sim_SampledONS(Phi,B,nodes);
    %     --- compare with optimum (only reliable for larger images):
    %     cs.set_original(signal); cs.show_recon;
    %
    % CASE 3: 1-dimensional signal
    % ------------------------------
    %     signal=TimeSignal.make_sinus(3);
    %     Phi=HartleyTrafo1(signal);
    %     B=IdTransform(signal);
    %     c=10; params.tol=0.01; nodes=B.nodesPDF(c,params);
    %     cs=CompressiveSensing_U();
    %     cs.set_original(signal);
    %     cs.sim_SampledONS(Phi,B,nodes);
    %
    % CASE 4: 3-dimensional signal
    % ------------------------------
    %     s2=Signal2D.make_fromImage('cameraman.bmp');
    %     L=32; videoflag=2; motionunit=[0.1,0.1];
    %     signal=Signal3D.make_CyclefromLowerDimSig(s2,videoflag,L,motionunit);
    %     Phi=Wavelet3D_mlab(signal);
    %     B=HartleyTrafo3(signal);
    %     -- Samling Nodes: c=10; nodes=B.nodesDefault(c);
    %     -- OR:  RL=60;nodes=B.nodesOnRadialLines(RL);
    %     cs=CompressiveSensing_U();
    %     cs.set_original(signal);
    %     cs.sim_SampledONS(Phi,B,nodes);
    %
    % --- 12.4.14: Treating the first iteration of the solvers separately (iter==1)
    %     yields for sampled ONS (SampleTrafo_U) a first iteration
    %     with small l2-error (but it may still have a large l1-norm).
    %
    
    
    properties (SetAccess=protected)
        APhi   %<Frame> not yet implemented: modified measurement matrix (saving memory in absence of fast transforms)
        nodes  %@<SampledNodes> used to test SampleTrafo_U
        
    end
    
    
    %% constructors and commands
    methods
        
        function obj=CompressiveSensing_U(hA,hPhi)
            % constructor obj=~(hA, [hPhi]) setting encoder and decoder
            % transforms
            if nargin <1
                hA=[];
            end
            if nargin <2 || isempty(hPhi)
                hPhi=[];
            end
            obj = obj@CSSparseSignal_U(hA,hPhi);
            
            obj.removeSparsityDefects=false;
            
        end
        
        
        function set_APhi(obj,hA,hPhi)
            % ~(hA,hPhi) set encoder A and decoder Phi
            set_APhi@CSSparseSignal_U(obj,hA,hPhi);
            obj.nodes=obj.A.nodes;
        end
        
        function set_APhiNodes(obj,Phi,B,nodes,signalsize)
            % ~(Phi,B,nodes) set sampled ONS
            % with decoder Phi, encoder B sampled at nodes.
            obj.requireOnly(isa(Phi,'FrameTrafo'),'local', 'decoder is FrameTrafo');
            obj.requireOnly(isa(B,'FrameTrafo'),'local', 'encoder is FrameTrafo');
            obj.requireOnly(isa(nodes,'SampledNodes'),'local', 'nodes type ok');
            obj.requireOnly(nargin>=4 || ~obj.xorig.isemptydata,'local', 'needs size of signal');
            if nargin <5
                signalsize=size(obj.xorig);
            end
            s=obj.xorig;
            if ~isequal(s.size,signalsize)
                s=obj.xorig.make();
                s.set_signal(zeros(signalsize));
            end
            if ~isequal(s,B.ts)
                B.set_signal(s);
            end
            if ~isequal(s,Phi.ts)
                Phi.set_signal(s); % use set_signal to be sure of a reset
                obj.xopt=[];
            end
            % construct sampled decoder:
            A=SampleTrafo_U(s);
            A.set_transform(B,nodes);
            % set encoder and decoder:
            obj.set_APhi(A,Phi);
            obj.ensureOnly(~isempty(obj.nodes),'local','nodes have been set');
        end
        
        function sim_SampledONS(obj,Phi,B,nodes)
            % ~(Phi,B,nodes) simulate measurement and reconstruction for a sampled
            % orthonormal system with decoder Phi, encoder B sampled at
            % nodes.
            obj.requireOnly(isa(Phi,'FrameTrafo'),'local', 'decoder is FrameTrafo');
            obj.requireOnly(isa(B,'FrameTrafo'),'local', 'encoder is FrameTrafo');
            obj.requireOnly(isa(nodes,'SampledNodes'),'local', 'nodes type ok');
            
            % set system
            obj.set_APhiNodes(Phi,B,nodes);
            % simulate measurement vector:
            obj.set_y(obj.A.sim_measurement());
            % reconstruct:
            obj.rec();
            % show result:
            obj.show_recon;
            % SampleTrafo22.synthesize and .recSynthesis have changed B.ts:
            B.ts=obj.xorig;
        end
        
        function compute_SampledONS(obj,Phi,B,nodes,y,signalsize)
            % ~(Phi,B,nodes,y) decoder Phi, encoder B sampled at nodes, meas. y
            obj.requireOnly(isa(Phi,'FrameTrafo'),'local', 'decoder is FrameTrafo');
            obj.requireOnly(isa(B,'FrameTrafo'),'local', 'encoder is FrameTrafo');
            obj.requireOnly(isa(nodes,'SampledNodes'),'local', 'nodes type ok');
            obj.requireOnly(isequal(length(y),length(nodes)),'local', 'y is measurement vector');
            
            % set system
            obj.set_APhiNodes(Phi,B,nodes,signalsize);
            % set measurement vector:
            obj.set_y(y(:));
            % reconstruct:
            obj.rec();
            % show result:
            obj.show_recon;
        end
        
        function computeRandomAPhi(obj,m)
            % ~(m,n) multiply Phi with random matrix randn(m,n).
            % used to optimize A (minimal t-averaged coherence)
            obj.require(nargin <2 || (~isempty(obj.y) && obj.y.numel==m),...
                'local','number of measurments ok');
            
            if ~exist('m','var')
                m=obj.y.numel;
            end
            n=obj.Phi.dim;
            assert(m<=n,'CS situation: # measurements << unknowns');
            
            obj.A=[];
            obj.Aphi=obj.Phi.RandnMatTimesDualFrame(m,n);
            obj.ensure(isa(obj.APhi,'Frame'),'local', 'APhi generated');
        end
        
        
        function set_nodesByMask(obj,mask)
            %~(mask) convert mask to nodes (!=0 & ~isnan)
            obj.requireOnly(isa(obj.A,'SampleTrafo_U'),'local',...
                'needs aSampleTrafo_U as encoder');
            obj.A.set_nodes(mask(~isnan(mask) & (mask~=0)));
            
        end
        
    end
    
    
    %% solvers
    methods (Access=protected)
                
        function Solve_PD(obj, PFstar)
            % ~(PFStar) fast iterative primal-dual algo applied to basis pursuit
            % plus additional regularisation term FoA given by a proximal
            % map PFStar.
            % min_x F(Ax)+ lambda *||x||_1
            % (original x must be a sparse signal!)
            
            function s= PG(z)
                % ~(z) proximal map for G is the soft threshold map.
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
            
            if ~obj.Phi.isTrafodone
                obj.Phi.dec;
            end
            p1=obj.Phi.frame_length;
            p2=obj.A.frame_dim;
            % x1 = zeros(p1,1);
            % start values: analyze seesm to to better than
            % adjointSynthesis in choosing a good start value:
            % x1=obj.Phi.adjointSynthesis(obj.A.recSynthesis(obj.y));

            x1=obj.Phi.analyze(obj.A.recSynthesis(obj.y));
            x1=x1(:);
            xdash=x1;
            xi=zeros(p2,1);
            
            normy=norm(obj.y,1);
            iter=1;
            
            label=['primal-dual (itmax=',num2str(maxIters),') ...'];
            multiWaitbar( label, 0);

            
            while  (iter <= maxIters)
                multiWaitbar( label, iter / maxIters);               
                
                % replace A by A o Phi and A' by Phi' o A'=Psi o A':
                % reconstruct in the natural basis from sparse vector xdash:
                z=obj.Phi.synthesize(xdash);
                
                % apply wiener filter on the original data
                % z = WienerFilter(z,3);
                
                if isempty(obj.A.frame)
                    % fast op
                    xdual=obj.A.recAnalysis(z);
                    zprimal=obj.A.recSynthesis(xi);
                else
                    % matrix op
                    xdual=obj.A.frame.data*z;
                    zprimal=obj.A.frame.data'*xi;
                end
                
                % here we need the adjoint of the synthesis op. Phi:
                xprimal=reshape(obj.Phi.adjointSynthesis(zprimal),[],1);
                
                xi=PFstar(xi+sig*xdual,sig);
                xn=x1;
                x1=PG(x1- tau*xprimal);
                xdash=x1+ theta*(x1-xn);
                
                iter=iter+1;
                % print norm of the solution at every 10 iter
                if mod(iter,10)==0
                    norm(x1)
                end
            end
            multiWaitbar( label, 'Close' );
            
            
            % final x1 is best approximation to solution
            obj.x=obj.xorig.make_like();
            obj.x.replace_signal(x1);
            
            obj.solutionProps.numIters = iter-1;
            obj.solutionProps.l1normSolution=norm(x1,1);
            obj.solver_paramstr=['[\sigma,\tau]=',vec2str([sig,tau],'%3.1g'),', '...
                '\lambda=', num2str(lambda,'%3.1g')];
            
            obj.postprocess();
        end
        
    end % Solve_PD
    
    methods
        
        function SolveBPDN_FISTA(obj)
            % ~() fast iterative shrinkage-thresholding algorithm:
            % accelerated proximal gradient method applied to basis pursuit
            % denoising:
            % min_x 1/2*||Ax-y||_2^2 + tau *||x||_1
            % (original x must be a sparse signal!)
            
            obj.require(~isempty(obj.y),'local', 'measurement vector set');
            obj.requireOnly(~isempty(obj.Phi.ts),'local','signal set in Phi');
            
            
            L=obj.params.L;
            tau=obj.params.FISTA.tau;
            gamma=tau/obj.params.L;
            OptTol=obj.params.OptTol;
            obj.solutionProps=struct;
            maxIters=obj.params.FISTA.maxIters;
            
            t1=1;
            
            if ~obj.Phi.isTrafodone
                obj.Phi.dec;
            end
            p=obj.Phi.frame_length;
            zhat = zeros(p,1);
            x1hat=zhat;
            normy=norm(obj.y);
            iter=1;
            
            label= ['FISTA (itmax=',num2str(maxIters),') ...'];
            multiWaitbar( label, 0);

            
            while  (iter <= maxIters)
                
                multiWaitbar( label, iter / maxIters); 

                % Since xhat are the coefficients in the frame Phi, e.g. wavelet,
                % reconstruct x in the natural basis from sparse vector xhat:
                z=obj.Phi.synthesize(zhat);
                
                if isempty(obj.A.frame)
                    % fast op
                    res=obj.A.recAnalysis(z)-obj.y;
                    zrec=obj.A.recSynthesis(res/L);
                else
                    % matrix op
                    res=obj.A.frame.data*z(:)-obj.y;
                    zrec=obj.A.frame.data'*res(:)/L;
                end
                
                % transform residuum backwards to Phi-frame	where it is sparse
                zhatdiff=obj.Phi.analyze(z(:)-zrec(:));
                
                % replacing A by AoPhi amounts to replacing A' by Phi'oA'.
                % here we need the adjoint of the synthesis op. Phi,
                % which is the same as Psi (analyze op.) only if orthogonal
                % zhatdiff=zhat(:)-reshape(obj.Phi.adjointSynthesis(zrec),[],1);
                
                x0hat=x1hat(:);
                % compute proximal mapping p_G of 1-norm operator
                x1hat  = (abs(zhatdiff) > gamma) .* ...
                    (zhatdiff - gamma.*sign(zhatdiff));
                
                t0=t1;
                t1=(1+sqrt(4*t0^2+1))/2;
                
                lambda=1+(t0-1)/t1;
                % next iteration of zhat:
                if iter==1
                    % for sampled ONS (SampleTrafo_U)
                    % the first correction is already quite close to the
                    % original
                    zhat=zhatdiff;
                else
                    zhat=x0hat+lambda*(x1hat(:)-x0hat);
                end
                
                iter=iter+1;
                
            end
            
            multiWaitbar( label, 'Close'); 

            
            % final x1hat is best approximation to solution
            obj.x=obj.xorig.make_like();
            obj.x.replace_signal(x1hat);
            
            normresrel=norm(res)/normy;
            obj.solutionProps.numIters = iter-1;
            obj.solutionProps.normresrel=normresrel;
            obj.solutionProps.l1normSolution=norm(x1hat,1);
            obj.solver_paramstr=['\gamma=',num2str(gamma,'%3.1g'),', '...
                '\tau=', num2str(tau,'%3.1g')];
            
            % transform from Phi-frame to natural basis:
            obj.postprocess();
            
        end %SolveBPDN_FISTA               
        
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
            maxIters=obj.params.Greedy.maxIters;
            
            tol=1e-9;
            maxit=100;
            if ~obj.Phi.isTrafodone
                obj.Phi.dec;
            end
            p=obj.Phi.frame_length;
            xhat = zeros(p,1);
            obj.activeSet=[];
            
            label=['orthogonal matching pursuit (itmax=',num2str(maxiters),') ...'];
            multiWaitbar( label, 0);   

            
            while  (iter <= maxIters)  && (normresrel > OptTol)
                
                multiWaitbar( label, iter / maxIters); 

                % Since xhat are the coefficients in the frame Phi, e.g. wavelet,
                % reconstruct x in the natural basis from sparse vector xhat:
                x0=obj.Phi.synthesize(xhat);
                
                if isempty(obj.A.frame)
                    % fast op
                    zdiff=obj.y-obj.A.recAnalysis(x0);
                    xdiff=obj.A.recSynthesis(zdiff);
                else
                    % matrix op
                    zdiff=obj.y-obj.A.frame.data*x0(:);
                    xdiff = obj.A.frame.data'*zdiff;
                end
                
                % transform residuum backwards to Phi-frame	where it is sparse
                xhatdiff=obj.Phi.analyze(xdiff);
                
                % maximum norm:
                [resmax,jmax]=max(abs(xhatdiff));
                % extend active set:
                obj.activeSet=sort([obj.activeSet,jmax]);
                
                % next iteration of xhat using new active set:
                [zhat,lsqr_flag,lsqr_relres,lsqr_iter] = ...
                    lsqr(@obj.AOp,obj.y,tol,maxit,[],[],xhat(obj.activeSet));
                xhat=zeros(p,1);
                xhat(obj.activeSet)=zhat;
                normresrel=resmax/normy;
                iter=iter+1;
                
            end
            
            multiWaitbar( label, 'Close'); 

            % final xhat is best approximation to solution
            obj.x=obj.xorig.make_like();
            obj.x.replace_signal(xhat);
            
            obj.solutionProps.numIters = iter-1;
            obj.solutionProps.normresrel=normresrel;
            obj.solutionProps.l1normSolution=norm(xhat,1);
            obj.solver_paramstr=['tol=',num2str(OptTol,'%3.1g')];
            
            
            obj.postprocess();
            
        end %SolveGreedy_OMP
        
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
            obj.requireOnly(~isempty(obj.Phi.ts),'local','signal set in Phi');
            
            if isempty(obj.params.epsilon)
                epsilon=0;
            else
                epsilon=obj.params.epsilon;
            end
            maxIters=obj.params.DR.maxIters;
            gamma=obj.params.DR.gamma;
            
            normy=max(eps,norm(obj.y));
            
            iter = 1;
            obj.solutionProps=struct;
            normresrel=Inf;
            normresrel0=Inf;
            OptTol=obj.params.OptTol;
            ok=true;
            
            if (obj.params.verbose)
                l1norm=zeros(axIters,1);
                rss=l1norm;
            end
            
            if ~obj.Phi.isTrafodone
                obj.Phi.dec;
            end
            p=obj.Phi.frame_length;
            % start value:
            xhat = zeros(p,1);      % coefficients in Phi-frame
            
            % one needs pseudo inverse, backslash (QR) does not work!
            % computational expensive for large frame.data !!!
            Ainv=[];
            if ~isempty(obj.A.frame)
                Ainv=pinv(obj.A.frame.data);
            end
            
            label=['DouglasRachford (itmax=',num2str(maxIters),') ...'];
            multiWaitbar( label, 0);   
           
            
            while  (iter <= maxIters)  %&& (normresrel > OptTol)
                
                multiWaitbar( label, iter / maxIters); 

                %while  (iter <=maxIters)
                
                % Since xhat are the coefficients in the frame Phi, e.g. wavelet,
                % reconstruct x in the naturaL basis from xhat:
                x=obj.Phi.synthesize(xhat);
                
                % residuum is measured in the natural basis:
                if isempty(obj.A.frame)
                    % fast op
                    res=obj.y-obj.A.recAnalysis(x);
                else
                    % matrix op
                    res=obj.y-obj.A.frame.data*x(:);
                end
                normres=norm(res);
                % normresrel0=normresrel;
                normresrel=normres/normy;
                
                if (obj.params.verbose) && mod(iter,20)==1
                    disp(normresrel);
                end
                
                % projection in the natural basis:
                rproj = max(1 -epsilon/normres,0).*res;
                
                % transform residuum backwards to Phi-frame
                if isempty(Ainv)
                    rprojhat=obj.A.recSynthesis(rproj);
                else
                    rprojhat=Ainv* rproj(:);
                end
                % correction in Phi-frame
                corr=obj.Phi.analyze(rprojhat);
                
                % result of proximal mapping P_F of first operator F
                % xhat1=P_F(xhat)=xhat+corr:
                xhat1 = xhat + corr(:);
                
                % compute proximal mapping p_G of 1-norm operator
                % xhatdiff is its argument
                xhatdiff = 2*xhat1 - xhat;
                % apply P_G (soft thresholding ) minus corr
                % xhat = P_G(2*xhat1 - xhat)-corr:
                if iter==1
                    % for sampled ONS (SampleTrafo_U)
                    % the first correction is already quite close to the
                    % original in the l2-norm but may still have a large
                    % sparsity or l1-norm:
                    xhat=corr(:);
                elseif isreal(xhatdiff(1))
                    xhat  = (abs(xhatdiff) > gamma) .* ...
                        (xhatdiff - gamma.*sign(xhatdiff)) - corr(:);
                else
                    xhat  = (abs(xhatdiff) > gamma) .* ...
                        (xhatdiff - gamma.*angle(xhatdiff)) - corr(:); % Soft-thresholding also valid for the complex case.
                end
                
                
                if obj.params.positivity
                    xhat1(xhat1 < 0) = 0;
                end
                
                if (obj.params.verbose)
                    l1norm(iter) = norm(xhat1,1);
                    activeSet = find(abs(xhat1) > eps);
                    fprintf('Iteration %d: |I| = %d, ||x||_1 = %g\n', iter, length(activeSet), l1norm(iter));
                    rss(iter) = normres;
                    plot(rss);drawnow
                end
                
                if obj.params.fullPath
                    obj.solutionProps.sols = [obj.solutionProps.sols xhat1];
                end
                
                iter = iter+1;
                % print norm of the solution at every 10 iter
                if mod(iter,10)==0
                    norm(xhat1)
                end
            end
            
            multiWaitbar( label, 'Close'); 

            
            % final xhat1 is best approximation to solution in Phi-frame
            obj.x=obj.xorig.make_like();
            obj.x.replace_signal(xhat1);
            
            normresrel=norm(res)/normy;
            obj.solutionProps.numIters = iter-1;
            obj.solutionProps.normresrel=normresrel;
            obj.solutionProps.l1normSolution=norm(xhat1,1);
            obj.solver_paramstr=['\gamma=',num2str(gamma,'%3.1g'),...
                ', \epsilon=',num2str(obj.params.epsilon,'%3.1g')];
            
            % transform from Phi-frame to natural basis:
            obj.postprocess();
            
        end %SolveBPDN_DouglasRachford
        
        function lssolution(obj)
            % ~() least-squares estimate of xhat once its support has been
            % obtained: y=A*x= A*Phi*xhat
            %
            % select smax largest components of the best approximation
            % of the solver where smax is given by Donoho's and Tanner's weak phase
            % transition curve.
            % parent version of CSsparseSignal must be redefined
            % because of change of frame (Phi).
            obj.require(isequal(numel(obj.x.xn),obj.Phi.frame_length),'local',...
                'obj.x is a coefficient vector of the frame obj.Phi');
            
            %ct=obj.ctrafo-1;
            ct=round(obj.s1w_Donoho())-1;
            obj.check(ct>=1,'local','sparsifying compression >=1');
            smax=ceil(numel(obj.x.xn)/ct); % max. reconstructible sparsity
            x0=sort(abs(obj.x.xn),'descend');
            thresh=max([obj.params.lambdaStop, eps, x0(smax)]);
            obj.activeSet = find(abs(obj.x.xn) > thresh);
            
            % use matlab's LSQR method to solve: ||A*Phi*xhat-y||_2=min
            % where Phi is restricted to obj.activeSet.
            xhat=obj.x.xn(obj.activeSet); % start value
            tol=1e-12;
            maxit=100;
            [xtilde,obj.solutionProps.lsqr_flag,obj.solutionProps.lsqr_relres,...
                obj.solutionProps.lsqr_iter] = ...
                lsqr(@obj.AOp,obj.y,tol,maxit,[],[],xhat);
            obj.x.xn=zeros(size(x0));
            obj.x.xn(obj.activeSet)=xtilde;
            
        end
        
        function y=AOp(obj,xtilde,transp_flag)
            % AOp(xhat,'notransp')==A*Phi*x, AOp(xhat,'transp')==(A*Phi)'*x
            % where Phi is rectricetd to obj.activeSet.
            % (function with this signature needed for lsqr)
            if strcmp(transp_flag,'transp')
                % y = (A*Phi)'*x: Range(A)-> activeSet
                assert(length(xtilde)==obj.A.frame_dim,'cond1');
                x1=obj.A.analyze(xtilde);
                y=reshape(obj.Phi.analyze(x1),[],1);
                y=y(obj.activeSet);
            else  % y = A*Phi*x: activeSet -> Range(A)
                % embed xtilde in the Hilbert space of activeSet:
                assert(length(xtilde)==numel(obj.activeSet),'cond2');
                xhat=zeros(obj.Phi.frame_length,1);
                xhat(obj.activeSet)=xtilde;
                x1=obj.Phi.synthesize(xhat);
                y=reshape(obj.A.synthesize(x1(:)),[],1);
            end
        end
        
        
    end
    
    %% queries
    methods
        
        function s= baseSignal(obj)
            % s=~() base signal to be reconstructed by solvers
            s=obj.xorig;
        end
        
        function str= runMode(obj)
            % recon mode
            Phi_str=['\Phi=',obj.sparsifyingFrameName];
            str=['DEC: ',Phi_str];
        end
        
        function nod=mask2nodes(obj,mask)
            % nod=~(mask) convert mask to nodes (!=0 & ~isnan)
            nod=SampledNodes(mask(~isnan(mask) & (mask~=0)),size(mask));
        end
        
        function sv= sparsity(obj)
            % s=~() rel. 0-norm of decomposition by Phi
            % where noise has been taken into account
            obj.Phi.set_signal(obj.xorig);
            obj.Phi.dec;
            COrig=obj.Phi.C2vec;  % original coefficients in Phi-frame
            obj.Phi.estimateNoise();
            s2=TimeSignal(COrig(:));
            s2.sigEst=max(obj.Phi.sigma,1e-6);
            sv=s2.norm(0)/numel(COrig);  % relative 0-norm
        end
        
        function sv= densityReconCoeff(obj)
            % s=~() density of reconstructed coeffients in Phi-frame
            % where noise has been taken into account
            obj.require(obj.isReconDone,'local',obj.msgDoTrafo);
            obj.Phi.set_signal(obj.x);
            obj.Phi.dec;
            CRec=obj.Phi.C2vec;  % recon. coeff. in Phi-frame
            
            obj.Phi.set_signal(obj.xorig);
            obj.Phi.dec;
            COrig=obj.Phi.C2vec; % original coefficients in Phi-frame
            obj.Phi.estimateNoise();
            sig=max(1e-6,obj.Phi.sigma);
            
            COrigRelevant=abs(COrig)>3*sig;
            CRecRelevant=abs(CRec)>3*sig;
            
            delta=abs(CRec(:)-COrig(:));
            sv=sum(CRecRelevant & COrigRelevant)/sum(COrigRelevant);
            
        end
        
        
    end
    
    
    %% tests
    methods
        
        
        
    end
    
    %% graphics
    methods
        
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='encoder and decoder are set ';
            ok= isa(obj.A,'FrameTrafo') && isa(obj.Phi,'FrameTrafo');
        end
    end
    
end

