classdef Frame < HMatrix
    % frame, aka dictionary
    % GT Stuttgart, 2014
    %
    % Examples:
    % F=Frame.test_operators();
    % F.isequiangular
    % M=[0,1,0.1; 1, 5,0.5];F=Frame(M); F.test_operators(F);
    % --- random frames:
    % M=randn(3,5);
    % F=Frame.test_operators(M);
    % M=randn(200,1000);
    % F=Frame.test_operators(M);
    % --- project a ONB to get a Parseval frame:
    % F=Frame(eye(3));
    % F.frame2project([0.5,0.33,1.5]);
    % F.test_operators(F);
    % --- transform to nearly Grassman frame (minimal coherence)
    % F.frame2Grass;
    % F.fig=2; F.test_operators(F);
    % --- improve coherence of random matrices:
    % M=randn(3,5);
    % M=randn(20,100);
    % F=Frame(M); F.graph_data;
    % F.frame2Grass; F.graph_data;
    % --- test frame2Grass in situations of high initial coherence
    % M=randn(10,100);
    % -- no special treatment of highly correlated pairs
    %    fails:
    % F=Frame(M);F.param_coh.cutoffcorr_mode=0; F.frame2Grass;
    % -- rotation of of highly correlated pairs, succeeds if minimal
    %    coherence (coherenceMin) is not too large (ok for N=32*32, m=N/16)
    % F=Frame(M);F.param_coh.cutoffcorr_mode=2; F.frame2Grass;
    %
    
    properties
        % data             (from HMatrix) synthesis matrix theta' of frame
        tol_deficit=1e-6   %@ (double) tolerance for decision about frame properties
        tol_recon=1e-6     %@ (double) tolerance for reconstruction
        fig=1              %@ (int) max # of open figures
        figopt
        param_coh          %@(struct)
        itmax=1000;
        verbose=true
    end
    
    properties (Hidden)
        isParsevalFlag=false %@ (boolean) is frame Parseval frame
    end
    
    %% commands and constructors
    methods
        
        function obj=Frame(hdata,u)
            % ~([v],[u]) constructor:hdata can be a sparse matrix or size,
            % u (string) are units
            if nargin==0
                % needed to allow empty calls from subclasses
                % and prevent their premature testing of invariant
                hdata=[100,100];
            end
            if ~exist('u','var')
                u='';
            end
            obj = obj@HMatrix(hdata,u);
            obj.descr='Frame';
            obj.set_param_coherence;
            obj.ensure(true,'local','verify invariant');
        end
        
        function set_data(obj,hdata)
            % ~(hdata) set analysis matrix matrix to hdata
            obj.requireOnly(ismatrix(hdata)&& size(hdata,2)>=size(hdata(1)),...
                'local', 'frame condition: frame vectors are columns of frame matrix');
            set_data@HMatrix(obj,hdata);
            obj.isParsevalFlag=[];
            obj.ensure(true,'local','verify invariant');
        end
        
        function set_param_coherence(obj)
            % parameters and defaults for frame2Grass
            obj.param_coh.itmax=1000;
            obj.param_coh.tolmin=0.01;
            obj.param_coh.tolmax=0.5;
            obj.param_coh.shrink=0.9;  % shrink factor of correlations
            obj.param_coh.t=0.9;  % fraction beyond top t% will be shrunk
            
            % vectors with correlations beyond cutoffcorr will be treated
            % specially
            obj.param_coh.cutoffcorr=0.75;
            % 0 no special measures (fastest)
            % 1 reject vectors with correlations beyond cutoffcorr
            % 2 rotate vectors ~ correlations beyond ... (best but slowest)
            obj.param_coh.cutoffcorr_mode=2;  % best results with rotation
            
        end
        
        
    end
    
    %% factories
    methods (Static)
       
        function F= make_Gaussian(m,n)
            % F=~(mn,) make a normalized Gaussian random matrix;
            % the normalization factor 1/sqrt(m)
            % enforces E(||m^(-1/2) A*x||)= ||x||^2.
           assert(m<=n,'frame condition');                      
           M=randn(m,n)/sqrt(m);
           F=Frame(M);
           F.descr='Gauss Frame';
        end
        
        function F= make_Bernoulli(m,n)
            % F=~(m,n) make a normalized Bernoulli random matrix with values 1,-1;
            % the normalization factor 1/sqrt(m)
            % enforces E(||m^(-1/2) A*x||)= ||x||^2.
           assert(m<=n,'frame condition');                      
           M=((-1).^randi(2,m,n))/sqrt(m);
           F=Frame(M);
           F.descr='Bernoulli Frame';
        end
        
        
    end                
        
    
    %% transform commands
    % change obj.data!
    methods
        
        function frame2normalized(obj)
            % normalize columns with a common factor, so that
            % the maximal column norm is <=1
            maxn=max(obj.norm_col(2));
            obj.data=obj.data/maxn;
        end
        
        function oldnorms=frame2unit(obj)
            % normalize all frame vectors to unit length
            oldnorms=obj.norm_col(2);
            for j=1:obj.length
                obj.data(:,j)=obj.data(:,j)/oldnorms(j);
            end
        end
        
        function frame2Times(obj,lambda)
            % ~(lambda) multiplies each frame vector by lambda
            obj.data=lambda*obj.data;
        end
        
        
        function frame2rotate(obj,alpha,k,l)
            % ~(alpha,[k,l]) Givens rotation by angle alpha (rad) in plane(k,l)
            obj.requireOnly(nargin <3 || (nargin==4 && k~=l && k<=obj.dim && l<=obj.dim),...
                'local','plane given by 2 different row indices');
            obj.requireOnly(obj.dim>1,...
                'local','min. dim for rotation is 2');
            if nargin <4
                k=1; l=2;
            end
            R=speye(obj.dim);
            R(k,k)=cos(alpha);
            R(k,l)=-sin(alpha);
            R(l,k)=sin(alpha);
            R(l,l)=cos(alpha);
            obj.data=R*obj.data;
            
        end
        
        function frame2project(obj,nv)
            % ~(frame2) orthogonal projection of frame to hyperplane normal to nv
            % --- projection of any frame will be frame with the same frame
            %     bounds
            % --- if frame is a ONB then its projection will be Parseval
            %
            obj.requireOnly(isvector(nv) && length(nv)==obj.dim,'local',...
                'dimension of frame equals dimension of column vector nv');
            nv=nv(:)/norm(nv,2);
            proj=@(v) v-(nv'*v)*nv;  % orthogonal projection
            % apply projection to each column of obj.data:
            obj.data=cell2mat(cellfun(proj,num2cell(obj.data,1),'UniformOutput',false));
            % compute ONB-coefficients of projection
            obj.data=orth(obj.data)'*obj.data;
            
        end
        
        function frame2parseval(obj)
            % ~() transforms the presetn frame to a Parseval frame#
            % with the frame matrix S we have M=S^(-1/2)*obj.data is
            % parseval
            obj.data= sqrtm(obj.frameOp) \ obj.data;
        end
        
        function [Res,idxremoved]=frame2Grass(obj,tol, t, shrink)
            % ~([tol,t, shrink]) transforms present frame to one with minimal
            % t-averaged coherence, which is a better indicator than the
            % coherence for the CS-reconstruction capability of th frame.
            % (in the ideal case convergence to a Grasssman matrix with coherence
            % given by obj.coherenceMin)
            % cf. M. Elad, Optimized Proj. for CS, IEEE Transactions 55(12) 2007.
            %
            % INPUT
            % obj.param_coh  ... (struct, class property)
            %
            % OUPUT:
            % Res  ... (it*3 matrix) stores evolution of optimization
            %          (1) theor. min. coherence
            %          (2) t-mean coherence (of top cross-correlations),
            %          (3) coherence=max(|cross correlations|)
            % idxremoved ... indices of removed frame vectors due to their high
            %                cross correlations the remainig vectors (cf.
            %                param_coh.cutoffcorr).
            % STOP CONDITIONS:
            % it >itmax
            % or first rel. error of coherence to min. coherence <tol
            % or second rel. error of coherence to t-mean coherence < tol2
            %
            % NOTES:
            % --- algo is non-linear, non convex: convergence is
            %     only assured if param_coh.shrink and .t are close to 1;
            %     by decreasing these 2 params, max. convegenxe speed can
            %     be increased.
            % --- convergence rate decreases with increasing obj.length/obj.dim
            %     and increasing initial coherence.
            %     In this case try to reduce param_coh.cutoffcorr e.g. to 0.6.
            % --- reducing param_coh.cutoffcorr will result in either rotating or
            %     removing nearly collinear vectors. Removed vectors can be
            %     recovered by using the index set idxremoved.
            %
            % variant to Elad algo:
            % algo handles high correlations badly, by modifying param_coh.cutoffcorr_mode
            % one can modify treatment of highly correlated vector pairs.
            %
            %
            
            
            function rotateWorstPairs(alpha)
                % alpha rotation angle max. pi -> correlation 0
                % rotate apart pairs with a too high correlation
                % compute current frame corresponding to G:
                if Gmodification_flag
                    Gmodification_flag=false;
                    [U,S,V]=svd(G);
                    % G positive, self-adjoint  => U==V
                    % also rank(G)==obj.dim
                    obj.data=S(1:N,1:N).^0.5*U(:,1:N)';
                    % transform to unit frame
                    obj.frame2unit;
                end
                % oldcoh=obj.coherence;
                for hj=1:length(rowWorst)
                    % find 2 vector of the standard basis which have
                    % max. correlation with problem vector
                    [~,idx]=sort(abs(obj.data(:,rowWorst(hj))),'descend');
                    row=rowWorst(hj);
                    col=colWorst(hj);
                    if row>col
                        % rotate one frame vector
                        % rotation matrix will have max. influence on the
                        % two max. correlated indices:
                        R=speye(obj.dim);
                        R(idx(1),idx(1))=cos(alpha);
                        R(idx(1),idx(2))=-sin(alpha);
                        R(idx(2),idx(1))=sin(alpha);
                        R(idx(2),idx(2))=cos(alpha);
                        % rotate in 2 directions:
                        v1=R*obj.data(:,col);
                        v2=R'*obj.data(:,col);
                        % check correlation with other not rotated vector:
                        corrv1=abs(dot(v1,obj.data(:,row)));
                        corrv2=abs(dot(v2,obj.data(:,row)));
                        corrG=abs(G(row,col));
                        % vold=obj.data(:,col);
                        if corrv1<corrv2 && corrv1<corrG
                            obj.data(:,col)=v1;
                        elseif corrv2<corrG
                            obj.data(:,col)=v2;
                        end                        
                        
                    end
                end
                
                % recalculate Gram matrix
                G=obj.gramOp;
            end
            
            %-------
            obj.require(true,'local','verify invariant');
            
            if ~obj.isequiangularPossible
                warning('there is no equiangular frame that large in a Hilbert space of this dimension');
            end
            
            if ~exist('tol','var') || isempty(tol)
                % max. rel. error to theoretical minimal coherence
                % interpolate between tolmin and tolmax depending on length
                % of signal:
                Li=min(100,max(10,obj.length));
                tol=interp1([10,100],[obj.param_coh.tolmin,obj.param_coh.tolmax],Li);
            end
            if ~exist('shrink','var') || isempty(shrink)
                % shrink factor of the top inner products
                shrink=obj.param_coh.shrink;
            end
            if ~exist('t','var') || isempty(t)
                % defines the top elements: all inner products larger than
                % fraction t of the total
                t=obj.param_coh.t;
            end
            
            
            modestr=['vectors with higher correlation than ',...
                num2str(obj.param_coh.cutoffcorr,'%3.2f'),' will be '];
            switch obj.param_coh.cutoffcorr_mode
                case 0
                    modestr=[modestr,' kept and treated like the others.'];
                case 1
                    modestr=[modestr,' rejected.'];
                case 2
                    modestr=[modestr,' rotated.'];
                    % case 3
                    %   modestr=[modestr,' rotated. Dito highest correlation pairs every 100 it.'];
                otherwise
                    error('unknown value of obj.param_coh.cutoffcorr_mode');
            end
            
            % second tolerance
            Li=min(100,max(10,obj.length));
            tol2=interp1([10,100],[0,0.01],Li);
            
            % transform to unit frame and save old norms
            colnorms=obj.frame2unit;
            G=obj.gramOp;
            
            mugrass=obj.coherenceMin;
            
            L=obj.length;
            N=obj.dim;
            tic;
            
            k=0;
            relerr1=1; % relative error to theoretical min. coherence
            relerr2=1; % rel. error to mean of worst quantile
            
            % number of actual top elements:
            % exclude L diagonal elements of the square Gram matrix,
            % but take at least 1
            quantile=min(L^2-L-1,round(t*(L^2-L)));  % L are on the diagonal
            
            idxremoved=[];  % removed vectors with high cross-correlation
            itmaxh= obj.param_coh.itmax;
            Res=zeros(itmaxh,3);
            kmin=min(10,itmaxh/10);
            Gmodification_flag=false;
            
            while k<=itmaxh && (k<kmin || (relerr1>=tol && relerr2>=tol2))
                
                k=k+1;
                % locate the high cross-correlations
                gg=sort(abs(G(:)));
                % determines cutoff from quantile
                cutoff=max(mugrass,gg(quantile));
                % determine top cross-correlations (off-diagonal elements)
                % (diagonal of normalized Gram contains only 1's)
                posworst=find(abs(G-diag(diag(G)))>=cutoff);
                
                % store 3 values:
                % (1) theor. min. coherence
                % (2) t-mean coherence (of top cross-correlations),
                % (3) coherence=max(|cross correlations|)
                
                coh= max(abs(obj.iifempty(G(posworst),0)));
                meancoh=mean(abs(obj.iifempty(G(posworst),0)));
                Res(k,:)=[mugrass,meancoh,coh];
                
                % relative errors of coherence w.r.t. to min. coherence
                % and t-mean-coherence
                relerr1=(coh-mugrass)/mugrass;
                relerr2=(coh-meancoh)/meancoh;
                
                if obj.verbose && (k==1 ||  mod(k,100)==1 || ...
                        k==itmaxh || relerr1<tol || relerr2<tol2 || toc>2)
                    if k==1
                        display(modestr);
                        fprintf(1,'%6s %6s %% %10s %14s  %14s \n',...
                                'it', 't', 'limit','mean-t-coh','coherence');
                    end
                    tic;
                    % output:
                    % loop index, fraction of top lements, min. coherence,
                    % t-mean coherence, coherence
                    fprintf(1,'%6i %6.2g %% %10.6f   %14.10f  %14.10f \n',...
                        [k,100*length(posworst)/(L^2-L), mugrass,meancoh, coh]);
                end
                
                if coh>obj.param_coh.cutoffcorr && obj.param_coh.cutoffcorr_mode~=0 ...
                        && mod(k,10)==1 %limit special treatment to very 10th it
                    % algorithm does not handle high cross-correlations well;
                    % there is a tendency of the number of nearly collinear pairs
                    % to grow.
                    if obj.param_coh.cutoffcorr_mode==1
                        % mode==1 is reject:
                        % simple method (used here): remove collinear pairs.
                        % alternative: rotate one of the nearly collinear vectors
                        %              of the frame, e.g. by a random Givens rotation by pi/2
                        outliers=abs(G-diag(diag(G)))>obj.param_coh.cutoffcorr +1e-6;
                        % symmetrize and locate high cross-correlations
                        outliers=outliers | outliers';
                        outliers=sum(triu(outliers),1);
                        
                        idx=find(outliers);
                        idxremoved=[idxremoved, idx];
                        G(idx,:)=[];
                        G(:,idx)=[];  % both times rows takes out only 1 of the colinear pair
                        
                        % recalculate changed size-dependent parameters:
                        L=size(G,1);
                        mugrass=sqrt((L-obj.dim)/(obj.dim*(L-1)));
                        quantile=min(L^2-L-1,round(t*(L^2-L)));  % L are on the diagonal
                        %cutoff=max(mugrass,min(0.9,gg(quantile)));
                        %posworst=find(abs(G-diag(diag(G)))>=cutoff);
                    else
                        % mode== 2 is rotate:
                        [rowWorst, colWorst]=find(abs(G-diag(diag(G)))>=obj.param_coh.cutoffcorr);
                        alpha=acos(mugrass*(0.5+0.1*rand));  % rotation angle max. pi -> correlation 0
                        rotateWorstPairs(alpha);
                    end
                                        
                else                    
                    % shrink the high cross-correlations
                    G(posworst)=G(posworst)*shrink;                    
                    Gmodification_flag=true;
                    % all tests to increase small correlations failed:
                    % 1.) grow small correlations by factor 1/shrink
                    % 2.) set all smaller correlations to threshold, e.g. .01*mugrass
                    % 3.) rotate frame vectors with small correlations
                    %     towards each other
                    
                    % the original Gram matrix had rank dim=N, if this was a frame;
                    % manipulations will generally have made new G of full rank k;
                    % reduce the rank back to N by truncating the singular
                    % values to the first N values:
                    
                    [U,S,V]=svd(G);
                    S(N+1:end,1+N:end)=0;
                    G=U*S*V';
                    
                    % Normalize the columns
                    
                    G=diag(1./sqrt(diag(G)))*G*diag(1./sqrt(diag(G)));
                    
                end
            end;
            
            Res=Res(1:k,:);
            
            % compute new resulting frame:
            [U,S,V]=svd(G);
            % G positive, self-adjoint  => U==V
            % also rank(G)==obj.dim
            obj.data=S(1:N,1:N).^0.5*U(:,1:N)';
            
            % recreate old column norms:
            colnorms(idxremoved)=[];
            if std(colnorms)>obj.tol_deficit
                for j=1:obj.length
                    obj.data(:,j)=obj.data(:,j)*colnorms(j);
                end
            end
            
        end % frame2grass
        
        
    end
    
    
    %% queries
    methods
        
        function n=length(obj)
            % n=~() number of vectors in frame
            n=obj.size(2);
        end
        
        function r=redundancy(obj)
            % redundancy of the frame
            r=obj.length/obj.dim;
        end
        
        function d=dim(obj)
            % d=~() dimension of space where frame is embedded
            d=obj.size(1);
        end
        
        function v=norm_row(obj,p)
            % v=(p) p-norm of rows
            nf=@(row) norm(nonzeros(obj.data(row,:)),p);
            v=arrayfun(nf,1:obj.size(1));
        end
        
        function v=norm_col(obj,p)
            % v=(p) p-norm of columns
            nf=@(col) norm(nonzeros(obj.data(:,col)),p);
            v=arrayfun(nf,1:obj.size(2));
        end
        
        function v=norm2(obj)
            % v=~(p) computes estimate of 2->2 operator norm of frame
            % useful for sparse or large matrices.
            v=normest(obj.data);
        end
        
        function v=norm(obj,p)
            % v=~(p) computes matrix norm p of frame
            % p=1: op.norm 1->1,
            % p=Inf: op. norm Inf->Inf
            % p='fro': Frobenius norm
            % p=2: op. norm 2->2 (largest singular value is hardest to calculate!)
            % bounds, which compute faster, are:
            % norm(.,2) <= sqrt(norm(.,1)*norm(.,Inf));
            % norm(.,2) <= norm(.,'fro');
            
            v=norm(obj.data,p);
        end
        
        function b=bounds(obj)
            % b=~() returns lower and upper frame bounds
            b(1)=eigs(obj.frameOp,1,'sm');
            b(2)=eigs(obj.frameOp,1,'lm');
        end
        
        function b=condest(obj)
            % b=~() condition number of frame matrix w.r.t. inversion   
            obj.requireOnly(obj.issquare,'local',...
                'frame is square');
            b=condest(obj.data);
        end
        
        function mumin=coherenceMin(obj)
            % mininmal mathematically possible coherence for the size of this frame
            [n,k]=size(obj.data);
            mumin=sqrt((k-n)/(n*(k-1)));
            
        end
        
        function [mu,pair]=coherence(obj, do_loop)
            % c=~([do_loop]) max. of off-diagonal Gram matrix with normalized column
            % vectors; result is sqrt((k-n)/(k*(k-1)))<=mu<=1.
            % Result is nearly zero for a ONB.
            % note: one may try flag do_loop if memory is not sufficient
            %       to compute the Gram matrix, but
            %       the run-time might be very long time
            obj.require(true,'local','verify invariant');
            if nargin <2
                do_loop = false;
            end
            
            if ~do_loop  % compute full gram matrix
                G=obj.gramOpNormalized-eye(obj.length);
                [mu,idx]=max(abs(G(:)));
                [pair(1),pair(2)]=ind2sub(size(G),idx);
            else
                % loop through Gram matrix
                [n,k]=size(obj.data);
                nv=obj.norm_col(2);
                mu=0;
                pair=[0,0];
                for j1=1:k
                    nv1=nv(j1);
                    for j2=1:k
                        if j1~=j2 % off-diagonal elements only
                            nh=nv1*nv(j2); % normalisation
                            if nh>0
                                h=abs(dot(obj.data(:,j1),obj.data(:,j2)))/nh;
                                if h>mu, mu=h; pair=[j1,j2]; end
                            end
                        end
                    end
                end
                mu=full(mu);  % if sparse
            end
            eps0=1e-9;
            obj.ensureOnly(mu<=1+eps0&& mu>=obj.coherenceMin-eps0,'local',...
                'theoretical bounds for frames');
        end
        
        function meancoh=coherence_tmean(obj,t)
            % c=~(t) t-averaged coherence, t in [0,1]
            obj.require(true,'local','verify invariant');
            if nargin <2
                t=obj.param_coh.t;
            end
            % normalize and save old column norms:
            colnorms=obj.frame2unit;
            % compute normalized Gram matrix
            G=obj.gramOp;
            
            % select top t%-cross correlations from Gram matrix:
            L=obj.length;
            quantile=min(L^2-L-1,round(t*(L^2-L)));  % L are on the diagonal
            gg=sort(abs(G(:)));
            % determines cutoff from quantile
            mugrass=obj.coherenceMin;
            cutoff=max(mugrass,min(0.9,gg(quantile)));
            % determine top cross-correlations (off-diagonal elements)
            % (diagonal of normalized Gram contains only 1's)
            posworst=find(abs(G-diag(diag(G)))>=cutoff);
            % compute t-mean coherence from top quantile:
            meancoh=mean(abs(obj.iifempty(G(posworst),0)));
            
            % recover old column norms
            for j=1:obj.length
                obj.data(:,j)=obj.data(:,j)*colnorms(j);
            end
            
            eps0=1e-9;
            obj.ensureOnly(meancoh<=1+eps0&& meancoh>=obj.coherenceMin-eps0,'local',...
                'theoretical bounds for frames');
        end
        
        function M=frameMat(obj)
            % M=~() returns the frame matrix containing the frame vectors
            % as column vectors.
            % do not confuse with frame operator frameOp (see below)
            M=obj.data;
        end
        
        function thetastar=synthMat(obj)
            % thetastar=~() returns synthesis matrix theta' of frame
            thetastar=obj.data;
        end
        
        function theta=analMat(obj)
            % theta=~() returns analysis matrix theta of frame
            theta=obj.data';
        end
        
        function S=frameOp(obj)
            % S=~() computes frame operator S=theta'*theta in matrix form;
            % do not confuse with frame matrix frameMat;
            % be careful: even if frame was sparse, result may be not sparse!
            % should be the identiy matrix for Parseval frames (wavelets and
            % curvelets)
            S=obj.data*obj.data';
        end
        
        function SS=frameCrossOp(obj,F2)
            % SS=~(F2) cross frame operator with frame F2
            % SS= theta_F'* theta_F2;
            % used to check duality or orthogonality of 2 frames
            % SS==0 if orthogonal frames, SS=1 if dual frames
            SS=obj.data*F2.data';
        end
        
        
        function G=gramOp(obj)
            % M=~() computes gram operator G=theta*theta' in matrix form
            % be careful: even if frame was sparse, result may be not sparse!
            % -- should be the identiy matrix for ONB (wavelets).
            % -- should be an orthogonal projection for Parseval frames
            % (curvelet), i.e. go^2==g, go'==go.
            G=obj.data'*obj.data;  % theta*theta'
        end
        
        function G=gramOpNormalized(obj)
            % M=~() computes normalized gram matrix, i.e.
            % cross-correlations
            fm=obj.data;
            p=2;
            % normalize columns of fm in p-norm:
            for col=1:size(fm,2)
                fm(:,col)=fm(:,col)/norm(nonzeros(fm(:,col)),p);
            end
            G=fm'*fm;
        end
        
        function D= dualFrame(obj)
            % D=~() computes canonical dual frame D= S^(-1)*theta'
            % if present frame was parseval, then the presetn frame is
            % returned (by reference, no copy)
            if islogical(obj.isParsevalFlag) && obj.isParsevalFlag
                D=obj;
            else
                D=Frame(obj.frameOp\obj.data);
            end
        end
        
        
    end
    
    
    %% frame type queries
    
    methods
        
        function ok= issquare(obj)
            ok=obj.size(1)==obj.size(2);
        end
        
        function ok=isbasis(obj)
            % ok ~() tests if synthesis matrix is a basis
            ok= obj.dim==obj.length && ...
                sprank(obj.data)==obj.dim;
        end
        
        function [ok,err]= isbiorthogonal(obj)
            %  ok=~() tests if frame is biorthogonal to its canonical dual
            %  frame (implied by frame being a basis)
            err= obj.iifempty(max(nonzeros(abs(obj.dualFrame'*obj.data-speye(obj.length)))),0);
            ok=err<obj.tol_deficit;
        end
        
        function [ok, err]= isdual2(obj,F2)
            % ok=~() tests if 2 frames are dual
            % by checking that their cross-frame product is the identity.
            obj.requireOnly(isa(F2,'Frame') && isequal(obj.size,F2.size),'local',...
                'same types');
            err=max(reshape(obj.frameCrossOp(F2)-eye(obj.dim),[],1));
            ok=err<obj.tol_deficit;
        end
        
        function [ok,err,angles]= isequiangular(obj)
            % ok=~() tests is frame matrix is equiangular
            % e.g. ONB are equiangular
            % optional: err ... numerical deviation
            %           angles  ... angles between frame vectors (rad)
            F2=Frame(obj.data);
            F2.frame2unit;  % needs to normalize all frame vectors
            G=F2.gramOp();
            % as G is self-adjoint, use only upper triangular part
            err=full(std(abs(G(triu(true(size(G)),1)))));
            ok=err<obj.tol_deficit;
            if nargout>2
                angles=acos(G(triu(true(size(G)),1)));
            end
        end
        
        function [ok, maxLength]= isequiangularPossible(obj)
            % ok= ~() is an equiangular frame of this size possible?
            maxLength=obj.dim*(obj.dim+1)/2;
            ok= obj.length<=maxLength;
        end
        
        function ok=isframe(obj)
            % ok ~() tests if synthesis matrix satisfies the frame condition, i.e.
            % if rank(theta')==#rows(theta')
            ok= obj.dim<=obj.length && ...
                sprank(obj.data)==obj.dim;
        end
        
        function [ok,err]= isgramOrthProj(obj)
            % ok=~() tests if Gram matrix is an orthogonal projection, i.e.
            % if it is self-adjoint and idempotent.
            % This is also equivalent to the frame being Parseval.
            G=obj.gramOp;
            err=zeros(2,1);
            err(1)=obj.iifempty(max(nonzeros(abs(G-G'))),0);
            err(2)=obj.iifempty(max(nonzeros(abs(G^2-G))),0);
            ok=all(err<obj.tol_deficit);
        end                
        
        function [ok,err]= isidempotent(obj)
            % ok=~() tests if A^2==A for synthesis matrix A
            obj.requireOnly(obj.issquare,'local',...
                'frame is square');
            err=obj.iifempty(max(nonzeros(abs(obj.data^2-obj.data))),0);
            ok=err<obj.tol_deficit;
        end
        
        function [ok,err]= isinvolution(obj)
            % ok=~() tests if A^2==I for synthesis matrix A
            obj.requireOnly(obj.issquare,'local',...
                'is square');
            err=obj.iifempty(max(nonzeros(abs(obj.data^2-eye(obj.dim)))),0);
            ok=err<obj.tol_deficit;
        end              
        
        function [ok,err]= isisometry(obj)
            % ok=~() tests if frame operator is an isometry
            % forall x with T:=theta=obj.data': ||x||^2==||Tx||==<x,T'*Tx>
            % is equivalent (polarisation) forall x,y:
            % <x,y>==<Tx,Ty>=<x,T'Ty> is equiv. (Riesz) T'T==Id.
            err=obj.iifempty(max(nonzeros(abs(obj.frameOp-speye(obj.dim)))),0);
            ok=err<obj.tol_deficit;
        end
        
        function [ok,maxn]=isnormalized(obj)
            % ok ~() checks if frame is normalized
            maxn=max(obj.norm_col(2));
            ok= maxn<=1+obj.tol_deficit;
        end
        
        function [ok, err]= isorthogonal2(obj,F2)
            % ok=~() tests if 2 frames are orthogonal
            % by checking that their cross-frame product is 0.
            obj.requireOnly(isa(F2,'Frame') && isequal(obj.size,F2.size),'local',...
                'same types');
            err=max(reshape(obj.frameCrossOp(F2),[],1));
            ok=err<obj.tol_deficit;
        end
        
        function [ok,err]= isparseval(obj)
            % ok=~() tests Parseval condition, i.e. if frame matrix is
            % identity
            err=obj.iifempty(max(nonzeros(abs(obj.frameOp-speye(obj.dim)))),0);
            ok=err<obj.tol_deficit;
            obj.isParsevalFlag= ok;
        end
        
        function [ok,err]= isselfadjoint(obj)
            % ok=~() tests if synthesis matrix is a self-adjoint matrix
            obj.requireOnly(obj.issquare,'local',...
                'is square');
            err=obj.iifempty(max(nonzeros(abs(obj.data-obj.data'))),0);
            ok=err<obj.tol_deficit;
        end
        
        function [ok,err,frameBound]=istight(obj)
            % [ok,err,bound]=~() tests if frame is tight with error err
            % and frame bound <<bound>> only reliable for small errors
            % for large erros use obj.bounds to compute frame bounds.
            S=obj.frameOp;
            frameBound=mean(diag(S));
            err=obj.iifempty(max(nonzeros(abs(S-...
                frameBound*speye(obj.dim)))),0);
            ok=err<obj.tol_deficit;
        end
        
        function [ok,err]=isuniform(obj)
            % ok= ~() checks if frame is uniform frame (all column vectors
            % have the same norm)
            err=std(obj.norm_col(2));
            ok= err<obj.tol_deficit;
        end
        
        function [ok,err]=isunit(obj)
            % [ok,d] ~() checks if frame is unit frame (all column vectors
            % have the norm 1); d(2)+1 is the mean norm of the frame.
            v=obj.norm_col(2);
            err(1)=std(v);
            err(2)=mean(v)-1;
            ok=all(abs(err)<obj.tol_deficit);
        end
        
        function [ok,err]= isunitary(obj)
            % ok=~() tests is frame matrix is unitary <=> frame is ONB
            % (orthonormal basis)
            [ok1,err(1)]=obj.isparseval;
            [ok2,err(2:3)]=obj.isunit;
            ok=ok1 && ok2;
        end
        
        
        
    end
    
    %% principal transform queries
    methods
        
        
        function C= analyze(obj,x)
            % C=~(signal) apply frame analysis operator to 2D image x
            % will result in frame coefficients <x,x_i> of signal x,
            % where x_i are the frame vectors (i=1:k)
            obj.requireOnly(numel(x)==obj.dim,'local','x is in frame''s Hilbert space');
            C=obj.data'*x(:);
        end
        
        function z=synthesize(obj,C)
            % xn=~(C) apply frame synthesis operator to frame coeffs C
            % z will reconstruct C if frame is Parseval, i.e.
            % if obj.frameOp == identiy
            obj.requireOnly(numel(C)==obj.length,...
                'local','C is in the frame''s range space');
            z=obj.data*C(:);
        end
        
        function z=recAnalysis(obj,C)
            % xn=~(C) reconstruct signal (vector) from frame coefficients C
            % is a left inverse of theta=obj.data'
            obj.requireOnly(numel(C)==obj.length,...
                'local','C is in the frame''s range space');
            if islogical(obj.isParsevalFlag) && obj.isParsevalFlag
                % use synthesis op.
                z=obj.reconViaParseval(C);
            else
                % we need the dual frame
                z=obj.reconViaDualFrame(C);
            end           
        end   
        
        function x=recSynthesis(obj,y)
            % xn=~(C) reconstruct signal (vector) from frame coefficients C
            % is a left inverse of theta=obj.data'
            obj.requireOnly(numel(y)==obj.dim,...
                'local','C is in the frame''s domain');
            x= obj.data \ y; % x= pinv(obj.data)* y;
                      
        end
        
        function z= reconViaDualFrame(obj,C)
            % z=~(C) reconstruct signal from frame coeffs. C using the dual frame
            obj.requireOnly(numel(C)==obj.length,...
                'local','C is in the frame''s range space');
            z=(obj.frameOp\obj.data)*C(:);
        end
        
        function z= reconViaParseval(obj,C)
            % z=~(C) reconstruct signal from frame coeffs using the synthesis op.
            % synthesis op. == dual frame if frame is Parseval (framMat =Id)
            obj.requireOnly(numel(C)==obj.length,...
                'local','C is in the frame''s range space');
            z=obj.data*C(:);
        end
        
        function z= reconViaConjGradient(obj,C)
            % z=~(C) reconstruct signal from frame coeffs. C using a conjugate
            % gradient method.
            % dual frame may have bad condition number, in this case
            % use conjugate gradient algo
            obj.requireOnly(numel(C)==obj.length,...
                'local','C is in the frame''s range space');
            % frame op S
            S= obj.data*obj.data';
            % initialisation
            uk=zeros(obj.size(1),1);
            rk= obj.data*C(:);   % r0=S(x)
            pk=rk;
            pkm=zeros(size(pk));
            
            err=1;
            it=0;
            itmin=4; % M=[0,1,0.1; 1,5,0.5]; yields erroneous recon if itmin=1
            canproceed=true;
            while canproceed && (err>obj.tol_recon && it<obj.itmax) || it<=itmin
                
                sk=S*pk;
                pksk=dot(pk,sk);
                lambdak=dot(rk,pk)/pksk;
                delk=lambdak*pk;
                % compute next it of uk, rk, pk, pkm
                uk=uk+delk;
                rk=rk-lambdak*sk;
                
                denom=dot(pkm,S*pkm);
                pk=sk-pk*norm(sk)/pksk;
                if denom~=0
                    pk=pk- pkm*dot(sk,S*pkm)/denom;
                end
                
                pkm=pk;
                it=it+1;
                err=norm(delk);
                canproceed=all(isfinite(pk));
            end % while
            z=uk;
            
            if ~canproceed
                [ok,maxn]=isnormalized(obj);
                if ~ok
                    warning('some column norms are bigger than 1. Try normalization');
                end
            end
        end
        
    end
    
    %% test functions
    methods (Static)
        
        function F=test_operators(M)
            % F=~([M]) test operations with synthesis matrix M.
            F=[];
            if nargin <1 || isempty(M)
                % M=sqrt(2/3)*[1, -1/2,  -1/2; 0,  sqrt(3)/2,  -sqrt(3)/2];
                % F=Frame(M);
                F=Frame(eye(3)); % ONB
                F.frame2project([1,1,1]); % orthogonal projection of a ONB is Parseval
            end
            if isempty(F)
                if isa(M,'Frame'), F=M; else F=Frame(M); end
            end
            
            props=[F.isframe,F.isparseval,F.isunitary];
            prepfigure(F.fig);
            
            subplot(1,3,1);
            F.graph_data(false);
            
            % testvector
            testv=ones(F.dim,1);
            tv=Signal2D();
            tv.xn=testv;
            
            if ismember(F.dim,[2,3]) && F.length<20
                hold all;
                scale=0;
                if F.dim==2
                    quiver(0,0,testv(1),testv(2),scale, 'color','red');
                else
                    quiver3(0,0,0,testv(1),testv(2),testv(3),scale, 'color','red');
                end
                legend('frame','testvector v','location','best');
            end
            
            middle_str1=',';
            middle_str2=middle_str1;
            if F.length>3
                middle_str2=',...,';
            end
            if F.dim>3
                middle_str1=',...,';
            end
            C= F.analyze(tv.xn);
            
            % coefficients of testv in standard basis and frame:
            coeff1_descr=vec2str(testv(1:min(F.dim,3)),'%2.0f');
            coeff1_descr=['(\langle v,e_i\rangle) = ', coeff1_descr(1:end-1),...
                middle_str1(1:end-1),']'];
            
            coeff2_descr=vec2str(C(1:min(F.length,3)),'%3.1f');
            coeff2_descr=['(\langle v,x_i\rangle) = ', coeff2_descr(1:end-1),...
                middle_str2(1:end-1),']'];
            
            xn=F.rec1(C);
            err=norm(tv.xn(:)-xn(:))/norm(xn(:));
            title({['Frame properties (frame,parseval,unitary)=',vec2str(props)],...
                [coeff1_descr,', ', coeff2_descr],...
                ['recon error(v) = ',num2tex(err,'%3.1e','none')]},'fontsize',12);
            
            subplot(1,3,2);
            F.graph_frameOp(false);
            subplot(1,3,3);
            F.graph_gramOp(false);
        end
        
    end
    
    
    %% graphics
    methods
        
        function graph_data(obj,open_new, maxsize)
            % show spy results in the open window (default)
            if nargin <2
                open_new=true;
            end
            if nargin <3
                maxsize=1e3;
            end
            if open_new
                prepfigure(obj.fig);
            end
            
            if ismember(obj.dim,[2,3]) && obj.length<20
                z=zeros(1,obj.length);
                scale=0;
                if obj.dim==2
                    quiver(z,z,obj.data(1,:),obj.data(2,:),scale);
                else
                    quiver3(z,z,z,obj.data(1,:),obj.data(2,:),obj.data(3,:),scale);
                end
                
                title(strtrim(obj.descr),'fontsize',12);
            else
                graph_data@HMatrix(obj,false,maxsize);
            end
            
        end
        
        function graph_frameOp(obj,open_new)
            % show spy results in the open window (default)
            if nargin <2 || open_new
                prepfigure(obj.fig);
            end
            
            imgsize=obj.dim^2;
            resize_factor=min(1,1e3/sqrt(imgsize));
            if resize_factor<1
                f=imresize(obj.frameOp,resize_factor);
            else
                f=obj.frameOp;
            end
            
            diag_dom=full(mean(reshape(abs(diag(f)),[],1))/...
                mean(abs(f(triu(true(size(f)),1)))));
            add_info=['diagonal dominance=',num2str(diag_dom,'%3.2g')];
            
            density=nnz(f)/numel(f);
            
            if diag_dom>1e3
                colormap('gray');
            else
                colormap('default');
            end
            
            if density<0.01
                spy(f);
            else
                imagesc([1,obj.dim],[1,obj.dim],f);
                if obj.length<10
                    set(gca,'XTick',1:obj.dim);
                    set(gca,'XTickLabel',1:obj.dim);
                    set(gca,'YTick',1:obj.dim);
                    set(gca,'YTickLabel',1:obj.dim);
                end
                
                colorbar; cbfreeze;
            end
            
            freezeColors;
            
            title({strtrim(obj.descr),...
                'frame matrix \theta^*\cdot\theta', add_info},...
                'fontsize',12);
            colormap('default');
        end
        
        function graph_gramOp(obj,open_new)
            % show spy results in the open window (default)
            if nargin <2 || open_new
                prepfigure(obj.fig);
            end
            
            imgsize=obj.dim^2;
            resize_factor=min(1,1e3/sqrt(imgsize));
            if resize_factor<1
                f=imresize(obj.gramOpNormalized,resize_factor);
            else
                f=obj.gramOpNormalized;
            end
            
            diag_dom=full(mean(reshape(abs(diag(f)),[],1))/...
                mean(abs(f(triu(true(size(f)),1)))));
            add_info=['diagonal dominance=',num2str(diag_dom,'%3.2g')];
            
            density=nnz(f)/numel(f);
            
            if diag_dom>1e3
                colormap('gray');
            else
                colormap('default');
            end
            
            if density<0.01
                spy(f);
            else
                imagesc([1,obj.length],[1,obj.length],f);
                if obj.length<10
                    set(gca,'XTick',1:obj.length);
                    set(gca,'XTickLabel',1:obj.length);
                    set(gca,'YTick',1:obj.length);
                    set(gca,'YTickLabel',1:obj.length);
                end
                colorbar;cbfreeze;
            end
            freezeColors;
            
            title({strtrim(obj.descr),...
                'correlation (normalized Gram \theta\cdot\theta^*)',...
                add_info},...
                'fontsize',12);
            colormap('default');
            
        end
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='frame vectors are the columns of the frame matrix.';
            ok= obj.size(2)>=obj.size(1);
        end
    end
    
end

