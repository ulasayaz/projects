classdef FrameTrafo <DC
    % n-dimensional frame transform (signal classes exist for n=1..3)
    %
    % example:
    %{
        signal=Signal2D.make_fromImage('cameraman.bmp');
        B=WalshHadamard2(signal);
        B.dec;
        B.graph_trafo;
    %
    % -- test transform and frame:
        res=B.test_framefeatures(); display(res);
        [ok,err1,err2]= f.test_DecRec(); disp(ok);
        [ok,err1,err2]= f.test_AnalyzeSynthsize(); disp(ok);
    %
    % --- UML diagram:
        exclude={'FourierTrafos2'};
        helpuml('FrameTrafo', 0, 4, exclude);
    
    %}
    
    properties (Access=public)
        
        ts      %@<SignalClass> n-dimensional signal (image wrapper)
        
        repfun  %@<function_handle> representation function for transformed signal
        use_repfun %@(boolean> use representation function in graphs of transformed signal
        energy_compressionlevel %@ double in [0,1], e.g. 0.999
        SNR     %@<real> simulated measurement noise in decibel, e.g. 30 dB
        
        tol     %@<real> numerical tolerance
        fig     %@<integer> graphical output
        figopt  %@<struct> (default empty struct) figure options
        fontsize %@(int> fontsize of titles of subplots
        colormap_active %@<string> color map for transform result
        verbose=true     %<logical> default true (use to supress additional output)
        
        % computed values
        C       %@<any> frame coefficents (decomposition of signal)
        frame   %@<Frame> synthesis matrix of frame can be replaced by a fast implicit operation (e.g. wavelets)
        frame_issparse=false %(logical)
        sigma=0  %@<real> noise added to measurements computed from SNR
        frameDepth %@<integer>
        
    end
    
    properties (SetAccess=protected)
        algo           %@<struct> finger print of algorithm
    end
    
    
    %% constructor and commands
    methods
        
        function obj=FrameTrafo(signal)
            % constructor obj=~()
            obj.requireOnly(nargin==0 || isempty(signal) || isa(signal,'SignalClass'),...
                'local', ' signal belongs to class SignalClass');
            if nargin==0
                signal=SignalClass();
            end
            obj.ts=signal;
            obj.C=[];
            obj.tol=1e-9;
            obj.set_algo;
            obj.repfun=@(W) W;  % identity by default
            obj.fig=1;
            obj.colormap_active='default';
            obj.figopt=struct;
            obj.figopt.pixsizeX=1000;
            obj.fontsize=12;
            obj.frameDepth=1;
        end
        
        function set_algo(obj)
            obj.algo.version=0.9;
            obj.algo.versiondate='1.2.2014';
            obj.algo.name=[];
            obj.algo.toolbox=[];
            obj.algo.mfilename=class(obj);
            obj.algo.author='G. Troll';
        end
        
        function reset_Trafo(obj)
            % reset trafo (e.g. after signal has changed)
            obj.C=[];
        end
        
        function set_content(obj,data)
            % ~(data) create new signal with data as content
            obj.ts=obj.ts.make_like();
            obj.ts.set_signal(data);
            obj.reset_Trafo;
            
            obj.ensureOnly(~obj.isTrafodone,'local', ...
                'decomposition is reset');
            obj.ensureOnly(isa(obj.ts,'SignalClass'),'local','2D signal created');
        end
        
        function set_signal(obj,hts)
            % ~(hts) set time signal to hts of class SignalClass
            obj.requireOnly(isa(hts,'SignalClass'),'local','ts is a signal');
            obj.ts=hts;
            obj.reset_Trafo;
            % do not check invariant because caller may be constructor
            obj.ensureOnly(~obj.isTrafodone,'local', ...
                'decomposition is reset');
        end
        
        function set_signalfun(obj,fun,hN,hR)
            % ~(fun,hN,[hR]) set a signal function and sample it at hN points
            obj.requireOnly( isa(fun,'function_handle')...
                && obj.isnatural(hN),'local',...
                'input is function plus sampling length');
            if nargin<3
                hN=100;
            end
            if nargin <4
                hR=2*[1,1];
            end
            if isscalar(hN)
                hN=hN*[1,1];
            end
            constructor= str2func(class(obj.ts));
            obj.ts=constructor(fun,hN,hR);
            obj.reset_Trafo;
            
            obj.ensureOnly(~obj.isTrafodone,'local', ...
                'decomposition is reset');
            obj.ensure(isa(obj.ts,'SignalClass'),'local','2D signal created');
        end
        
        function set_0padding(obj,padd)
            % ~(padd) symm. 0-padding (left and right) by padd;
            % calling with padd=0 removes old padding;
            
            obj.requireOnly(all(arrayfun(@obj.isnatural0,padd)),'local',...
                'padd is pair of non-negative integer');
            obj.ts.set_0padding(padd);
            
            obj.reset_Trafo;
            
            obj.ensureOnly(~obj.isTrafodone,'local', ...
                'wavelet decomposition is reset');
        end
        
        function set_dwtmode(obj, modestr)
            % ~(modestr) border extension mode
        end
        
        function set_basisname(obj,bname)
            % redefined in subclasses
        end
        
        function set_C(obj,cC)
            %~(cC) set transform coefficients ( to save or test recon)
            obj.requireOnly(isvector(cC),'local','cC is vector');
            obj.C=obj.vec2C(cC);
        end
        
        
    end
    
    
    %% queries
    methods
        
        function bn= basisname(obj)
            % bn=~() basis name
            bn=obj.frame.descr;
        end
        
        function str=transformName(obj)
            % name of frame transform
            str=[obj.algo.name,' (',obj.basisname,')'];
        end
        
        function w2=clone(obj)
            % w2=~() clone object (shallow)
            constructor = str2func(class(obj));
            w2=constructor();
            w2.set_signal(obj.ts);
            w2.fig=obj.fig;
            w2.figopt=obj.figopt;
            w2.repfun=obj.repfun;
            w2.use_repfun=obj.use_repfun;
            w2.energy_compressionlevel =obj.energy_compressionlevel;
            w2.tol=obj.tol;
            
            w2.algo=obj.algo;
            
            w2.C=obj.C;
            w2.frame=obj.frame;
            w2.sigma=obj.sigma;
        end
        
        
        function ok=isTrafodone(obj)
            % ok=~(): is frame decomposition of signal available?
            ok=~isempty(obj.C);
        end
        
        function nv=frame_length(obj)
            % n=~() dimension of transform space
            % must be redefined for some subclasses
            nv=obj.frame.length;
        end
        
        function nv=frame_dim(obj)
            % n=~() dimension of Hilbert space
            % must be redefined for some subclasses
            nv=obj.frame.dim;
        end
        
        function nv=signal_dim(obj)
            % n=~() geometric dimension of signal           
            nv=obj.ts.dim;
        end
        
        function r=frame_redundancy(obj)
            % redundancy of the frame used (1 for orthogonal wavelets
            % coeff, max. ~7.8 for 2nd gen. curvelets)
            r=obj.frame_length/obj.frame_dim;
        end
        
        function nv=frame_norm(obj)
            % l2-norm of transform coefficients
            % tests in this form only isometries (e.g. ONS)
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            nv=1;
        end
        
        function F=get_frame(obj)
            % needed because redefined in SampleTrafo2
            F=obj.frame;
        end
        
        function yn=sample(obj,nodes)
            % Y=~() measure transform ftrafo at at nodes
            cC=obj.C2vec;
            yn=cC(nodes.get(1));
        end
        
        function cC=embedSample(obj,nodes,yn)
            % embed samples into image of correct size
            % ~ is right inverse of sample;
            cC=zeros(obj.frame_length,1);
            % set values at nodes
            cC(nodes.get(1))=yn;
        end
        
        function L=deepestlev(obj)
            % deepest level
            L=NaN;
        end
        
        function n=nnz(obj)
            % n=~() number of non-zero elements in the transformation
            % coefficients
            n=nnz(obj.C2vec);
        end
        
        function [sd,p]=sparsityDefect(obj,rs,cC)
            % d=~(rs) rel. sparsity defect of transformation coefficients
            % for relative sparsity rs in [0,1].
            % a.k.a. lp/error of best s/term approximation 
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || isempty(rs) || (rs>=0 && rs<=1),...
                'local','rs is rel. sparsity'); 
            if nargin<3
                cC=obj.C;
            end
            p=1;
            x=reshape(obj.C2vec(cC),[],1);            
            if nargin <2 || isempty(rs)  
                % compute for all s    
                xstar=sort(abs(x),'ascend').^p;
                sd=flipud(cumsum(xstar).^(1/p))/length(x);
            else                
                % absolute sparsity:
                s=max(1,floor(rs*numel(x)));               
                xstar=sort(abs(x),'descend');                
                sd=norm(xstar(s:end),p)/length(x);                     
            end           
        end
 
        
    end
    
    %% transforms
    
    methods
        
        function dec(obj)
            % decomposition (analysis) yields frame coeffs. obj.C
            % can be redefined in subclasses to implement faster
            % or more memory efficient decomposition schemes, e.g.
            % for wavelets.
            % ---> in redefinitions datatype of result must be right for obj.rec,
            % i.e. need not be matrix or vector!
            obj.C=reshape(obj.frame.analyze(obj.ts.xn),obj.N);
        end
        
        function xn=rec(obj,cC)
            % signal reconstruction from frame coefficients C;
            % can be redefined in subclasses to implement faster
            % or more memory efficient decomposition schemes, e.g.
            % for wavelets.
            % ---> result in redefinitions is always a  matrix
            % of the same dimensions as obj.N):
            
            xn=reshape(obj.frame.synthesize(cC),obj.N);
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) apply analyis op. Psi (is one-to-one for a frame)
            % operates mathematically like obj.dec, BUT
            % redefinitions in subclasses must guarantee that:
            % --> INPUT: vector
            % --> OUTPUT: matrix, which can be used as input of synthesize.
            yn=reshape(obj.frame.data'*x(:),obj.N); % cf. frame.analyze 
        end
        
        function xn= synthesize(obj,y)
            % yn=~(x) apply synthesis op. Phi (is onto for a frame)
            % operates mathematically like obj.dec, BUT
            % redefinitions in subclasses must guarantee that:
            % ---> INPUT:  matrix as returned by analyze;
            % ---> OUTPUT:  matrix of the dimensions given by obj.N:
            %xn=reshape(obj.frame.synthesize(y),obj.N);
            xn=obj.frame.data*y(:);   % cf. frame.synthesize(y);
        end
        
        function yn=adjointSynthesis(obj,x)
            % yn=~(x) apply adjoint of synthesis operator Phi to vector x
            % returns vector
            yn=obj.analyze(x);
        end
        
        function x=recAnalysis(obj,yn)
            % x=~(yn) left inverse of analyze
            % for implicit operators (e.g. wavelets) this will be redefined in SampleTrafo
            x=obj.frame.recAnalysis(yn);
        end
        
        function y=recSynthesis(obj,xn)
            % y=~(xn) right inverse of synthesize
            % for implicit operators (e.g. wavelets) this will be redefined in SampleTrafo
            y=reshape(obj.frame.recSynthesis(xn),obj.N);
        end
        
        function y=sim_measurement(obj)
            % simulate a measurement vector y by applying the frame
            % is usually a subsampling,i.e. synthesis (frame.data*signal)
            % for a signal whose length is the frame length.
            obj.require(numel(obj.ts.xn)==obj.frame.length,...
                'local','C is in the frame''s range space');
            y= obj.frame.synthesize(obj.ts.xn(:));
            if ~isempty(obj.SNR) && obj.SNR>0
                % add SNR [dB] noise to measurement
                obj.sigma = std(y(:))*10^(-obj.SNR/20);
                y=y+obj.sigma*randn(size(y));
            end
            obj.ensureOnly(numel(y)==obj.frame.dim,'local', 'result has same dim as frame');
        end
        
        function ss=measurement2graph(obj,yn)
            % ss=~(yn) transforms measurement vector to output form.
            % In this class it just wraps it into a 1d signal;
            % method is however redefined in SampleTrafo.
            % yn ... measurement vector, e.g. output of sim_measurement
            obj.requireOnly(isvector(yn) && obj.frame.dim==length(yn),...
                'local', 'yn is measurement vector');
            ss=obj.TimeSignal(yn);
            ss.repfun=obj.ftrafo.repfun;
        end
        
        function vC=C2vec(obj,cC)
            % vc= ~([cC]) convert frame coefficients to vector form
            % must be redefined for those subclasses, where C is not a
            % matrix but e.g. a cell array.
            obj.requireOnly(isnumeric(obj.C),'local','column vector defined for input. Redefine method?');
            if nargin <2
                vC = obj.C(:);
            else
                vC=cC(:);
            end
            obj.ensureOnly(isvector(vC),'local','result is column vector');
        end
        
        function nc=numelC(obj,cC)
            % nc=~([cC]) number of elements in decomposition
            % must be redefined for those subclasses, where C is not a
            % matrix but e.g. a cell array.
            if nargin <2
                nc=numel(obj.C);
            else
                nc=numel(cC);
            end
        end
        
        function cC=vec2C(obj,vC)
            % convert vector to form usable for reconstruction rec
            cC=vC;
        end
        
        function [mat, levelshown]=C2graph(obj, lev)
            % mat=~() convert coefficient vector to format suitable for
            % graphical output (usually rearrangement of elements, e.g.
            % for wavelet tansform);
            % will be redefined in child classes
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            levelshown=[];
            mat=obj.C;
        end
        
        function computeFrame(obj,hN, cutoff)
            % ~([hN,cutoff]) compute frame matrix (matrix form of decomposition)
            % calling the algorithm implemented in dec (using linearity
            % only, i.e. row(!) vector v'_i of frame matrix is the result of
            % applying the linear dec operation on the canonical basis
            % e_i, so that with frameMat=[v_i] we get dec(x)== frameMat*x(:)).
            %
            % Enter hN if the frame is to be calculated of a particular size without respect to
            % the signal obj.ts and if it is also a basis;
            % if no hN is given, the dimensions are taken from the signal obj.ts and its
            % transform obj.C.
            % cutoff determines at which size coefficients are to be
            % regarded as 0.
            % Remark: long run-time !
            % ex.: WavePath; w=Wavelet2D_wlab();  w.set_basisname('db1');
            %      hN=16; w.computeFrame(hN); w.graph_trafo;
            %
            
            obj.requireOnly(nargin<2 ||~any(hN==0),'local', 'need positive size');
            if ~exist('hN','var')
                hN=obj.N;
            end
            if  any(hN==0)
                dim=length(obj.ts.origin);
                hN=max(16, 2^(7-dim));   % default size decreases with dim
            end
            
            if nargin <3
                cutoff=1e-10;
            end
            
            if length(hN)~=length(obj.ts.origin)  % ts.size has length 2 if empty
                hN=hN(1)*ones(1,length(obj.ts.origin));
            end
            
            n=prod(hN);  % dimension of signal space
            
            % smallest size of transform to ensure
            % in the case of matlab's wavelet toolbox
            % that all transformed signals have the same size k.
            obj.set_dwtmode('per');
            
            % get dimensions from transform of some basis vector of signal
            % space;
            % for basis transforms (frame is basis): k==n (e.g. wavelets
            % with periodic boundary conditions)
            % but for general frames (e.g. curvelets): k>n
            % We assume here that the parameters of the decomposition
            % have been chosen, s.t. we get the same k for all basis
            % vectors ej.
            SP=obj.ts.make();  % create correct subclass of SignalClass
            ej=SP.make_delta(floor(n/2),hN,0); % some basis vector
            obj.set_signal(ej);
            obj.dec;
            k=obj.numelC; %numel(obj.C2vec);  % dimension of decomposition space
            
            % obj.frame=Frame([k,n]); % k frame vectors of dim n
            if obj.frame_issparse
                obj.frame=FrameSparse();
                % generate the three vectors i_s, j_s, n_s
                % required to create a sparse matrix efficiently
                i_s=cell(n,1);
                j_s=i_s; n_s=i_s;
            else
                obj.frame=Frame([n,k]);
                obj.frame.data=zeros(n,k);
            end
            
            wait_handle = waitbar(0,...
                ['computing frame (dim=',num2str(n),') ...']);
            
            for i=1:n
                
                waitbar(i / n);
                
                ei=SP.make_delta(i,hN,0); % new basis vector
                obj.set_signal(ei);
                % decompose ei; result is colum of analysis op.,
                % i.e. row of synthesis op (=frame):
                obj.dec;
                % set row vector of frame to obj.C:
                cC=obj.C2vec;
                if numel(cC)==k
                    if ~obj.frame_issparse
                        obj.frame.data(i,:)=cC(:);
                    else
                        j_s{i}=find(abs(cC)>cutoff); % non-zeros
                        i_s{i}=i*ones(length(j_s{i}),1);
                        n_s{i}=cC(j_s{i});
                    end
                else
                    error('frame redundancy k/n is not constant');
                end
            end
            close(wait_handle );
            
            if obj.frame_issparse
                i_s=cell2mat(i_s);
                j_s=cell2mat(j_s);
                n_s=cell2mat(n_s);
                obj.frame.sparse(i_s,j_s,n_s); % instead of filling sparse matrix elementwise
            else
                % other case already handled in i-loop above
            end
            
            obj.frame.descr=['frame ',obj.basisname,...
                '(',obj.title2, ', ',obj.algo.toolbox,')'];
            obj.ensure(~isempty(obj.frame),'local','frame generated');
        end
        
    end % computeFrame
    
    %% filters
    
    methods
        
        function estimateNoise(obj)
            % estimate noise by measuring the standard deviation of the first
            % diagonal fluctutation of an orthogonal(!) wavelet transform of NoisyImage.
            % observe correction factor for noise when thresholding:
            % used as threshold: thresh= obj.frame_norm *4*obj.sigma
            try
                w=Wavelet2D_mlab();
            catch
                w=Wavelet2D_wlab();
            end
            w.set_signal(obj.ts);
            w.set_basisname('db2');
            w.dec;
            w.estimateNoise;
            obj.sigma=w.sigma;
        end
        
        function hardthresh1_inplace(obj,thresh)
            % filter all coefficients obj.C
            obj.C(abs(obj.C)<thresh)=0; % operation on place to save memory
        end
        
        function hardthresh2_inplace(obj,thresh1, thresh2)
            % filter with 2 thresholds thresh1<thresh2
            obj.requireOnly(thresh1<=thresh2,'local', 'thresh1<=thresh2');
            % redefine in subclasses
            obj.C(abs(obj.C)<thresh2)=0;
        end
        
        function softthresh1_inplace(obj,thresh)
            % simple soft filter: reduce all coeff. obj.C by size
            obj.C(abs(obj.C)<thresh)=0;
            % now also reduce other coefficients
            obj.C=obj.C-sign(obj.C)*thresh; % operation on place to save memory
        end
        
        function [signal2, thresh, fCoeff]=hardthresh1(obj,thresh)
            % [signal2,thresh.fCoeff]=~([thresh]) hard threshold transform coefficients  
            % e.g. usage for denoising: obj.hardthresh() uses noise estimator.
            % input:
            %   thresh  ... threshold to be used (if missing automatic
            %               denoising threshold is used).
            % output:
            %   signal2 ... reconstructed signal after thresholding
            %               coefficients obj.C.
            %   thresh  ... threshold used
            %   fCoeff  ... thresholded coefficients
            % comments:
            % -- cf. hardthresh2 using 2 thresholds for different kinds of
            %    coefficients.
            % -- hardthresh1 is used by sparseApprox; 
            % -- subclasses redefine hardthresh1_inplace if necessary.
            % -- observe correction factor obj.frame_norm for noise (e.g.
            %    curvelets).
            %
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin < 2 || isempty(thresh)
                obj.estimateNoise();
                thresh=4*obj.frame_norm *obj.sigma;
            end
            
            oldC=obj.C; % must be clone
            obj.hardthresh1_inplace(thresh);
            
            % reconstructed signal
            yn = obj.rec(obj.C) ;
            
            % renormalisation prevents calculation of PSNR, hence commented
            % out:
            %           yn = yn-min(yn(:));
            % 			yn = yn/max(yn(:));
            
            if nargout >2
                fCoeff=obj.C;
            end
            obj.C=oldC;   % restore
            signal2=obj.ts.clone();
            signal2.replace_signal(yn);
            signal2.signalname=[obj.ts.signalname,', hard thresh. (',obj.basisname,...
                ', ',num2str(thresh,'%3.1e'),')'];
            signal2.colormap_active=obj.ts.colormap_active;
            
        end
        
        function [signal2, thresh1, thresh2,fCoeff]=hardthresh2(obj,thresh1, thresh2)
            % [signal2,thresh1,thresh2]=~([thresh1,thresh2]) hard threshold transform coefficients
            % (cf. hardthresh1) using 2 thresholds for different kinds of
            % coefficients; subclasses redefine hardthresh2_inplace if necessary.
            % use cases:
            % --- denoise: obj.hardthresh2() uses noise estimator
            % --- sparsify: hardthresh1 is used by sparseApprox; 
            % observe correction factor obj.frame_norm for noise
            % propagation.
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if  nargin<2 || isempty(thresh1) || nargin <3
                obj.estimateNoise;
            end
            if nargin<2 || isempty(thresh1)
                thresh1=4*obj.frame_norm *obj.sigma;
            end
            if nargin < 3
                thresh2=3*obj.frame_norm *obj.sigma;
            end
            assert(thresh1>=thresh2,'thresh1 is for finest details, thresh2 for coarsest');
            
            oldC=obj.C; % must be clone
            obj.hardthresh2_inplace(thresh1,thresh2);
            
            % reconstructed signal
            yn = obj.rec(obj.C) ;
            % renormalisation prevents calculation of PSNR, hence commented
            % out:
            %           yn = yn-min(yn(:));
            % 			yn = yn/max(yn(:));
            if nargout >3
                fCoeff=obj.C;
            end
            obj.C=oldC;   % restore
            signal2=obj.ts.make_like();
            signal2.replace_signal(yn);
            signal2.signalname=[obj.ts.signalname,', 2 hard thresh. (',obj.basisname,...
                ', ',vec2str([thresh1,thresh2],'%3.1e'),')'];
            signal2.colormap_active=obj.ts.colormap_active;
            
        end
        
        function [signal2, thresh, fCoeff]=softthresh1(obj,thresh)
            % [yn,thresh]=~(thresh) soft threshold transform coefficients,
            % e.g. to desnoise signal: obj.softthresh() uses noise
            % estimator;
            % cf. hardthresh1,2
            % subclasses redefine softthresh1_inplace if necessary.
            % observe correction factor obj.frame_norm for noise
            % propagation:
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin < 2 || isempty(thresh)
                obj.estimateNoise();
                thresh=4*obj.frame_norm *obj.sigma;
            end
            
            oldC=obj.C; % must be clone
            obj.softthresh1_inplace(thresh);
            
            % reconstructed signal
            yn = obj.rec(obj.C) ;
            
            if nargout >2
                fCoeff=obj.C;
            end            
            obj.C=oldC;   % restore
            
            signal2=obj.ts.make_like();
            signal2.replace_signal(yn);
            signal2.signalname=[obj.ts.signalname,', soft thresh. (',obj.basisname,...
                ', ',num2str(thresh,'%3.1e'),')'];
            signal2.colormap_active=obj.ts.colormap_active;
            
            
        end
        
        function thresh=computeThreshold(obj,cs)
            % thresh=~(cs) computes treshold at compression rate cs
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(cs>1,'local','compression rate >1');
            
            sC=sort(abs(obj.C2vec),'descend');            
            % threshold is chosen s.t. 1/cs of all elements are bigger
            thresh=sC(ceil(numel(sC)/cs));
        end
        
        function [signal2,BSE,s,ortho_defect,fCoeff]=sparseApprox(obj,cs,p)
            % [signal2,BSE,s]=~(c_s,p) best p-norm s-term approx. at compression c_s
            % signal2 ... approximation of signal obj.ts by best s-terms.
            % BSE ... lp-error of best s-term approximation (sparsity defect)
            % s ... sparsity of filtered transform
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(cs>1,'local','compression rate >1');
            if nargin < 2 || isempty(cs)
                cs=16;
            end
            if nargin <3 || isempty(p)
                p=2;  % lp-norm
            end
            sC=sort(abs(obj.C2vec),'descend');
            obj.check(sC,'all(local>=0)','sC is non-neg. vector');
            
            % threshold is chosen s.t. 1/cs of all elements are bigger
            thresh=sC(ceil(numel(sC)/cs));
            
            if nargout<5
                signal2=obj.hardthresh1(thresh);
            else % return also filtered coefficients
                [signal2,~,fCoeff]= obj.hardthresh1(thresh);
            end
            signal2.signalname=[obj.ts.signalname,', hard thresh. (',obj.basisname,...
                ', c=',num2str(cs,'%3.1f'),')'];
            
            if nargout>=2   % ||W-F(W)||/||W||
                BSE=(sum((sC(sC<thresh)).^p))^(1/p);   % needs original values
            end
            if nargout>=3
                % sparsity s should be ca. numel(obj.C)/cs;
                s=sum(sC<thresh);  % needs non.neg. values
            end
            if nargout>=4
                % for orthognal transforms ortho_defect=0
                img_norm=obj.ts.norm(2);
                ortho_defect=abs(norm(sC,2)-img_norm)/img_norm; % needs original values
            end
            
        end
        
        function [nodes,mask]=nodesPDF(obj,c, params)
            % nodes=~(c) nodes sampled with compression rate c
            obj.requireOnly(c>1,'local','compression rate in ]0,1[');
            obj.requireOnly(~obj.ts.isemptydata,'local','non-empty signal');
            DN = obj.ts.N; % data Size
            p0 = 1/c;% [0.25];% undersampling factor
            if nargin<3 || ~isfield(params,'deg')
                params.deg=5; % Variable density polymonial degree
            end
            if ~isfield(params,'iter')
                params.iter=10;
            end
            % generate variable density random sampling
            pdf = genPDF(DN,p0,params);	% generates the sampling PDF
            % generate sampling mask:
            mask = genSampling(pdf,params);
            % compute nodes:
            nodes=SampledNodes(find(mask),DN,obj);
        end
        
        
    end
    
    %% frame operations
    
    methods
        
        function cC=dec2UsingFrame(obj)
            % cC=~() decompose using matrix representation of frame (x_i).
            % result should be numerically the same as that of dec;
            % cf. test_frame.
            % obj.frame can be given explicitly or can be the result of computeFrame
            % if dec does not use the frame in some subclass like wavelet.
            obj.require(true,'local','verify invariant');
            cC=reshape(obj.frame.analyze(obj.ts),obj.N);
        end
        
        function s=rec2UsingFrame(obj, cC)
            % s=~(cC) reconstruct signal from frame coefficients cC
            % result should be numerically the same as that of rec.
            % obj.frame can be given explicitly or can be the result of computeFrame
            % if dec does not use the frame in some subclass like wavelet.
            obj.require(isvector(cC),'local','argument is vector');
            s=reshape(obj.frame.rec1(cC),obj.N);
        end
        
        function M=synthMat(obj)
            % thetastar=~() returns synthesis matrix theta' of frame
            % i.e. frame matrix
            M=obj.frame.synthMat;
        end
        
        function theta=analMat(obj)
            % theta=~() returns analysis matrix theta of frame
            % i.e. adjoint frame matrix
            theta=obj.frame.analMat;
        end
        
        function F2=RandnMatTimesDualFrame(obj,m,n)
            % F2=(m,n) multiply random matrix randn(m,n) with dual frame
            % used e.g. for compressive sensing problems
            % calculates A*S^(-1)*theta' with a random matrix A.
            obj.require(rem(sqrt(n),1)==0,'local', 'n is square number');
            F2=Frame(randn(m,n)*obj.frame.dualFrame.data);
        end
        
    end
    
    
    methods (Hidden)
        
        function ss=N(obj)
            % sample size (of image)
            ss=obj.ts.N;
        end
        
        function xn=sdata(obj)
            % xn=~() signal data
            xn=obj.ts.xn;
        end
        
        function ss=K(obj)
            % size of frame coefficients in rectangle form
            % squares
            ss=sqrt(obj.frame.length)*ones(1,2);
            if ~all(rem(ss,1)==0)
                % rectangle closest to square
                s=obj.frame.length;a
                % find all factors:
                D=[1,unique(cumprod(perms(factor(s)),2))];
                % find factor closest to root
                pos=find(D>ss,'first');
                ss=[D(pos-1),D(pos)];
            end
        end
        
        function ne=numel(obj)
            % # elements in signal
            ne=obj.ts.numel;
        end
        
        function hR=Rs(obj)
            % sampled rectangle
            hR=obj.ts.Rs;
        end
        
        function v=padd(obj)
            % zero padding size
            v=obj.ts.padd;
        end
        
    end
    
    %% static
    methods (Static)
        
        function str=msgDoTrafo()
            str='decomposition available (call dec).';
        end
        
    end
    
    
    %% tests
    methods
        
        
        function [ok,err1,err2]= test_DecRec(obj)
            % [ok,res]= ~() tests if dec is the inverse of rec
            obj.require(~obj.ts.isemptydata,'local','signal has content');
            
            oldsignal=obj.ts;
            
            obj.dec;
            y=obj.rec(obj.C);
            err1=abs(max(abs(y(:)-obj.ts.xn(:))));
            cC1=obj.C2vec;
            
            ts2=obj.ts.make_like();
            ts2.set_signal(y);
            obj.set_signal(ts2);
            
            obj.dec;
            err2=abs(max(abs(obj.C2vec)-abs(cC1)));
            
            ok=max(err1,err2) <=obj.tol;
            
            obj.ts=oldsignal;
            
        end
        
        function [ok,err1,err2]= test_AnalyzeSynthsize(obj)
            % [ok,res]= ~() tests if analyze is the inverse of synthesize
            obj.require(~obj.ts.isemptydata,'local','signal has content');
            
            cC1= obj.analyze(obj.ts.xn);
            xn2= obj.synthesize(cC1);
            err1=abs(max(xn2(:)-obj.ts.xn(:)));
            
            cC2=obj.analyze(xn2);
            err2= max(abs(cC1(:)-cC2(:)));
            
            ok=max(err1,err2) <=obj.tol;
            
        end
        
        function res=test_framefeatures(obj,hN)
            % ok=~([cmdstr,n]) test some frame properties such as unitarity;
            % returns frame features as structure res.
            obj.requireOnly(nargin<2 || isnumber(hN),'local', 'needs frame size');
            dim=length(obj.ts.origin);
            if ~exist('hN','var') || isempty(hN)
                hN=max(16, 2^(7-dim));   % default size decreases with dim
            end
            if length(hN)~=dim
                hN=hN(1)*ones(1,dim);
            end
            if isempty(obj.frame) || ~isequal(obj.frame.size(1),prod(hN))
                obj.computeFrame(hN);
            end
            
            fprintf('.');
            res=struct;
            res.name=[obj.algo.name,' (',obj.basisname(),')'];
            res.test_size=hN;
            res.isunitary=obj.frame.isunitary;
            if ~res.isunitary
                res.isbasis=obj.frame.isbasis;
            else
                res.isbasis=true;
            end
            
            fprintf('.');
            res.redundancy=obj.frame_redundancy;
            if ~res.isunitary
                res.isparseval=obj.frame.isparseval;
            else
                res.isparseval=true;
            end
            
            fprintf('.');
            if ~res.isparseval
                [res.istight,~,FB]=obj.frame.istight;
                if res.istight
                    res.framebound=FB;
                end
            else
                res.istight=true;
            end
            
            fprintf('.');
            if res.redundancy==1 && obj.frame.issquare
                res.isinvolution=obj.frame.isidempotent;
            else
                res.isinvolution=false;
            end
            signal=obj.ts.make_like();
            if length(hN)==1
                signal.set_signal(randn(hN,1));
            else
                signal.set_signal(randn(hN));
            end
            obj.set_signal(signal);
            obj.dec;
            y=obj.rec(obj.C);
            err=max(abs(y(:)-obj.ts.xn(:)));
            
            res.isrecInverseOfdec=err<obj.frame.tol_deficit;
            disp('.');
        end
        
        function test_frame(obj)
            % ~() test if computeFrame generates the same transform as dec
            % obj.frame of method obj.computeFrame;
            % as the runtime of computeFrame is long, frame should be
            % precalculated and entered as a parameter of this method.
            obj.require(~isempty(obj.frame),'local','needs frame (e.g. computed by computeFrame)');
            obj.require(obj.isnaturalArray(sqrt(size(obj.frame))),'local',...
                'frame operates on square matrix');
            
            w2=obj.clone;
            w2.set_dwtmode('per');
            hN=sqrt(size(obj.frame));
            imagesize=[hN(1),hN(1)];
            if any(obj.N==0) || obj.ts.density<0.1
                % replace by test image if image does not fit frame
                hts=obj.ts.make_standardExample();
                w2.set_signal(hts);
            elseif ~isequal(hN,obj.N)
                hts=obj.ts.clone;
                hts.resize(imagesize);
                w2.set_signal(hts);
            end
            
            w2.dec;
            
            prepfigure(1);
            hts=w2.ts;
            % apply analysis operator to signal and compare with
            % decompistion result of dec:
            diff=reshape(obj.frame.analMat*hts.xn(:)-w2.C2vec,w2.K);
            maxdiff=max(abs(diff(:)));
            meandiff=mean(abs(diff(:)));
            imagesc(real(diff));   % use real in case of complex frames
            colorbar;
            title({[upper(hts.get_signalname),...
                ': difference \Delta between matrix operation and decomposition operation'],...
                [w2.basisname,'(',w2.title2,') ',w2.algo.toolbox],...
                ['\mu(\Delta)=',num2tex(meandiff,'%3.1e','none'),', max(|\Delta|)=',...
                num2tex(maxdiff,'%3.1e','none')]},...
                'fontsize',obj.fontsize);
            
        end
        
    end
    
    %% graphics
    methods (Hidden)
        
        function str=add_signalname(obj)
            str=[];
            sn=obj.ts.signalname;
            if ~isempty(sn)
                str=[', ',sn];
            end
        end
        
    end
    
    methods
        
        function tit2= title2(obj,varargin)
            tit2=[];
        end
        
        function cC=dec2graph(obj,cC)
            % cC=~(cC) post process decomposition obj.dec for graphical output
            % redefined e.g. for fft as fftshift
        end
        
        function graph_signal(obj, open_new)
            % show signal in the open window (default)
            obj.requireOnly(isa(obj.ts,'SignalClass'),'local','signal exists');
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            obj.ts.graph_signal(false);
        end
        
        function graph_recon(obj,open_new)
            % ~() show signal reconstructed from transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin<2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            RS=obj.ts.make_like();
            RS.set_signal(reshape(obj.rec(obj.C),obj.ts.N));
            PSNR=obj.ts.PSNR(RS);
            RS.signalname={['Recon by ',obj.basisname,': PSNR=',stat2str(PSNR,'%3.1f'),' dB'],...
                obj.ts.signalname};
            RS.graph_signal(false);
        end
        
        function show_recon(obj)
            % ~() show signal reconstructed from transform in new window
            obj.graph_recon(true);
        end
        
        function graph_trafo(obj, open_new, maxsize)
            % ~() show frame of transform (computed by computeFrame) in open window
            obj.require(~isempty(obj.frame),'local','needs frame (e.g. computed by computeFrame)');
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            if nargin <3
                maxsize=1e3;
            end
            obj.frame.graph_data(false,maxsize);
            colormap(obj.colormap_active);
            SNR_str=[num2str(obj.SNR,'%3.0f'),' dB'];
            if isempty(obj.SNR)
                SNR_str='no noise';
            end
            title([obj.basisname,' - ',SNR_str],'fontsize',12);
            
        end
        
        function graph_distribution(obj, open_new,y)
            % ~() show distribution of frame in open window
            obj.require(~isempty(obj.frame),'local','needs frame (e.g. computed by computeFrame)');
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            if nargin <3
                y=[];
            end
            obj.frame.graph_distribution(false);
            SNR_str=[num2str(obj.SNR,'%3.0f'),' dB'];
            if isempty(obj.SNR)
                SNR_str='no noise';
            end
            title([obj.basisname,' - ',SNR_str],'fontsize',12);
            
        end
        
        function graph_sparsityDefect(obj,open_new,c)
            % ~() sparsity defect (l2-error of best s-term approx) of transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            if nargin<3
                c=[];
            end
           s=TimeSignal(reshape(obj.C2vec,[],1));
           s.signalname=[obj.algo.name,obj.add_signalname,' ',vec2str(obj.ts.size)];
           s.graph_sparsityDefect(false,c);
            
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal class';
            ok= isa(obj.ts,'SignalClass');
            if ok
                descr='frame empty or exists as frame object';
                ok=isempty(obj.frame) || isa(obj.frame,'Frame');
            end
        end
    end
    
    
end

