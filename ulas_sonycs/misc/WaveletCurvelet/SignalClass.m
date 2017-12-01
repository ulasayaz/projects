classdef SignalClass<DC
    % parent class representing n-dim. signals
    %
    % G. Troll, Berlin 2013
    %
    
    
    properties
        xn         %@<matrix> discrete signal data
        sigEst=0   %@<double) estimated noise std
        lp_norm=2  %@<integer> norm used in this object
        padd=0     %@<pair> (int), default (0,0), length of 0-padding
        paddDir    %@<string> padd direction
        signalname %@<string)
        winfun     %@<function_handle>, default empty, windowing function, e.g. @chebwin
        
        amp_quantity %<string> quantity of amplitude
        amp_unit     %<string> unit of amplitude
        fig
        colormap_active     %@<string> color map for matrix output
        colormap_freeze     %@<logical> freeze colormap?
        colorlimits         %@<pair> color limits used for caxis
        ison_colorbar       %@<logical> w/o colorbar
        repfun  %@<function_handle> representation function for graphical output of signal
        hide_outliers=false %@<logical> hide outliers in graph_signal?
        borderSize=0        %@<integer> size of border in pixels
    end
    
    properties (SetAccess=protected)
        units   %@<cellstr> e.g. '[mm]', default is {}
        scale   %@<PosScales> n, t , x, y        
        origin       %@<pair> origin of coord. system => domain=[z-Rs/2, z+Rs/2]
    end
    
    
    
    %% constructor and commands
    methods
        
        function set_scale(obj,sc)
            % e.g. set_scale(PosScales.n)
            obj.requireOnly(isa(sc,'PosScales'),'local','enumeration type');
            obj.scale=sc;
            obj.ensureOnly(isa(obj.scale,'PosScales'),'local','enumeration type');
        end
        
        function replace_signal(obj,xn)
            % without reset
            obj.requireOnly(isnumeric(xn) ,'local','input is numeric');
            obj.xn=xn;
        end
        
        function set_signal(obj,xn,hR)
            % ~(xn,[hR]) set signal to matrix xn
            obj.requireOnly(isnumeric(xn) ,'local','input is matrix');
            obj.reset_Signal;
            obj.xn=xn;
        end
        
        function set_automaticStyle(obj)
            % ~() set graphics style adapted to data;
            % redefine for signals of different dimensionality
            
        end
        
        function reset_Signal(obj)
            % reset signal trafo (e.g. after signal has changed)
            % redefined in subclasses
            
        end
        
        function set_winfun(obj,wf)
            % ~(wf) set windowing function to wf, e.g. wf=@chebwin
            obj.requireOnly( isempty(wf)|| isa(wf,'function_handle'),'local',...
                'input is empty or function handle');
            obj.winfun=wf;
        end
        
    end
    
    %% transformation commands
    
    methods
        
        function set_mean0(obj)
            % subtracts mean from signal, new signal has mean 0.
            obj.xn=obj.xn-mean(obj.xn(:));
        end
        
        function normalize(obj,do_meannormalize)
            % ~() normalize signal to peak-valley 1 and either mean value 0
            % or min value 0.
            if nargin <2
                do_meannormalize=true;
            end
            if ~obj.isemptydata()
                obj.xn=obj.xn/obj.pv; % peak-valley normalized to 1
                if do_meannormalize
                    obj.xn=obj.xn-meanNaN(obj.xn(:)); % mean value normalized to 0
                else
                    obj.xn=obj.xn-min(obj.xn(:));
                end
            end
        end
        
        function padd2dyadic(obj, pval, direction)
            % ~([pval] pad signal to next dyadic hypercube with value pval
            obj.require(obj.numel>0,'local','signal is nonempty');
            if ~exist('pval','var')
                pval=0;
            end
            if ~exist('direction','var')
                direction='both';  % both, post, pre
            end
            ss=max(2.^ceil(log2(obj.size)))*ones(1,length(squeeze(obj.size)));
            obj.padd=ss-obj.size;
            if strcmp(direction,'both')
                padsize=floor(obj.padd/2);
                obj.xn=padarray(obj.xn,padsize,pval,direction);
                rest=padamount-2*padsize;
                if max(rest)>0
                    obj.xn=padarray(obj.xn,rest,pval,'post');
                end
            else
                obj.xn=padarray(obj.xn,obj.padd,pval,direction);
            end
            obj.ensureOnly(obj.isdyadicHypercube(),'local',...
                'signal dimensions are equal and dyadic');
        end
        
        function crop(obj,ss)
            % ~([ss]) crop to size ss (default lower dyadic bound)
            if nargin <2
                ss=min(2.^floor((log2(obj.size))));
            end
            ssold=size(obj.xn);
            if isscalar(ss) && ~isvector(obj.xn)
                ss=ss*ones(1,length(ssold));
            end
            ss=min(ss,ssold);
            sn=obj.signalname;
            obj.set_signal(wkeep(obj.xn,ss));  % central part
            obj.signalname=sn;
        end
        
        function shiftSignalValues(obj,amount)
            % ~(amount) shift signal values by amount
            obj.xn=obj.xn-amount;
        end
        
        function addNoise(obj,sig)
            %~(sig) add Gaussian noise with std sig
            obj.requireOnly(nargin>1 && sig>=0,'local','sig is non-negative');
            obj.xn=obj.xn+sig*randn(size(obj.xn));
        end
        
        function outliers2NaN(obj,q, nstd)
            % ~(q) replace outliers by NaN
            % outlier condition cf. function outliers
            obj.requireOnly(q>0 && q<1,'local','0<q<1');
            if nargin <3
                nstd=3;
            end
            [thresh1, thresh2]=obj.outliers(q, nstd);
            
            obj.xn(obj.xn>thresh2)=NaN;
            obj.xn(obj.xn<thresh1)=NaN;
            
        end
        
        function sparsify(obj, ct)
            % ~(N,ct) hardthreshold to sparsity obj.numel/ct
            obj.requireOnly(ct>1,'true','sparsity compression >1');
            ss=sort(abs(obj.xn(:)),'descend');
            thresh=ss(ceil(numel(ss)/ct));
            obj.xn(abs(obj.xn)<thresh)=0;
            descr=[' sparsified c_t=',num2str(ct,'%3.1g')];
            obj.signalname=strrep(obj.signalname,descr,''); % prevent doubling
            obj.signalname=[obj.signalname,descr];
        end
        
    end
    
    %% queries
    methods
        
        function s2=clone(obj)
            % clone object
            s2=obj.make_like();  % make same class
            s2.xn=obj.xn;
        end
        
        function s2=make_like(obj)
            % s2=~() clone all but signal content
            s2=obj.make();  % make same class
            s2.units =obj.units;
            s2.origin=obj.origin;
            s2.signalname=obj.signalname;
            s2.hide_outliers=obj.hide_outliers;
            s2.colormap_active=obj.colormap_active;
            s2.colorlimits=obj.colorlimits;
            s2.colormap_freeze=obj.colormap_freeze;
            s2.scale=obj.scale;
            s2.padd=obj.padd;
            s2.paddDir=obj.paddDir;
            s2.winfun=obj.winfun;
        end
        
        function ok=isemptydata(obj)
            %~() tests if signal content is empty
            ok=isempty(obj.xn);
        end
        
        function ss=size(obj,k)
            % s=~(): sample size
            ss=size(obj.xn);
            if length(ss)==2
                ss=[ss,1];
            end
            if nargin>1
                ss=ss(k);
            end
        end
        
        function ss=size_dyadic(obj,k)
            % s=~(): ceiling of log2 of sample size
            if nargin>1
                ss=ceil(log2(size(obj.xn,k)));
            else
                ss=ceil(log2(size(obj.xn)));
            end
        end
        
        function ok= isdyadicHypercube(obj)
            % ok=~()check if all dimensions are equal and dyadic
            ok=obj.numel==0; % must admit empty case for diverse invariants
            if ~ok
                ok=all(2.^ceil(log2(obj.size))-obj.size==0);
                ok=ok && all(obj.size == obj.size(1));
            end
        end
        
        function n=numel(obj)
            % number of elements in 2D signal (image)
            n=numel(obj.xn);
        end
        
        function n=dim(obj)
            % n=~() geometric dimension of signal
            n=sum(obj.N>1);
        end
        
        function s=sparsity(obj,heps)
            % s=~([heps)] sparsity yof signal
            if nargin <2
                heps=eps;
            end
            s=sum(abs(obj.xn(:) )>heps);
        end
        
        function [sd,p]=sparsityDefect(obj,rs)
            % sd=~([rs]) s-sparsity defect (lp-error of best s-term approx.) of signal
            % with optional parameter relative sparrsity rs=s/obj.numel.
            % sigma_s(x)_p=inf{||x-z||; z \in C^N s-sparse}=(\sum_{i>s} |x_i^*|^^p)^{1/p}
            obj.requireOnly(rs>=0 && rs<=1,'local','rs is rel. sparsity');
            p=obj.lp_norm;
            if nargin <2 || isempty(rs)
                % compute for all s
                x=sort(abs(obj.xn(:)),'ascend').^p;
                sd=flipud(cumsum(x).^(1/p));
            else
                cc=numel(obj.xn);
                % absolute sparsity:
                s=max(1,floor(rs*cc));
                xstar=sort(abs(obj.xn(:)),'descend');
                sd=norm(xstar(s:end),p);
            end
            
        end
        
        function activeSet=supp(obj, heps)
            % as=~([eps]) support of signal
            if nargin <2
                heps=eps;
            end
            activeSet=find(abs(obj.xn) > heps);
        end
        
        function n=nnz(obj)
            % number of non-zero elements
            n=nnz(obj.xn);
        end
        
        function d= density(obj)
            % density of non-zero elements
            d= nnz(obj.xn)/numel(obj.xn);
        end
        
        function sn=get_signalname(obj)
            % name of signal
            if isempty(obj.signalname)
                sn= '';
                if isprop(obj,'xt')
                    sn=func2str(obj.xt);
                end
            else
                sn=obj.signalname;
            end
            if ~isempty(sn)
                if ~iscell(sn)
                    if ~strfind(sn,'_{')  % do not replace latex code
                        sn=strrep(sn,'_','\_');  % to avoid subscript
                    end
                end
            end
        end
        
        function sn=shortname(obj)
            % short name of signal (cut off after 1st non-alphabetic char)
            sn=obj.signalname;
            sn=regexp(sn,'[a-zA-Z0-9]*','once','match');
        end
        
        function [ok_pv,ok_mv,err]= isnormalized(obj)
            % [ok_pv, ok_meanv] =~() tests if signal has been normalized
            % ok_pv ... peak-valley normalized to 1
            % ok_mv ... mean value normalized to 0
            err(1)= abs(obj.pv-1);
            err(2)= abs(obj.mean);
            ok_pv= err(1) <=1e-6;
            ok_mv= err(2) <1e-6;
        end
        
        function s= winfun2str(obj)
            % string identifuer of winfun
            if ~isempty(obj.winfun)
                s=func2str(obj.winfun);
            else
                s='rect';
            end
        end
        
        function ts2=apply_winfun(obj)
            % ts2=~() return windowed signal;
            obj.requireOnly(~isempty(obj.winfun),'local','winfun is set');
            weights=obj.eval_winfun;
            ts2=obj.make_like();
            ts2.replace_signal(weights.*obj.xn);
            ts2.signalname=[obj.signalname,' ',obj.winfun2str,'-windowed'];
        end
        
        function ts2=apply_winfunInverse(obj)
            % ts2=~() return windowed signal;
            obj.requireOnly(~isempty(obj.winfun),'local','winfun is set');
            weights=1./obj.eval_winfun;
            ts2=obj.make_like();
            ts2.replace_signal(weights.*obj.xn);
            ts2.signalname=[obj.signalname,' ',obj.winfun2str,'-unwindowed'];
        end
        
        function weights= eval_winfun(obj)
            % calculate weights from windows function
            % implement in subclasses
            weights=ones(obj.N);
        end
        
        function [smat,BSE,s]=signal2sparse(obj,ct)
            % [smat,BSE,s]=~(c_s) compute best s-sterm approx. at compression factor c_s
            % smat ...best s-term approx. of signal
            % BSE ... error of best s-term approximation (sparsity defect)
            % s ... sparsity of filtered curvelet trafo
            % apply e.g. to a test signal created by w= Signal2D.make_sinus(5,0,0.1).
            obj.requireOnly(ct>1,'local','compression rate >1');
            if nargin < 2 || isempty(ct)
                ct=4;
            end
            
            sC=sort(abs(obj.xn(:)),'descend');
            obj.check(sC,'all(local>=0)','sC is non-neg. vector');
            
            % threshold is chosen s.t. 1/ct of all elements are bigger
            thresh=sC(ceil(numel(obj.xn)/ct));
            
            smat=zeros(size(obj.xn));
            
            smat(abs(obj.xn)>=thresh)=obj.xn((abs(obj.xn)>=thresh));
            
            if nargout>=1   % ||W-F(W)||/||W||
                BSE=norm(sC(sC<thresh),obj.lp_norm);   % needs original values
            end
            if nargout>=2
                % sparsity s should be ca. numel(obj.C)/ct;
                s=sum(sC<thresh);  % needs non.neg. values
            end
            
        end
        
        function nd=signal2nodes(obj)
            % nd=~() find all well defined nodes of signal
            nd=SampledNodes(find(~isnan(obj.xn)),obj.N,obj);
        end
        
        function [nodes,mask]=nodesPDF(obj,c, params)
            % nodes=~(c) nodes sampled with compression rate c
            obj.requireOnly(c>1,'local','compression rate in ]0,1[');
            obj.requireOnly(~obj.isemptydata,'local','non-empty signal');
            if ~exist('params','var')
                params=struct;
            end
            if ~isfield(params,'tol')
                % 5% error
                params.tol=0.05;
            end
            DN = obj.N; % data Size
            p = 1/c; % undersampling factor
            % uniform density
            pdf= p*ones(DN);
            assert(abs(sum(pdf(:))-p*obj.numel)<1,'definition of pdf');
            % generate sampling mask:
            K = sum(pdf(:));  % number of sampled nodes (positive outcomes)
            ctol= max(1,K*params.tol); % max. absolute deviation
            mask = zeros(size(pdf));
            itmax=10;
            it=0;
            while abs(sum(mask(:)) - K) > ctol && it<=itmax
                % Bernoulli trials with success prob. given by pdf
                mask = rand(size(pdf))<pdf;
                it=it+1;
            end
            % compute nodes:
            nodes=SampledNodes(find(mask),obj.N,obj);
        end
        
        function r=PSNR(obj,yn)
            % r=~(yn) ... peak signal to noise ratio of reconstructed signal yn in dB
            % obj.xn is assumed to be the ideal signal
            % used to test quality of similarity of yn to signal xn
            
            function MSE= fitImages(img1)
                % mean square error
                pvxn=max(img1(:))-min(img1(:));
                pvyn=max(yn(:))-min(yn(:));
                if pvxn<pvyn % avoid fitting an image to a constant image, yielding MSE=0
                    p = polyfit(img1(:),yn(:),1);
                    MSE=1/numel(yn)*sum((yn(:)-p(1)*img1(:)-p(2)).^2);
                else
                    p = polyfit(yn(:),img1(:),1);
                    MSE=1/numel(yn)*sum((img1(:)-p(1)*yn(:)-p(2)).^2);
                end
            end
            
            obj.requireOnly(numel(yn)==numel(obj.xn),'local','compatible test signal');
            if isa(yn,'SignalClass')
                yn=yn.xn;
            end
            % mean square error
            % minimized by linear rescaling
            
            if obj.borderSize>0
                ns=obj.size-2*obj.borderSize;
                yn=wkeep(yn,ns);
                xnc=wkeep(obj.xn,ns);
                MSE=fitImages(xnc);
            else
                MSE=fitImages(obj.xn);
            end
            pv=obj.pv; % take always pv from same signal
            r=20*log10(pv/sqrt(MSE));
        end
        
        function ts=PSNRquantiles(obj,yn)
            % ts=~(yn) ... PSNR varying over quantiles of differences
            % obj.xn is assumed to be the ideal signal
            % addresses question: how does PSNR change when only part
            % of the differences are considered?
            
            if isa(yn,'SignalClass')
                yn=yn.xn;
            end
            delta=sort(abs(obj.xn(:)-yn(:)),'ascend'); % sorted differences
            L= numel(delta);
            pv=obj.pv; % peak to valley of original signal
            quantiles=1:100; % 1:100 are percentiles
            psnr_vals=zeros(size(quantiles));
            for j=quantiles
                del1=delta(1:ceil(L*j/100));
                MSE=1/numel(del1)*sum(del1.^2);  % mean square error without fit
                psnr_vals(j)=20*log10(pv/sqrt(MSE));
            end
            ts=TimeSignal(psnr_vals,L,1);
            ts.set_origin(-1/max(quantiles));
            ts.set_scale(PosScales.quantile);
            ts.units='[of ordered differences]';
            ts.amp_quantity='PSNR';
            ts.signalname=...
                ['PSNR vs. quantiles considered (original: ',...
                obj.get_signalname,')'];
            
        end
        
        function r=SSIM(obj,yn,K,window,L)
            % Structural similarity index of comparing yn with the present image obj.xn
            % obj.xn is assumed to be the ideal signal
            r=[];  % only valid for 2d signals
        end
        
        function xdiff=diff(obj, signal2)
            % xd=~(signal2) difference with signal2 (mean value corrected)
            xdiff=obj.clone;
            xdiff.xn=obj.xn-signal2.xn;
            mv_correction=xdiff.mean;
            xdiff.xn=xdiff.xn-mv_correction;
            
            % 3.1g does not work with num2tex
            xdiff.signalname=['difference signal (\mu=',num2tex(xdiff.meanAbs,'%3.1e','none'),...
                ', \sigma=',num2tex(xdiff.stdAbs,'%3.1e','none'),')'];
        end
        
        function pn=norm(obj,p)
            % lp-norm of signal as vector (p=0 admitted)
            if nargin <2
                p=obj.lp_norm;
            end
            if ~isequal(p,0)
                pn=norm(obj.xn(:),p);
            else
                pn=sum(abs(obj.xn(:))>3*obj.sigEst);
            end
        end
        
        function d=dist(obj,s2,p)
            % d=~(s2) lp-distance of present signal to signal s2
            obj.requireOnly(isequal(class(obj),class(s2)),'local', 'same data type');
            if nargin <3
                p=obj.lp_norm;
            end
            d=norm(obj.xn(:)-s2.xn(:),p);
        end
        
        function m=mean(obj)
            % m=~() mean value of signal
            m=meanNaN(obj.xn(:));
        end
        
        function m=min(obj)
            m=min(obj.xn(:));
        end
        
        function m=max(obj)
            m=max(obj.xn(:));
        end
        
        function m=meanAbs(obj)
            % m=~() mean value of signal
            m=meanNaN(abs(obj.xn(:)));
        end
        
        function m=std(obj)
            % m=~() mean value of signal
            m=stdNaN(obj.xn(:));
        end
        
        function m=stdAbs(obj)
            % m=~() mean value of signal
            m=stdNaN(abs(obj.xn(:)));
        end
        
        function s= pv(obj)
            % s= peak to valley of signal
            s=max(obj.xn(:))-min(obj.xn(:));
        end
        
        function s= pv_density(obj)
            % s= peak to valley of signal using only density peaks
            r=findDensityPeaks(obj.xn(:));
            s=max(r.peakx)-min(r.peakx);
        end
        
        function m= stat(obj)
            % m= ~() basic statistical info about signal
            m=struct;
            m.min=obj.min;
            m.mean=obj.mean;
            m.max=obj.max;
            m.pv=obj.pv;
            [m.outliers1, m.outliers2]=obj.outliers(0.001);
        end
        
        function [thresh1, thresh2]=outliers(obj,q, nstd)
            % ~(q) bracket maximally q-fraction of outliers
            % outlier condition is lying outside of mu +/- nstd*sig
            obj.requireOnly(q>0 && q<1,'local','0<q<1');
            if nargin <3
                nstd=3;
            end
            sig=stdNaN(obj.xn(:));
            mu=meanNaN(obj.xn(:));
            M=sort(obj.xn(~isnan(obj.xn)));
            ml=length(M);
            qq=1-0.5*q;
            
            thresh2=max(mu+nstd*sig,M(min(ml,floor(ml*qq ))));
            thresh1=min(mu-nstd*sig,M(max(1,ceil(ml*0.5*q))));
            
        end
        
        function [s2,params]= denoise(obj, mra,params)
            % s2=~([mra,params]) denoise signal; must be redefined!
            % mra is of class MultiResTransform, e.g. mra=Wavelet2D_mlab()
            % example for 2d signal: s2=obj.denoise(Wavelet2D_mlab());
            obj.requireOnly(isa(mra,'MultiResTransform'),'local',...
                'mra is a subclass of MultiResTransform, e.g. Wavelet2D_mlab');
           if ~exist('params','var') || ~isstruct(params)
               params=struct;
           end
           if ~isfield(params,'hard')
               params.hard=true;
           end
           if ~isfield(params,'thresh1')
               params.thresh1=[];
           end
           if ~isfield(params,'thresh2')
               params.thresh2=[];
           end
           if ~isfield(params,'numberOfThresholds')
               params.numberOfThresholds=1;
           end
           if ~isfield(params,'basisname')
               if ~isa(mra,'Curvelet2_clab')
                   params.basisname='db3';
               else
                   params.basisname='dctg2';
               end
           end
                
           mra.set_signal(obj);
           
           try
               mra.set_basisname(params.basisname);
           catch
               warning('basisname not admitted. Using default.');
           end
                      
           mra.dec;
           thresh1=[];
           thresh2=[];
           if params.hard && params.numberOfThresholds==1
               [s2, thresh1]=mra.hardthresh1(params.thresh1);
           elseif params.hard && params.numberOfThresholds>1
               [s2, thresh1,thresh2]=mra.hardthresh2(params.thresh1,params.thresh2);
           elseif ~params.hard
               [s2, thresh1]=mra.softthresh1(params.thresh1);
           end
           s2.colormap_active=obj.colormap_active;
           s2.fig=obj.fig;
           params.thresh1_used=thresh1;
           params.thresh2_used=thresh2;
           
            
        end
        
    end
    
    methods (Hidden)
        
        function ss=N(obj)
            % sample size with padding included
            ss=size(obj.xn);
        end
        
    end
    
    
    %% factories
    
    methods
        
        function s2= make_delta(obj,pos,hN,sig)
            % constructor s2=~([pos,hN,sig]]), delta at pos in matrix of
            % size hN with noise std sig;
            % used to compute representation matrices of linear maps.
            
            if nargin <4 || isempty(sig)
                sig=0;
            end
            if nargin <3 || isempty(hN)
                hN=256;
            end
            if length(hN)~=length(obj.origin)  % obj.size has lenght 2 if empty
                hN=hN(1)*ones(1,length(obj.origin));
            end
            if nargin <2 || isempty(pos)
                pos= floor(hN/2);
            end
            if length(pos)>1
                parg= num2cell(pos);
                pos=sub2ind(hN,parg{:}); % convert to linear index
            end
            
            hR=2*ones(1,length(hN));horigin=hR/2; %  => D=[-Rs/2,Rs/2]
            xdata=zeros(hN);
            xdata(pos)=1;
            if sig>0
                xdata=xdata+sig*randn(size(xdata));
            end
            
            s2=obj.make();
            s2.set_signal(xdata,hR);
            
            s2.signalname=['\delta(x,y)',s2.noise_str(sig)];
            s2.set_scale(PosScales.x);
            s2.origin=horigin;
        end
        
    end
    
    methods (Static)
        
        
        function ts=make()
            % contructor with universal name (callable by parent)
            ts=SignalClass();
        end
        
        function s2= make_standardExample()
            % standard example to be redefined in sublasses
            s2=SignalClass();
        end
        
        function ns=noise_str(sig)
            % text output for noisy signals
            if sig==0
                ns=[];
            else
                ns=['+\epsilon_\sigma, \sigma=',...
                    num2tex(sig,'%3.2e','none')];
            end
        end
    end
    
    %% graphics
    methods
        
        function graph_kde(obj,open_new)
            %  ~([open_new]) plot kerndel density estimator of
            %  distribution of signal values
            if nargin <2 || open_new
                prepfigure(obj.fig);
            end
            
            if license('test','Statistics_Toolbox')
                ksdensity(obj.xn(:));
            else
                params.open_new=false;
                kde(obj.xn(:),params);
            end
            title(['kde of ',obj.get_signalname],'fontsize',12);
            ylabel('density');
            ystr=obj.amp_quantity;
            if isempty(ystr)
                ystr='intensity of signal';
            end
            xlabel(ystr);
        end
        
        function graph_signalWindowed(obj,open_new)
            % ~([open_new]) show windowed signal in the open window by default
            obj.requireOnly(~isempty(obj.winfun),'local','winfun is set');
            if nargin <2 || open_new
                prepfigure(obj.fig);
            end
            
            ts2=obj.apply_winfun();
            ts2.graph_signal(false);
            
        end
        
        function graph_signalInverseWindowed(obj,open_new)
            % ~([open_new]) show windowed signal in the open window by default
            obj.requireOnly(~isempty(obj.winfun),'local','winfun is set');
            if nargin <2 || open_new
                prepfigure(obj.fig);
            end
            
            ts2=obj.apply_winfunInverse();
            ts2.graph_signal(false);
            
        end
        
        function graph_sparsityDefect(obj,open_new,c)
            % ~() plot sparsity defect (lp-error of best s-term approx) of
            % signal
            if nargin <2 || open_new
                prepfigure(obj.fig);
            end
            if nargin<3
                c=[];
            end
            p=obj.lp_norm;
            sd1=obj.sparsityDefect()/numel(obj.xn);
            nvals=numel(sd1);
            semilogx((1:nvals)/nvals,sd1);
            
            if ~isempty(c)
                line([1,1]/c,ylim,'color','r');
                text(1/c,mean(ylim),['c=',num2str(c,'%3.1f')],...
                    'HorizontalAlignment','left');
                idx_sdef=find((1:nvals)/nvals>1/c,1,'first');
                sdef=sd1(idx_sdef);
                line(xlim,[sdef,sdef],'color','g');
                text(1/nvals,sdef,['\sigma_{1/c}/dim(x)=',num2str(sdef,'%3.1e')],...
                    'HorizontalAlignment','left','VerticalAlignment','bottom');
                
            end
            hylim=ylim;
            hylim(2)=1.01;
            ylim(hylim);
            xlabel('rel. sparsity s/dim(signal)');
            ylabel(['\sigma_s(x)_{',num2str(p),'}_{',num2str(p),'}/dim(x)']);
            titstr{1}=['l_{',num2str(p),'}-error of best s-term approx.'];
            titstr{2}=obj.get_signalname;
            title(titstr,'fontsize',12);
        end
        
    end
    
end


