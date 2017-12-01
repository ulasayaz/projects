classdef Wavelet1D < MultiResTransform1
    % 1-dim. decimated wavelet transformations with graphical output
    % (analogous to FourierTrafos1D)
    % features:
    %    --- both discrete multi-resolution analysis (.dec)
    %    and continuous (.cwt)
    %    --- noise filter (.graph_noisefilter)
    %    
    % Example: 
    %====================================
    % --- discrete wavelet trafo (multi-resolution analysis):
    %   w= Wavelet1D(); w.show_test;
    %
    %   signal=TimeSignal.make_chebyshev(3);
    %   w= Wavelet1D(signal);
    %   w.show_trafos_selection();
    % --- continuous wavelet trafo:
    %   w=Wavelet1D.test_cwt();
    %
    % omega0=5;sig=0.2;hN=512;
    % signal=TimeSignal.make_sinus(omega0,sig,hN);
    % w= Wavelet1D(signal);
    % --- discrete wavelet trafo (multi-resolution analysis):
    % w.dec;
    % w.show_resultd;
    % --- continuous wavelet trafo:
    % w.set_scales(1:40);
    % w.cwt;
    % w.show_resultc;
    %
    % --- rectangle signal (Gibbs for FT)
    % sig=0.05;
    % w= Wavelet1D(TimeSignal.make_rect(30,sig,128));
    % w.graph_compare_compression;
    % b=w.bestwavelet(w.ts);disp(b);
    % w.set_basisname(b);
    % w.set_dwtmode('zpd'); % 0-padding at border
    % w.dec;
    % w.show_noisefilter;
    % -- Test with Chebyshevs:
    % Wavelet1D.test_bestwavelet(0.99);
    %
    %
    % GT Berlin, 2013
    %
    
    properties (Access=public)
        
        L    %@<vector> wavelet decomposition:  bookkeeping vector
        CWTcoeffs %@<matrix> coeffs of CWT
      
    end
    
    properties (SetAccess=private)
        
        scales  %@ (vector) scale vector for CWT (cont. wavelet trafo)
        
    end
    
    
    %% constructor and commands
    methods
        
        function obj=Wavelet1D(signal)
            % constructor
            if nargin==0
                signal=TimeSignal();
            end
            obj = obj@MultiResTransform1(signal);
            obj.requireOnly(license('test','Wavelet_Toolbox'),'local','needs wavelet toolbox');
            obj.set_basisname('db1');
            %obj.dwtmode_active=dwtmode('status','nodisp');  % default border extension
            % default border extension: sp1 (smooth padding) works well for
            % smooth signals, 'sym' works well for photos.
            % photos:
            obj.dwtmode_active=dwtmode('sym','nodisp');
            obj.scales=[];
            
        end
        
        
        function set_dwtmode(obj, modestr)
            % ~(modestr) border extension mode cf. help dwtmode
            % e.g. 'sym' (symmetric extension), 'zpd' (zero-padding), etc.
            obj.dwtmode_active=modestr;
            dwtmode(modestr,'nodisp');
        end
        
        function set_algo(obj)
            set_algo@MultiResTransform1(obj);
            obj.algo.name='DWT'; % decimated wavelet trafo
            obj.algo.toolbox='Wavelet\_Toolbox';
            
        end
        
        function reset_Trafo(obj)
            % reset trafo (e.g. after signal has changed)
            obj.C=[]; obj.L=[];
        end                
        
        function set_signalfun(obj,xt,hN,hT)
            % ~(xt,hN,[hT]) set a signal function and sample it at hN points
            obj.requireOnly( isa(xt,'function_handle')...
                && obj.isnatural(hN),'local',...
                'input is function plus sampling length');
            if nargin<3
                hN=100;
            end
            if nargin <4
                hT=2;
            end
            obj.ts=TimeSignal(xt,hN,hT);
            obj.reset_Trafo;
            obj.ensureOnly(~obj.isTrafodone,'local', ...
                'wavelet decomposition is reset');
        end
        
        
        function interp1(obj, upsamp, ipmethod)
            % ~(upsm, [ipm]) interpolate by upsampling by a factor of upsm
            % using method ipm (default is 'spline').
            %
            obj.requireOnly(obj.isnatural(upsamp),'local',...
                'upsamp is a natural number');
            obj.requireOnly(nargin<3 || ischar(ipmethod),'local',...
                'ipmethod is a string denoting a method');
            if nargin < 3 || isempty(ipmethod)
                ipmethod='spline';
            end
            obj.ts.interp1(upsamp, ipmethod);
            
            obj.reset_Trafo;
            obj.ensureOnly(~obj.isTrafodone,'local', ...
                'wavelet decomposition is reset');
        end
        
        function set_scales(obj,hs)
            % ~(hs) set scale vector for CWT
            obj.requireOnly(isvector(hs),'local',' is scale vector');
            obj.scales=hs;
            obj.CWTcoeffs=[];
            obj.ensureOnly(~obj.isCWTdone,'local','cwt reset');
        end                
        
        function dec(obj)
            % discrete multiresolution analysis
            obj.require(obj.N>0,'local','non-empty sample size');
            [obj.C,obj.L] = wavedec(obj.xn,obj.deepestlev(),obj.basisname);
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wc)
            %yn=~(cC) % signal reconstruction from wavelet decomp wc
            % wavelet recon from matlab's wavelet toolbox
            obj.requireOnly(~isempty(obj.L),'local','needs recon info');
            yn = waverec(wc,obj.L,obj.basisname) ;
        end
        
        function cwt(obj)
            % continuous wavelet trafo
            obj.require(obj.N>0,'local','non-empty sample size');
            obj.CWTcoeffs=cwt(obj.xn,obj.get_scales, obj.basisname);
            obj.ensureOnly(obj.isCWTdone,'local', 'wavelet trafo done');
        end
    end
    %% queries
    
    methods
        
        function sL=length(obj)
            sL=length(obj.ts.xn);
        end
        
        function sampleSize=size(obj)
            % s=~(): sample size with padding inmcluded
            sampleSize=size(obj.ts.xn);
        end
                
        function L=level(obj)
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            L= size(obj.L,1)-2;
        end        
        
        function [phi,psi,xval] = wavefun(obj,iter)
            % [phi,psi,xval]=~(iter) returns scaling and wavelet function
            [phi,psi,xval] = wavefun(obj.basisname,iter);
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x
            yn=wavedec(x,obj.deepestlev(),obj.basisname);
        end
        
        function xn= synthesize(obj,y)
            % xn=~(y) synthesize
            obj.requireOnly(~isempty(obj.L),'local','needs recon info');
            xn=waverec(y,obj.L,obj.basisname) ;
        end
        
        function coeff=detcoef(obj,level)
            % coeff=~(level) extract detail coefficients from wavelet transform
            % type ='d' for diagonal fluctutation.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(level>=1 && level<=obj.level,'local',...
                'level is below computed level');
            coeff=detcoef(obj.C,obj.L,level);
        end
        
        function A= appcoef(obj)
            % coeff=~() extract deepest level approximate coefficients
            % from wavelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            A= appcoef(obj.C,obj.L,obj.basisname,obj.level);
        end
        
        function ok=isTrafodone(obj)
            % ok=~(): is result of DWT available?
            ok=~isempty(obj.C) && ~isempty(obj.L);
        end
        
        function ok=isCWTdone(obj)
            % ok=~(): is result of CWT available?
            ok=~isempty(obj.CWTcoeffs);
        end
                        
        function s= get_scales(obj)
            % s=~(), get scale vector from property or calculate one
            if ~isempty(obj.scales)
                s=obj.scales;
            else
                s=1:min(obj.N, max(20,1:ceil(obj.N/10)));
            end
        end
        
    end
    
    methods (Hidden)
        
        
        function sL=N(obj)
            % sample length with padding included
            sL=length(obj.ts.xn);
        end
        
        function sTL=Ts(obj)
            % sampling time interval
            sTL=obj.ts.Ts;
        end
        
        function v=xn(obj)
            % signal vector
            v=obj.ts.xn;
        end
        
        function v=padd(obj)
            % zero padding length
            v=obj.ts.padd;
        end
        
        
    end
    
    
    
    methods
        %% coordinate transformations
        
        
        
        %% filters
        
        
        
        function [yn, thresh]=noisefilter_shiftinvariant(obj, thresh)
            % [yn,thresh,sd]=~(thresh,[as]) denoise using shift-inv. wavelet thresholding
            % apply e.g. to a test signal created by w= TimeSignal.make_sinus(5,0.1).
            % use ndwt and indwt (undecimated WT) to ensure shift invariance
            obj.require(obj.isTrafodone,'local',obj.msgDoDWT);
            
            if nargin < 2 || isempty(thresh)
                obj.estimateNoise();
                thresh=4.5*obj.sigma;
            end
            
            hC=ndwt(obj.xn, obj.level, obj.basisname);
            % threshold all detail signals
            for j=1:obj.level
                d=hC.dec{j+1};
                d(abs(d)<thresh)=0;
                hC.dec{j+1}=d;
            end
            % reconstruct
            yn=indwt(hC);
            
        end
        
        
        
    end
    
    
    %% static
    methods (Static)
        
        function str=msgDoDWT()
            str='Wavelet decomposition available (call wavedec).';
        end
        
        function f= make(ts)
            % make an object from a time signal
            assert(isa(ts,'TimeSignal'),'ts is time signal');
            f=Wavelet1D();
            f.set_ts(ts);
        end
        
        %% test
        function w=test_dwt(xdata,hN)
            % testing the wavelet trafo of this class
            if nargin ==0 || isempty(xdata)
                xdata=@(t) sin(5*t);
            end
            if nargin <2
                hN=128;
            end
            
            %  choose del0 ~=0 to construct a signal with x[N]~=x[0]
            %  del0=0; del0=0.1;
            del0=pi/4;
            hT=2*pi+del0;
            s=TimeSignal(xdata,hN,hT);
            w=Wavelet1D(s);
            w.fig=2;  % max. 2 open windows
            
            w.show_dwt_trafos({'db1','db2','db3'});
            
        end
        
        function w=test_cwt(xdata, hN, hT, wvname)
            % testing the wavelet trafo of this class
            % with a delta and a simple harmonics as signal
            if nargin ==0 || isempty(xdata)
                if nargin <2 || isempty(hN)
                    hN=256;
                end
                % delta function
                xdata=zeros(hN,2);
                xdata(floor(hN/2),1)=1;
                % simple harmonics
                del0=pi/4;
                hT=2*pi+del0;
                nv=linspace(1,hT,hN);
                xf=@(t) sin(5.*t);
                xdata(:,2)= xf(nv(:));
            else
                xdata=xdata(:);
            end
            if nargin <3
                wvname='haar';
            end
            
            cc=size(xdata,2)+1;
            w=Wavelet1D();
            w.set_basisname(wvname);
            w.fig=2;
            
            prepfigure(w.fig);
            subplot(2,cc,1);
            w.graph_wavefun([],false);
            
            for j=1:cc-1
                w.set_signal(TimeSignal(xdata(:,j)));
                w.scales=1:min(w.N,floor(min(200,w.N/2)));
                
                subplot(2,cc,j+1);
                w.graph_signal(false);
                
                subplot(2,cc,cc+j+1);
                w.cwt;
                w.graph_cwt(false);
                title('CWT (max at max slope of signal!)');
            end
        end
        
        
        function w=test1_distortion(omega)
            % test noisefilter deformation of a simple harmonics
            assert(nargin <1 || isnumeric(omega),'omega is angular freq');
            
            % testing the noisefilterfilter
            if nargin ==0 || isempty(omega)
                omega=5;
            end
            
            prepfigure(2);
            
            xdata=@(t) sin(omega*t);
            hN=128;
            del0=[0,0.1,pi/4];
            s=TimeSignal(xdata,hN,2*pi);
            w=Wavelet1D(s);
            minxn=min(w.xn);
            maxxn=max(w.xn);
            wxn=maxxn-minxn;
            mxn=0.5*(minxn+maxxn);
            common_ylim=[mxn-0.6*wxn, mxn+0.6*wxn];
            
            for j=1:3
                subplot(3,3,j);
                %  choose del0 ~=0 to construct a signal with x[N]~=x[0]
                %  del0=0; del0=0.1; del0=pi/4;
                hT=2*pi+del0(j);
                s=TimeSignal(xdata,hN,hT);
                w=Wavelet1D(s);
                w.wavedec;
                w.estimateNoise();
                thresh=4.5*w.sigma;
                
                yn=w.noisefilter_simple(thresh);
                err_max=max(abs(yn-w.xn));
                err_med=medianNaN(abs(yn-w.xn));
                
                plot(1:w.N,w.xn,'--k', 1:w.N, yn,':r','Linewidth',2);
                title({['deviation from full period \delta=',...
                    num2str(del0(j),'%3.1f')]...
                    },...
                    'fontsize',12);
                legend('Original','filtered','location','best');
                ylim(common_ylim);
                
                subplot(3,3,j+3);
                plot(1:w.N,yn-w.xn);
                title({['error[max, median])= ',...
                    vec2str([err_max, err_med],'%3.1e')]...
                    }, 'fontsize',12);
                
                subplot(3,3,j+6);
                plot(w.C);
                xlabel('n');
                title({['Wavelet Trafo,' w.basisname,],['thresh=',...
                    num2str(thresh)]},...
                    'fontsize',12);
            end
            suptitle(['Rauschfilterverzerrung einer einfachen Harmonischen - ',...
                'Method : noisefilter\_simple'], 12);
        end
        
        function test_bestwavelet(enthresh)
            % what happens when signal gets smoother (more often
            % differentiable)
            % Test with Chebyshev polynomials
            if nargin <1
                enthresh=0.99;
            end
            c=5;
            w=Wavelet1D();
            hN=512;
            step=15; % 3
            prepfigure(w.fig);
            for j=1:c
                subplot(2,c,j);
                hts=TimeSignal.make_chebyshev(step*(j-1)+1,1,0,hN);
                w.set_signal(hts);
                w.graph_signal(false);
                subplot(2,c,j+c);
                w.show_compare_compression(enthresh);
            end
            suptitle('increasing order of polynomial increases order of best wavelet',14);
            
        end
        
    end
    
    %% graphics
    methods
        
        function graph_wavefun(obj, iter, open_new)
            % plot wavelet in open window (default)
            if nargin <2 || isempty(iter)
                iter=20;
            end
            if nargin <3 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            [~,psi,xval] = obj.wavefun(iter);
            plot(xval,psi); %axis([0 1 -1.5 1.5]);
            title(['Wavelet ',obj.basisname, ', it=',num2str(iter)],...
                'fontsize',obj.fontsize);
        end
        
        function graph_cwt(obj,open_new)
            % show result using the continuous Wavelet trafo in open win
            obj.require(obj.isCWTdone,'local','Wavelet decomposition available (call cwt).');
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            
            scales_used=obj.get_scales;
            imagesc(1:obj.N,scales_used,obj.CWTcoeffs);
            axis xy;
            xlabel('n');
            ylabel('scales');
            colormap jet; colorbar('location','east');
            title(['Continuous Wavelet Trafo ', obj.basisname],...
                'fontsize', 12);
            
        end
        
        function show_resultc(obj)
            % show result using the continuous Wavelet trafo in a new window
            obj.require(obj.isCWTdone,'local','Wavelet decomposition available (call cwt).');
            obj.require(obj.fig~=0,'local', 'figure output activated');
            
            prepfigure(obj.fig);
            subplot(2,1,1);
            obj.graph_signal(false);
            subplot(2,1,2);
            obj.graph_cwt(false);
            
        end
        
        function graph_compare_compression(obj, open_new, enthresh)
            % compare compression capabilty of Daubechies
            if nargin <3
                enthresh=0.99;
            end
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            nmax=20;
            cfactor=ones(1,nmax);
            hN=obj.N;
            for j=1:nmax
                obj.set_basisname(['db',num2str(j)]);
                obj.dec;
                ens=sort(obj.C.^2,'descend');
                ennorm=sum(ens);
                enmap=cumsum(ens)/ennorm;
                idx99=find(enmap>=enthresh,1,'first');
                cfactor(j)=hN/idx99;
            end
            [bestn,bestidx]=max(cfactor);
            plot(1:nmax,cfactor,bestidx,bestn,'or');
            text(bestidx,bestn,num2str(bestidx),'VerticalAlignment','bottom');
            legend('#','best','location','best');
            legend('#','best','location','best');
            title({'Compression',obj.ts.signalname},'fontsize',12);
            xlabel('Daubechies dbn');
            ylabel(['compression factor at ',num2str(100*enthresh,'%3.2f'),'% point']);
            
        end
        
        function show_quality_compression(obj, p, wvnames)
            % compare PSNR and RBSE as a function of the compression factor
            % PSNR ... peak signal to deviation ratio of reconstructed signal
            % RBSE ... rel. error of best s-term approximation
            if nargin <3
                wvnames={'db1','db2','db3','db4'};
            end
            if nargin <2
                p=1;
            end
            p_str=num2str(p);
            jmax=20;
            fontsize=12;
            
            L=length(wvnames);
            clevel=0.99;
            b= obj.bestwavelet(obj.ts);
            if ~ismember(b,wvnames)
                wvnames{L+1}=b;
            end
            PSNR=zeros(jmax,length(wvnames));
            PSNR(1,:)=NaN;
            RBSE=PSNR;
            cfactor=1:jmax;
            
            for k=1:length(wvnames)
                obj.set_basisname(wvnames{k});
                obj.dec;
                for j=2:jmax
                    [yn,RBSE(j,k)]=obj.sparseApprox(cfactor(j),p);
                    PSNR(j,k)=obj.ts.PSNR(yn);
                end
            end
            
            prepfigure(obj.fig);
            
            subplot(2,2,1);
            obj.graph_signal(false);
            
            subplot(2,2,2);
            semilogy(cfactor, numel(obj.ts.xn)./cfactor);
            title('sparsity','fontsize', fontsize);
            xlabel('compression factor');
            ylabel('s');
            
            subplot(2,2,3);
            plot(cfactor,PSNR);
            hleg=legend(wvnames,'location','best');
            htitle = get(hleg,'Title');
            set(htitle,'String','Wavelets')
            xlabel('compression factor');
            ylabel('PSNR [dB]');
            title({'PSNR (peak signal to deviation ratio of reconstructed signal)',...
                ['(Wavelet filters cf. legend)']},'fontsize',fontsize);
            
            subplot(2,2,4);
            plot(cfactor,RBSE);
            hleg=legend(wvnames,'location','best');
            htitle = get(hleg,'Title');
            set(htitle,'String','Wavelets')
            xlabel('compression factor');
            ylabel(['RBSE=\sigma_s(W)_{',p_str,'}/||W||_',p_str]);
            title({['RBSE (rel. best s-term approx. error in l_{',p_str,'}-norm)'],...
                ['is in l_1-norm a bound of reconstruction error in compressive sensing']},...
                'fontsize',fontsize);
            
        end
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 1d';
            ok= isa(obj.ts,'TimeSignal');
        end
    end
    
end




