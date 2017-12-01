classdef SparseRep< DC
    % sparse representations using Wavelets
    %
    % Example
    % 1.) create signal:
    % --- small image: signal=Signal2D.make_fromImage('cameraman.bmp');
    % --- large image: signal=Signal2D.make_fromImage('parrots.tiff'); signal.crop();
    % --- large image: signal=Signal2D.make_fromImage('Houses8.tiff'); signal.crop();
    % --- large image: signal=Signal2D.make_fromImage('river13.tiff'); signal.crop();
    % 2.) create a sparisfying transform:
    % --- matlab tolbox: w= Wavelet2D_mlab(signal); w= Wavelet2WP_mlab(signal);
    % --- wavelab toolbox: WavePath; % if not yet done
    %                      w= Wavelet2D_wlab(signal);    
    % 3.) create main object:
    %     sr=SparseRep(w);
    %     sr.figopt.pixsizeX=1000;
    %   a) optional: plot phase transition curves: sr.graph_PhaseTransitionDonoho;
    %                                              sr.graph_SufficientCP_Rauhut;
    % 3.) compute quality bounds for selected signal:
    % --- might need to add path for curvlets: addpathSplittingSolvers(false); 
    %     sr.show_quality_frameSet({1,2,{0,0,1,2}});
    % --- results for periodic boundary conditions are worse than for
    %     symmetric ones:
    %     sr=SparseRep(Curvelet2_clab());
    % 4.) variant: use all images in all subdirs
    %     sr.set_imgfolder([],true);  % 
    %     sr.show_quality_images({1,2,{0,1,2}});
    %     sr.show_Qdistribution();
    %
    
    properties (Access=public)
        frameSet %@(cellstr) wavelet/curvelet names to be studied
        w       %@(MultiResTransform2) determines type set by frameSet
        cspecial  %@(vector) special compression ratios of interest
        PSNR    %@(matrix) distribution of peak signal to noise ratios
        SSIM    %@(matrix) structural similarity index between original and recon
        BSE     %@(matrix) error of best s-term approximation in p-norm
        imgfiles %@(cellstr) names of image files to be studied
        
        algo    %@ (struct) finger print of algorithm
        fig     %@ (int) number of open figures
        figopt %@ (int,default empty struct) figure options
        fontsize %@(int) fontsize of titles of subplots
        verbose %@(boolean)
        do_ssim %@(boolean) compute also quality parameter SSIM
        maxThumnails %@(int) max. number of images shown as matrix
        default_wvnames
        toolbox_pref %@(char) preeferred toolbox
        
        pretransformFFT  %@(bolean
        
    end
    
    
    
    
    %% commands and constructors
    methods
        
        function obj=SparseRep(w)
            % obj=~(w) construct from wavelet w
            obj.requireOnly(isempty(w) || isa(w,'MultiResTransform2')|| ...
                isa(w,'FourierTrafos2') || isa(w,'Fourier2') ,...
                'local','w is of type MultiResTransform2');
            if nargin <1
                w=[];
            end
            obj.w=w;
            obj.fig=1;
            obj.set_algo;
            obj.verbose=true;
            obj.do_ssim=true;
            obj.fontsize=12;
            obj.figopt=struct;
            obj.maxThumnails=12;  % 3x4
            obj.cspecial=[8,16];  % 8x and 16x compression
            obj.default_wvnames={'db1','db2','db3','bior4.4','dctg2','bicubic'}; % Daubechies, biorthogonal, curvelet
            obj.frameSet=obj.default_wvnames;
            obj.set_imgfolder(fullfile(get_projectpath(),'testdata','images')); % image files
            if license('test','Wavelet_Toolbox')
                obj.toolbox_pref='Wavelet_Toolbox';
            else
                obj.toolbox_pref='wavelab';
            end
            obj.pretransformFFT=false;
        end
        
        function set_algo(obj)
            obj.algo.version=0.91;
            obj.algo.versiondate='13.01.2014';
            obj.algo.variant=[];
            obj.algo.name='';
            obj.algo.mfilename=class(obj);
            obj.algo.author='G. Troll';
        end
        
        function set_imgfolder(obj, imgfolder, isrecursive)
            % ~(imgfolder,isrec) creates full filenames of all non-empty files in imgfolder
            % imgfolder need not exist. In this case the current directory
            % is used;
            if nargin <2 || isempty(imgfolder)
                imgfolder=fullfile(get_projectpath(),'testdata','images');
                if ~exist(imgfolder,'dir')
                    imgfolder=cd();
                end
            end
            if nargin <3
                isrecursive=false;
            end
            
            if isrecursive
                obj.imgfiles=dirrec(imgfolder);
            else
                recfiles=dir(imgfolder);
                recfiles=recfiles([recfiles(:).bytes]>0);  % remoce empty files
                [~,idx]=sort([recfiles(:).bytes],'descend');
                recfiles=recfiles(idx);
                extendfull=@(k) fullfile(imgfolder,recfiles(k).name);
                obj.imgfiles=arrayfun(extendfull,1:length(recfiles),'UniformOutput',false);
            end
        end
        
        function set_w(obj,w)
            % ~(w) set wavelet object (determining type created from frameSet).
            obj.requireOnly(isa(w,'MultiResTransform2'),...
                'local','w is of type MultiResTransform2');
            obj.w=w;
        end
        
        function add_frame(obj,fn)
            % ~(fn) add frame(s) fn to frameSet
            obj.frameSet=union(obj.frameSet,fn);
        end
        
        function remove_frame(obj,fn)
            % ~(fn) remove frame(s) fn to frameSet
            obj.frameSet=setdiff(obj.frameSet,fn);
        end
        
        
    end
    
    
    %% queries
    methods
        
        function tsa=ts(obj)
            % tsa=~() signal (image) used
            tsa=obj.w.ts;
        end
        
        function nm=signalname(obj)
            % nm=~() signal name
            nm=obj.ts.signalname;
        end
        
        function ig=img(obj)
            % I=~() image data (matrix)
            ig=obj.ts.xn;
        end
        
        function wvl=make_wavelet2D(obj)
            % wvl=~() make wavelet choosing availbale toolbox
            if license('test','Wavelet_Toolbox')
                wvl=Wavelet2D_mlab(obj.ts);
            else
                wvl=Wavelet2D_wlab(obj.ts);
            end
        end
        
        function crop_images(obj,maxsize,fmt)
            % ~(fmt) crop to nearest square power of 2, resize to maxsize and
            % and save in format fmt
            
            if nargin <3
                fmt='png';
            end
            if nargin<2
                maxsize=Inf;
            end
            obj.set_imgfolder();
            for j=1:length(obj.imgfiles)
                
                try
                    s2=Signal2D.make_fromImage(obj.imgfiles{j});
                    s2.crop(2.^floor(log2(s2.size)));
                    if size(s2.xn,1)>maxsize
                        s2.resize(maxsize);
                    end
                    
                    [~,fn]=fileparts(obj.imgfiles{j});
                    s2.imwrite([fn,'.',fmt],fmt);
                catch
                    continue;
                end
                
            end
            
        end
        
        
        function smax=s1w_RauhutEq924(obj,c)
            % smax=~(c) max. sparsity for a number of measurements with compression c.
            % in the case of large M: m=2*s_m*ln(e*M/s_m) and nonuniform
            % recovery with prob. p close to 1
            % c_m=M/m; c_p=M/sp;
            % source: [MICS] eq. 9.24
            M=obj.ts.numel;
            m=M/c;  % # measurements
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
        end
        
        function smax=s0s_Rauhut_p293(obj,c)
            % smax=~(c) max. sparsity for a number of measurements with compression c.
            % in the case of delta_2s=0.6129, uniform recovery with prob. p close to 0
            % c_m=M/m; c_p=M/sp;
            % Source [MICS] p.293
            M=obj.ts.numel;
            m=M/c;  % # measurements
            e=exp(1);
            a=27.434;
            fun=@(s) m-2*a*s.*log(e*M./(2*s))-a*log(2);
            % fun assumes min at M; sparsity limit smax is in [1, argmin]
            argmin=M/2;
            fmin=fun(argmin);
            % [argmin,fmin] = fminbnd(fun,1,M);
            if fmin<0
                smax=floor(fzero(fun,[1,argmin]));
            else % no restriction as to sparsity: all vectors can be reconstructed
                smax=M;
            end
        end
        
        
        function smax=s1w_Donoho(obj,c)
            % smax=~(c) max. sparsity for a number of measurements with compression c.
            % weak phase transition for Gaussians according to Donoho et al.
            % in the case of m/M -> 0.
            % recovery with prob. p close to 1 if s < smax.
            % failure with high prob. if s>smax
            % c_m=M/m; c_p=M/sp;
            % source: [MICS] notes p.305
            M=obj.ts.numel;
            m=M/c;  % # measurements
            smax=0.5*m./log(M./m);
        end
        
        function smax=s1s_Donoho(obj,c)
            % strong phase transition for Gaussians according to Donoho et al.
            % in the case of m/M -> 0.
            % recovery with prob. p close to 1 if s < smax.
            % failure with high prob. if s>smax
            % c_m=M/m; c_p=M/sp;
            % source: [MICS] notes p.305
            M=obj.ts.numel;
            m=M/c;  % # measurements
            e=exp(1);
            spi=sqrt(pi);
            smax=m./(2*e*log(M./(spi*m)));
        end
        
    end
    
    %% static
    methods (Static)
        
        function b=iscurvelet(name)
            % b=~(name) is <<name>> a curvelet?
            b=ismember(name,{'dctg2','dctg2fft'});
        end
        
        function b=isfourier(name)
            % b=~(name) is <<name>> a fourier transform?
            b=ismember(name,{'fft2'});
        end
        
        function b=isinterp(name)
            % b=~(name) is <<name>> a fourier transform?
            b=ismember(name,{'bilinear','bicubic'});
        end
        
        function w=make_trafo(wvname, ts, toolbox_pref)
            % ~(w) set wavelet object (determining type created from frameSet).
            assert(ischar(wvname) && isa(ts,'Signal2D'), 'types ok');
            if nargin<3
                toolbox_pref=[];
            end
            if SparseRep.iscurvelet(wvname)
                w=Curvelet2_clab(ts);
            elseif SparseRep.isfourier(wvname)
                ts2=ts.clone;
                ts2.xn=ts2.xn-mean(ts2.xn(:));
                w=Fourier2(ts2);              
            elseif SparseRep.isinterp(wvanem)
                w=IdTransform(ts);               
            else
                if strcmpi(toolbox_pref,'Wavelet_Toolbox')
                    w=Wavelet2D_mlab(ts);
                else
                    w=Wavelet2D_wlab(ts);                
                end
            end
            w.set_basisname(wvname);
            
        end
        
        function sfn=shorten(fullfn,Lmax)
            % sfn=~(fullfn,Lmax) shorten full filename to name of max. Lmax chars
            [~,sfn]=fileparts(fullfn);
            sfn=regexp(sfn,'[a-zA-Z0-9]*','once','match');
            if nargin >1
                sfn= sfn(1:min(Lmax,length(sfn)));
            end
        end
        
    end
    
    
    %% graphics
    methods
        
        function show_sparseApprox(obj, cs, p)
            % ~(c_s,p) show best s-sterm approx. at compression factor c_s
            % in new window
            obj.require(obj.fig~=0,'local', 'figure output activated');
            obj.require(obj.w.isDWTdone,'local',obj.w.msgDoDWT);
            
            if nargin < 2 || isempty(cs)
                cs=16;
            end
            if nargin <3 || isempty(p)
                p=2;  % lp-norm
            end
            p_str=num2str(p);
            
            [s2,obj.BSE]=obj.w.sparseApprox(obj,cs,p);
            obj.PSNR=obj.ts.PSNR(s2.xn);
            lp=obj.ts.norm(p);                        
            
            s3=obj.ts.clone;
            s3.replace_signal(yn-obj.ts.xn);           
            lp_diff=s3.norm(p);
            lp_quot=lp_diff/lp;
            
            method_str=['(',obj.w.basisname,', c_s=',num2str(cs),')'];
            result_str=['||W-F_s(W)||_{p=',p_str,'}=', num2tex(obj.BSE,'%3.1e','none'),...
                ', PSNR=',num2str(obj.PSNR,'%3.0f')];
            
            prepfigure(obj.fig,obj.algo,obj.figopt);
            subplot(1,3,1);
            obj.w.graph_signal(false);
            
            subplot(1,3,2);
            s2.graph_signal(false);
            title({['Best s-term filter W^{-1}\circ F_s\circ W, ',method_str],...
                result_str}, 'fontsize',obj.fontsize);
            
            subplot(1,3,3);
            s3.graph_signal(false);
            res2_str=['||I-J||_{p=',p_str,'}/||I||_p=', num2tex(lp_quot,'%3.1e','none'),...
                ];
            title({['Difference Image ',method_str],res2_str}, 'fontsize',obj.fontsize);
        end
        
        
        
        function graph_SufficientCP_Rauhut(obj,open_new)
            % source: [MICS] thm. 9.27
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            M=obj.w.numel;
            c=linspace(1,20);
            m=M./c;
            
            mspecial=M./obj.cspecial;
            s1w_special=arrayfun(@obj.s1w_RauhutEq924,obj.cspecial);
            s1s_special=arrayfun(@obj.s0s_Rauhut_p293,obj.cspecial);
            
            s1w=arrayfun(@obj.s1w_RauhutEq924,c);
            s1s=arrayfun(@obj.s0s_Rauhut_p293,c);
            
            semilogy(c,M./s1w,c,M./s1s,...
                obj.cspecial,M./s1w_special,'o',...
                obj.cspecial,M./s1s_special,'d'...
                );
            opt.format='%4.1f'; opt.fontsize=12; opt.outputY=true;
            textvec(obj.cspecial,M./s1w_special,opt);
            textvec(obj.cspecial,M./s1s_special,opt);
            ylabel('c_p=M/s_p','fontsize',obj.fontsize);
            
            xlabel('c=M/m','fontsize',obj.fontsize,'fontsize',obj.fontsize);
            title('Strong (UR) and weak (NUR) sufficient conditions on c_p (Rauhut)','fontsize',...
                obj.fontsize);
            hleg=legend('non-uniform c^W_1',' c^S_{p\rightarrow 0}, \delta_{2s}\leq 0.61',['c\in',vec2str(obj.cspecial)],...
                ['c\in',vec2str(obj.cspecial)],...
                'location','best');
            set(hleg,'FontSize',obj.fontsize);
            
        end
        
        function graph_PhaseTransitionDonoho(obj, do_plotrho, open_new)
            % Gaussian polytope transitions
            % source: Donoho [5] Thm. 1.4-5, also cf.[MICS] notes p.305
            % rho:= smax/m vs c:=M/m
            % c_s=2*e*c*log(c/pi^0.5)
            % c_w= 2*c*log(c)
            if nargin <2
                do_plotrho=false;
            end
            if nargin <3 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            M=obj.w.numel;
            cmax=20;
            cmin=4;
            
            c=linspace(1,cmax);
            m=M./c;
            
            mspecial=M./obj.cspecial;
            s1w_special=arrayfun(@obj.s1w_Donoho,obj.cspecial);
            s1s_special=arrayfun(@obj.s1s_Donoho,obj.cspecial);
            
            s1w=arrayfun(@obj.s1w_Donoho,c);
            s1s=arrayfun(@obj.s1s_Donoho,c);
            
            s1w(c<cmin)=NaN;
            s1s(c<cmin)=NaN;
            
            
            if do_plotrho
                plot(c,m./s1w,c,m./s1s,...
                    obj.cspecial,mspecial./s1w_special,'o',...
                    obj.cspecial,mspecial./s1s_special,'d'...
                    );
                opt.format='%4.1f'; opt.fontsize=obj.fontsize; opt.outputY=true;
                textvec(obj.cspecial,mspecial./s1w_special,opt);
                textvec(obj.cspecial,mspecial./s1s_special,opt);
                ylabel('1/\rho=m/s','fontsize',obj.fontsize);
                
            else
                plot(c,M./s1w,c,M./s1s,...
                    obj.cspecial,M./s1w_special,'o',...
                    obj.cspecial,M./s1s_special,'d'...
                    );
                opt.format='%4.1f'; opt.fontsize=obj.fontsize; opt.outputY=true;
                textvec(obj.cspecial,M./s1w_special,opt);
                textvec(obj.cspecial,M./s1s_special,opt);
                ylabel('M/s','fontsize',obj.fontsize);
                
            end
            
            xlim([cmin,cmax]);
            xlabel('c=M/m','units','normalized','Position',[0.95,-0.08],...
                'fontsize',obj.fontsize);
            title({'Strong (UR) and weak (NUR) phase transition curves c_{1/2}=M/s_{1/2}',...
                'in the limit m/M \rightarrow 0 (Donoho)'}, 'fontsize',obj.fontsize);
            hleg=legend('non-uniform','unif. recovery',['c\in',vec2str(obj.cspecial)],...
                ['c\in',vec2str(obj.cspecial)],...
                'location','best');
            set(hleg,'FontSize',obj.fontsize);
            
            
        end
        
        
        function graph_mVScp(obj,ctest, open_new)
            % m-min vs cp=s/M   for compression ctest=M/m;
            % according to Rauhut
            if nargin <2
                ctest=16;
            end
            if nargin <3 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            M=obj.w.numel;
            m=M/ctest;  % # measurements
            e=exp(1);
            spi=sqrt(pi);
            fun1=@(s) (m-2*e*s.*log(M./(spi*s)))/M;
            fun2=@(s) (m-2*s.*log(e*M./s))/M;
            s=1:10:M;
            
            % fun2 has min at s=M
            smax2=M;
            % fun1 has min at s=M/(e*spi)
            smax1=M/(e*spi);
            plot(s/M,fun1(s),s/M,fun2(s),...
                smax1/M,fun1(smax1),'s',...
                smax2/M,fun2(smax2),'s',...
                1,fun1(M),'o', 1, fun2(M),'o');
            line([0,2],[0,0],'color','red');
            xlabel('sparsity compression s/M','fontsize',obj.fontsize);
            ylabel('(m-m_{min})/M','fontsize',obj.fontsize);
            title(['m_{min}=min(# measurements) at compression c=',num2str(ctest)]);
            hleg=legend('uniform rec','non-unif. rec.','min. UR', 'min. NUR');
            set(hleg,'FontSize',obj.fontsize);
            
        end
        
        
        function show_quality_frameSet(obj, layout, supp_output)
            % ~([layout],[supp_output]) compare quality of hardthresholding
            % recon of several trafos in new window
            % specified in obj.frameSet of type obj.w (e.g. decimated or
            % undecimated);
            % compare PSNR and BSE as a function of the compression factor
            % PSNR ... peak signal to deviation ratio of reconstructed signal
            % BSE ... error of best s-term approximation in p-norm
            % layout={2,1,{0,0,1,2}} to suppress first 2 figures
            % supp_output ... supplementatry output (0 ->none,
            % -1->difference signal, +1 -> reconstructed signal).
            
            if nargin <2 || isempty(layout)
                layout={2,2,{1,2,3,4}};  % {3,2,{1,2,[3,4],[5,6]}}
            end
            if length(layout)<3
                layout{3}={1,2,3,4};
            end
            use_ksdensity=false;
            if nargin <3 || isempty(supp_output)
                supp_output=-1;  % difference signal
            end
            
            p=2;  %l_1 norm
            
            rows=layout{1}; cols=layout{2};
            range=layout{3};
            
            p_str=num2str(p);
            cmax=20; % max compression
            cmin=4;
            goalj=16;
            
            L=length(obj.frameSet);
            Lmax=min(L,21);  % Lmax=21 cylces 3x through all colors
            
            if ~isa(obj.w,'Curvelet2_clab')
                wvl=obj.w;
                cvl=Curvelet2_clab(obj.ts);
            else
                wvl=obj.make_wavelet2D();
                cvl=obj.w;
            end
            
            obj.PSNR=NaN*zeros(cmax,L);
            obj.BSE=obj.PSNR;
            if obj.do_ssim
                obj.SSIM=obj.PSNR;
            else
                obj.SSIM=[];
            end
            s=obj.PSNR;
            
            cvec=1:cmax;  % measurement compression
            
            s1wD=arrayfun(@obj.s1w_Donoho,cvec);
            s1sD=arrayfun(@obj.s1s_Donoho,cvec);
            
            s1wD(cvec<cmin)=NaN;
            s1sD(cvec<cmin)=NaN;
            
            c1wD=obj.ts.numel./s1wD;
            c1sD=obj.ts.numel./s1sD;  % s-compression
            
            
            if supp_output~=0
                obj.algo.variant='reconstruction';
                if supp_output==-1
                    refsignal=obj.img;  % -> difference signal
                else
                    refsignal=zeros(size(obj.img)); % -> reconstructed signal
                end
                hfig1=prepfigure(obj.fig,obj.algo,obj.figopt);
                %sd=factor_subplots(L);
                sd(1)=ceil(sqrt(L));
                sd(2)=ceil(L/sd(1));
                sd=sort(sd,'ascend');
            end
            levels=zeros(1,L);
            wvlwithdwtmodes=cell(1,L);
            s=zeros(length(cmin:cmax),L);
            
            for k=1:L  % loop over trafos
                if obj.verbose && supp_output==0 % default progress indicator
                    fprintf('.');
                end
                fn=obj.frameSet{k};
                w1=wvl.clone;
                if obj.iscurvelet(fn)
                    w1=cvl.clone;
                elseif obj.isfourier(fn)
                    w1=FourierTrafos(obj.ts);
                elseif obj.isinterp(fn)
                    w1=IdTransform(obj.ts);
                    w1.set_dwtmode(obj.frameSet{k});
                    fn='interp';
                end
                wvlwithdwtmodes{k}=[fn,' (',w1.dwtmode_active,')'];
                w1.set_basisname(fn);  % all have same type as w1
                levels(k)=w1.deepestlev;
                
                % $$$ problens with copy and clome
                % signal is changed during loop if obj.pretransformFFT=true
                ts1=obj.ts.clone;                
                if obj.pretransformFFT
                   % apply fft to signal  
                   ts1.xn=ts1.xn-ts1.mean;
                   obj.ts.xn=obj.ts.xn-obj.ts.mean;
%                    f=fft2(obj.ts.xn);                                      
%                    w1.ts.xn=real(f)-imag(f);                                                          
                   w1.ts.xn=dct2(obj.ts.xn); 
                end
                w1.dec;
                for cj=cmin:cmax
                    % uses the same decomposition (condition: obj.C
                    % must not be modified by sparseApprox!
                    if ~obj.isinterp(obj.frameSet{k})
                        [signal2,obj.BSE(cj,k),s(cj,k)]=w1.sparseApprox(c1wD(cj),p);
                        ts1.borderSize=0;
                    else
                        [signal2,obj.BSE(cj,k)]=w1.gridApprox(cj,p);
                        ts1.borderSize=ceil(sqrt(cj)/2);
                    end
                    yn=signal2.xn;
                    if obj.pretransformFFT
                        % apply inverse transform (hartely is its own inverse)
%                         f=fft2(yn);
%                         yn=real(f)-imag(f);    
                        yn=idct2(yn);
                    end
                    obj.PSNR(cj,k)=ts1.PSNR(yn,use_ksdensity);
                    if obj.do_ssim
                        obj.SSIM(cj,k)=ts1.SSIM(yn);
                    end
                    if cj==goalj && supp_output~=0
                        figure(hfig1); % bring to front
                        subplot(sd(1),sd(2),k);
                        imagesc(abs(yn-refsignal)); colormap(gray);
                        title({['error frame ',obj.frameSet{k}, ...
                            ', PSNR=',num2str(obj.PSNR(cj,k),'%3.1f'),' dB'],...
                            ['c=',num2str(goalj),', c_1=', num2str(c1wD(goalj),'%3.1f'),...
                            ', ',w1.title2(levels(k))]},...
                            'fontsize',obj.fontsize);
                        
                        refresh;  % redraw figure (works line progress bar)
                        %  subplot(sd(1),sd(2),k+sd(2));
                        %  [d,err]=ksdensity(abs(yn(:)-w1.ts.xn(:)));
                        %  plot(err,d);
                        %  title('error distribution');
                    end
                end
                
            end % for k (loop over trafos)
            if obj.verbose && supp_output==0  % default progress indicator
                disp('!');
            end
            
            levels=unique(levels);
            
            obj.algo.variant='quality';
            prepfigure(obj.fig, obj.algo,obj.figopt);
            set(0,'DefaultAxesLineStyleOrder','-|--|:');
            
            if range{1}>0
                subplot(rows,cols,range{1});
                obj.ts.graph_signal(false);
            end
            
            if range{2}>0
                subplot(rows,cols,range{2});
                obj.graph_PhaseTransitionDonoho([],false);
                grid on;
                title('sparsity compression c_{1/2}=M/s_{1/2} required (Donoho)',...
                    'fontsize', obj.fontsize);
            end
            
            if range{3}>0
                subplot(rows,cols,range{3});
                semilogx(cvec,obj.PSNR(:,1:Lmax));
                hleg=legend(wvlwithdwtmodes,'location','best');
                set(hleg,'FontSize',obj.fontsize);
                htitle = get(hleg,'Title');
                %  wvl.basisname is enumrated in legend body  not in legend title!
                set(htitle,'String',[wvl.algo.name,'(',wvl.algo.toolbox,')']);
                xlabel('c=M/m','units','normalized','Position',[0.95,-0.08],...
                    'fontsize',obj.fontsize);
                ylabel('PSNR [dB]','fontsize',obj.fontsize);
                tit_add='';
                if use_ksdensity
                    tit_add=[tit_add,'(ksd)'];
                end
                set(gca,'XTick',cmin:obj.cspecial(2));
                grid on;
                line(obj.cspecial(1)*ones(1,2),ylim,'color','black','LineStyle',':');
                line(obj.cspecial(2)*ones(1,2),ylim,'color','black','LineStyle',':');
                xlim([cmin,cmax]);
                title({['PSNR (ideal recon vs. original "', w1.ts.shortname,'")'],...
                    ['at weak c_{1/2}(c) (',w1.title2(levels,false),')']},...
                    'fontsize',obj.fontsize);
            end
            
            if range{4}>0
                subplot(rows,cols,range{4});
                if isempty(obj.SSIM)
                    semilogx(cvec,obj.BSE(:,1:Lmax)./sqrt(s(:,1:Lmax)*obj.ts.numel));
                    ylabel(['\sigma_s(x'')_{',p_str,'}\cdot(sN)^{-1/2}']);
                    title({[' \sigma_s(x'')_{',p_str,'}\cdot(sN)^{-1/2}, at weak c_{1/2}(c)'],...
                        '(error term for sparsity defects)',...
                        },'fontsize',obj.fontsize);
                else
                    semilogx(cvec,obj.SSIM(:,1:Lmax));
                    ylabel('SSIM');
                    title('Strucural similarity (SSIM)','fontsize',obj.fontsize);
                end
                
                
                hleg=legend(wvlwithdwtmodes,'location','best');
                set(hleg,'FontSize',obj.fontsize);
                htitle = get(hleg,'Title');
                %  wvl.basisname is enumrated in legend body not in its title!
                set(htitle,'String',[wvl.algo.name,'(',wvl.algo.toolbox,')']);
                xlabel('c=M/m','units','normalized','Position',[0.95,-0.08],...
                    'fontsize',obj.fontsize);
                set(gca,'XTick',cmin:obj.cspecial(2));
                grid on;
                line(obj.cspecial(1)*ones(1,2),ylim,'color','black','LineStyle',':');
                line(obj.cspecial(2)*ones(1,2),ylim,'color','black','LineStyle',':');
                xlim([cmin,cmax]);
                
            end
            
            obj.algo.variant=[];
            
        end
        
        
        function show_quality_images(obj, layout, supp_output)
            % ~([imagefolder]) compare (in new window) quality of hardthresholding
            % recon of several images where the same sparsifying trafo is
            % used from obj.w.
            % compare PSNR and BSE as a function of the compression factor
            % PSNR ... peak signal to deviation ratio of reconstructed signal
            % BSE ... error of best s-term approximation in p-norm
            % layout={2,1,{0,0,1,2}} to suppress first 2 figures
            % supp_output ... supplementatry output (0 ->none,
            % -1->difference signal, +1 -> reconstructed signal).
            % -2 or +2 shows them in separate figures.
            
            obj.requireOnly(~isempty(obj.imgfiles),'local','there are image file names');
            
            if nargin <2 || isempty(layout)
                layout={2,2,{[1,2],3,4}};  % {1,2,{0,1,2}}
            end
            if length(layout)<3
                layout{3}={[1,2],3,4};
            end
            if nargin <3 || isempty(supp_output)
                supp_output=1;  % show reconstructed image
            end
            
            p=1;  %l_1 norm
            L=length(obj.imgfiles);
            Lmax=min(L,21);  % Lmax=21 cylces 3x through all colors
            use_ksdensity=false;
            
            rows=layout{1}; cols=layout{2};
            range=layout{3};
            
            p_str=num2str(p);
            cmax=20; % max compression
            cmin=4;  % min compression
            goalj=16;
            
            obj.PSNR=NaN*zeros(cmax,L);
            obj.BSE=obj.PSNR;
            s=obj.PSNR;
            if obj.do_ssim
                obj.SSIM=obj.PSNR;
            else
                obj.SSIM=[];
            end
            
            cvec=1:cmax;  % measurement compression
            
            Ltn=L;
            if supp_output~=0
                obj.algo.variant='reconstruction';
                if abs(supp_output)==1
                    hfig1=prepfigure(obj.fig,obj.algo,obj.figopt);
                    % sd=factor_subplots(L);
                    Ltn=min(L,obj.maxThumnails);
                    sd(1)=ceil(sqrt(Ltn));
                    sd(2)=ceil(Ltn/sd(1));
                    sd=sort(sd,'ascend');
                else
                    % output into separate figures
                end
            end
            levels=zeros(1,L);
            shortfn=cell(1,L);
            
            for k=1:L  % loop over images
                if obj.verbose && (supp_output==0 || k>Ltn)  % default progress indicator
                    fprintf('.');
                end
                
                fn=obj.imgfiles{k};
                shortfn{k}=obj.shorten(fn,8);
                try
                    s1=Signal2D.make_fromImage(fn);
                    if ~isa(obj.w,'Wavelet2D_mlab')
                        s1.crop();  % crop to dyadic square
                    end
                catch
                    % skip non-image files
                    continue;
                end
                obj.w.set_signal(s1);
                w1=obj.w;
                
                s1wD=arrayfun(@obj.s1w_Donoho,cvec);
                s1sD=arrayfun(@obj.s1s_Donoho,cvec);
                s1wD(cvec<cmin)=NaN;
                s1sD(cvec<cmin)=NaN;
                c1wD=obj.ts.numel./s1wD;
                c1sD=obj.ts.numel./s1sD;  % s-compression
                
                levels(k)=w1.deepestlev;
                w1.dec;
                
                if supp_output~=0
                    if supp_output<0
                        refsignal=w1.ts.xn;  % -> difference signal
                    else
                        refsignal=zeros(size(w1.ts.xn)); % -> reconstructed signal
                    end
                end
                
                
                for cj=cmin:cmax
                    [signal2,obj.BSE(cj,k),s(cj,k)]=w1.sparseApprox(c1wD(cj),p);
                    yn=signal2.xn;
                    obj.PSNR(cj,k)=w1.ts.PSNR(yn,use_ksdensity);
                    if obj.do_ssim
                        obj.SSIM(cj,k)=w1.ts.SSIM(yn);
                    end
                    if cj==goalj && supp_output~=0 && (isempty(obj.maxThumnails) ||...
                            k<= obj.maxThumnails)
                        if abs(supp_output)==1
                            % show image as thumbnail in matrix subplots
                            figure(hfig1); % bring to front
                            subplot(sd(1),sd(2),k);
                        else % extra image for each, but avoid accumulation
                            % by using prepfigure
                            hfig1=prepfigure(1,struct('mfilename',shortfn{k}),obj.figopt);
                            subplot(1,2,1);
                            imagesc(obj.img); colormap(gray); axis off;
                            title(shortfn{k},'fontsize',obj.fontsize);
                            subplot(1,2,2);
                        end
                        imagesc(abs(yn-refsignal)); colormap(gray);  axis off;
                        adddata=[num2str(obj.PSNR(cj,k),'%3.1f'),' dB'];
                        if obj.do_ssim
                            addata=[adddata,', ',num2str(obj.SSIM(cj,k),'%3.2f')];
                        end
                        
                        title(['R ',shortfn{k},' (',addata,')'],...
                            'fontsize',obj.fontsize);
                        refresh; % redraw figure (works line progress bar)
                        %  subplot(sd(1),sd(2),k+sd(2));
                        %  [d,err]=ksdensity(abs(yn(:)-w1.ts.xn(:)));
                        %  plot(err,d);
                        %  title('error distribution');
                    end
                end % for cj
                
                
            end  % for k (loop over images)
            if supp_output~=0
                suptitle(['Recon: ',w1.basisname,' (',...
                    w1.algo.name,'-',w1.algo.toolbox,')'...
                    ', c=',num2str(goalj),', c_1=', num2str(c1wD(goalj),'%3.1f')],14);
            end
            if obj.verbose && k>Ltn
                disp('!');
            end
            
            levels=unique(levels);
            
            obj.algo.variant='quality';
            prepfigure(obj.fig, obj.algo,obj.figopt);
            set(0,'DefaultAxesLineStyleOrder','-|--|:');
            
            if range{1}>0
                subplot(rows,cols,range{1});
                obj.graph_PhaseTransitionDonoho([],false);
                % plot(c, c1wD,'-r', c, c1sD,':b', obj.cspecial, c1wD(obj.cspecial),'o');
                %
                % textvec([9,16],round([c1wD(9),c1wD(16)]),[],[],true,[]);
                title('sparsity compression c_{1/2}=M/s_{1/2} required (Donoho)',...
                    'fontsize', obj.fontsize);
                % legend('weak c_w(A)','strong c_s(A)',...
                % ['c=',vec2str(obj.cspecial)],'location','best');
                %  xlabel('c=M/m');
                %  ylabel('M/s');
            end
            
            if range{2}>0  % PSNR
                subplot(rows,cols,range{2});
                
                semilogx(cvec,obj.PSNR(:,1:Lmax));
                
                % long legend in general, keep small fontsize
                hleg=legend(shortfn(1:Lmax),'location','southwest');
                htitle = get(hleg,'Title');
                % quote common trafo w1.basisname in legend title
                set(htitle,'String',[w1.algo.toolbox,' (',w1.algo.name,...
                    '): ',w1.basisname]);
                xlabel('c=M/m','units','normalized','Position',[0.95,-0.08],...
                    'fontsize',obj.fontsize);
                ylabel('PSNR [dB]', 'fontsize',obj.fontsize);
                tit_add='';
                if use_ksdensity
                    tit_add=[tit_add,'(ksd)'];
                end
                grid on;
                line(obj.cspecial(1)*ones(1,2),ylim,'color','black','LineStyle',':');
                line(obj.cspecial(2)*ones(1,2),ylim,'color','black','LineStyle',':');
                xlim([cmin,cmax]);
                set(gca,'XTick',cmin:obj.cspecial(2));
                iif  = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
                n2s=@(n) iif(2*floor(n/2)==n,num2str(n),true,'');
                lab=arrayfun(n2s,cmin:obj.cspecial(2),'UniformOutput',false);
                set(gca,'XTickLabel',lab);
                title({'PSNR (ideal recon vs. original)',...
                    ['at weak c_{1/2}(c) (',w1.title2(levels,false),')']},...
                    'fontsize',obj.fontsize);
            end
            
            if range{3}>0 % SSIM
                subplot(rows,cols,range{3});
                if isempty(obj.SSIM)
                    semilogx(cvec,obj.BSE(:,1:Lmax)./sqrt(s(:,1:Lmax)*obj.ts.numel));
                    ylabel(['\sigma_s(x'')_{',p_str,'}\cdot(sN)^{-1/2}']);
                    title({[' \sigma_s(x'')_{',p_str,'}\cdot(sN)^{-1/2}, at weak c_{1/2}(c)'],...
                        '(error term for sparsity defects)',...
                        },'fontsize',obj.fontsize);
                else
                    semilogx(cvec,obj.SSIM(:,1:Lmax));
                    ylabel('SSIM');
                    title('Structural similarity (SSIM)','fontsize',obj.fontsize);
                end
                % long legend in general, keep small fontsize
                hleg=legend(shortfn(1:Lmax),'location','southwest');
                htitle = get(hleg,'Title');
                % quote common trafo w1.basisname in legend title
                set(htitle,'String',[w1.algo.toolbox,' (',w1.algo.name,...
                    '): ',w1.basisname]);
                xlabel('c=M/m','units','normalized','Position',[0.95,-0.08],...
                    'fontsize',obj.fontsize);
                
                grid on;
                line(obj.cspecial(1)*ones(1,2),ylim,'color','black','LineStyle',':');
                line(obj.cspecial(2)*ones(1,2),ylim,'color','black','LineStyle',':');
                xlim([cmin,cmax]);
                set(gca,'XTick',cmin:obj.cspecial(2));
                iif  = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
                n2s=@(n) iif(2*floor(n/2)==n,num2str(n),true,'');
                lab=arrayfun(n2s,cmin:obj.cspecial(2),'UniformOutput',false);
                set(gca,'XTickLabel',lab);
                
            end
            
        end
        
        function show_orig_images(obj,eachinto1figure)
            % ~() show originals of images in new window
            % to compare with recons
            
            if nargin <2
                eachinto1figure=false;
            end
            L=length(obj.imgfiles);
            if ~isempty(obj.maxThumnails)
                L=min(L,obj.maxThumnails);
            end
            
            obj.algo.variant='originals';
            if ~eachinto1figure
                prepfigure(1,obj.algo,obj.figopt);
            end
            %sd=factor_subplots(L);  % not so good if sd is not approx. square
            sd(1)=ceil(sqrt(L));
            sd(2)=ceil(L/sd(1));
            sd=sort(sd,'ascend');
            
            for k=1:L  % loop over images
                fn=obj.imgfiles{k};
                shortfn=obj.shorten(fn,8);
                if ~eachinto1figure
                    subplot(sd(1),sd(2),k);
                else
                    prepfigure(1,struct('mfilename',shortfn),obj.figopt);
                end
                img1=Signal2D.make_fromImage(fn);
                imagesc(img1.xn); colormap(gray); axis off;
                title(shortfn, ...
                    'fontsize',obj.fontsize);
            end
            
            suptitle([num2str(L),' Originals'],14);
            obj.algo.variant=[];
        end
        
        function graph_distribution(obj, probvar, unit, open_new)
            % ~() show distribution density of probvar in open window( default)
            % e.g. probvar=obj.PSNR, SSIM or BSE
            obj.requireOnly(~isempty(probvar),'local','probvar has been computed');
            if nargin <3
                unit=[];
            end
            if nargin <4 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            probvarname=inputname(2);
            L=length(obj.cspecial);
            mu=zeros(1,L);
            sig=mu;
            hold all;
            
            for k=1:L
                mu(k)=meanNaN(probvar(obj.cspecial(k),:));
                sig(k)=stdNaN(probvar(obj.cspecial(k),:));
                ksdensity(probvar(obj.cspecial(k),:));
                % hist(probvar(obj.cspecial(k),:));
            end
            
            grid on;
            
            title({[probvarname,' densities - recon with ',obj.w.basisname,...
                ' (',obj.w.algo.name,')'],...
                ['\mu=',vec2str(mu,'%3.2f'),deblank([' ',unit]),...
                ', \sigma=',vec2str(sig,'%3.2f'),deblank([' ',unit]),...
                ', ',num2str(size(probvar,2)),' images']},...
                'fontsize',obj.fontsize);
            xlabel([probvarname, deblank([' ',unit])],...
                'fontsize',obj.fontsize);
            ylabel('prob. density');
            
            n2s=@(x) num2str(x);
            lab=arrayfun(n2s,obj.cspecial,'UniformOutput',false);
            hleg=legend(lab,'location','best','fontsize',obj.fontsize);
            set(hleg,'FontSize',obj.fontsize);
            htitle = get(hleg,'Title');
            set(htitle,'String','compression c');
            
        end
        
        
        function show_Qdistribution(obj)
            % ~() show distribution desnity of PSNR and, if possible, SSIM
            obj.requireOnly(size(obj.PSNR,1)>max(obj.cspecial),'local',...
                'PSNR matrix includes data for special compression rates');
            
            obj.algo.variant='Quality Recon';
            prepfigure(obj.fig,obj.algo,obj.figopt);
            if ~isempty(obj.SSIM)
                subplot(1,2,1);
            end
            PSNR=obj.PSNR; % necessary to retrieve name in graph_distribution
            obj.graph_distribution(PSNR,'[dB]');
            
            if ~isempty(obj.SSIM)
                subplot(1,2,2);
                SSIM=obj.SSIM;
                obj.graph_distribution(SSIM);
            end
            
            obj.algo.variant=[];
            
        end
        
        function graph_sortedCoeff(obj,func,open_new)
            % ~() energy map of trafo in open window (default)
            % e.g. func.f=@(x) x; func.str= '|.|'; func.descr='modulus';
            %      func.f= @(x) cumsum(x.^2); func.desc='||.||_2'; func.descr='energy map';
            obj.require(isstruct(func) && isa(func.f,'function_handle'),'local',...
                'func.f is function, func.descr string');
            if nargin <3 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            L=length(obj.frameSet);
            Lmax=min(L,21);  % Lmax=21 cylces 3x through all colors
                        
            
            for k=1:Lmax
                
                w1=SparseRep.make_trafo(obj.frameSet{k}, obj.ts);
                w1.dec;
                
                graph_sortedCoeff(w1,func,false);
                hold all;
                
            end
            
            tnames=obj.frameSet(1:Lmax);
            
            title([func.descr,' of size-ordered coeff. ("',obj.signalname,'")'],...
                'fontsize',obj.fontsize);
            
            hleg=legend(tnames,'location','best');
            set(hleg,'FontSize',obj.fontsize);
            htitle = get(hleg,'Title');
            %  wvl.basisname is enumrated in legend body  not in legend title!
            set(htitle,'String',w1.algo.toolbox);
            
        end
        
        
        
        
    end
    
end