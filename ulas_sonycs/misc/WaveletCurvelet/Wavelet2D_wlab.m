classdef Wavelet2D_wlab <MultiResTransform2
    % 2-dim. decimated wavelet transformations using the open source toolbox wavelab
    % set pathes by calling WavePath;
    %
    % features:
    %    --- both discrete multi-resolution analysis (.dec)
    %    --- noise filter (.show_noisefilter)
    %
    %
    % Example:
    % ====================================
    % --- choose signal (image):
    %     signal1=Signal2D.make_exp_symm();
    %     signal2=Signal2D.make_fromImage('cameraman.bmp');
    % --- discrete wavelet trafo (multi-resolution analysis):
    %     w=Wavelet2D_wlab(); w.show_test;
    % --- db2 (not db1) sensitive to discont. in 1st derivative of signal exp_symm
    %     w=Wavelet2D_wlab(signal1);
    %     .... db2,3 show edges at border and along image axis:
    %     w.show_trafos_selection({'db1','db2','db3'});
    % --- noise filter
    %     omega0=5;sig=0.2;hN=512;
    %     w= Wavelet2D_wlab(Signal2D.make_sinus(omega0,sig,hN));
    %     w.dec;
    %     w.estimateNoise; disp(w.sigma);
    %     w.show_resultd;
    %     w.show_noisefilter;
    %
    % --  photos have many edges where 1st derivative is discontinuous
    %     w= Wavelet2D_wlab(signal2);
    %     .... limit deepest level of decomposition: w.set_deepestlev(2);
    %     .... compute decomposition:  w.dec;
    %     .... show results: w.show_resultd;
    %     w.figopt.pixsizeX=1000;
    %     w.graph_trafo;
    %     w.show_trafo_components(1);
    %
    
    properties
        
        % computed values:
        qmf    %@<vector> Orthonormal QMF Filter for Wavelet Transform
        dqmf   %@<vector> decomposition filter for biorthogonal wavelets
    end
    
    properties (SetAccess=protected)
        lev     %@ (int) actual level of wavelet decomp
        OW_mlab2wlab, OW_wlab2mlab   % name maps for orthogonal wavelet families
        BOW_mlab2wlab, BOW_wlab2mlab % % name maps for biorthogonal wavelet families
    end
    
    %% constructor and commands
    methods
        function obj=Wavelet2D_wlab(signal)
            global WAVELABPATH
            if nargin==0
                signal=Signal2D([]);
            end
            obj = obj@MultiResTransform2(); % cannot set signal before filter)
            obj.requireOnly(~isempty( WAVELABPATH),'local', ...
                'needs  WAVELABPATH (call WavePath.m).');
            % cf. MakeONFilter.m:
            
            orth_wlab={'Daubechies','Haar','Coiflet','Symmlet'};
            orth_mlab={'db','haar','coif','sym'};
            obj.OW_mlab2wlab=containers.Map(orth_mlab,orth_wlab);
            obj.OW_wlab2mlab=containers.Map(orth_wlab,orth_mlab);
            
            biorth_wlab={'CDF'};
            biorth_mlab={'bior'};
            obj.BOW_mlab2wlab=containers.Map(biorth_mlab,biorth_wlab);
            obj.BOW_wlab2mlab=containers.Map(biorth_wlab,biorth_mlab);
            
            obj.set_basisname('db1');
            obj.set_signal(signal);  % set again to padd if required
            
            % default (only?)border extensionof wavelab is periodic:
            obj.dwtmode_active='ppd';
            obj.wcL=2;
            
            % detail signal:
            detail_labels={'v','h','d'};
            detail_idx=num2cell(1:3);
            obj.detail_str2idx=containers.Map(detail_labels,detail_idx);
            obj.detail_idx2str=containers.Map(detail_idx,detail_labels);
            
        end
        
        function set_signal(obj,hts)
            % ~(hts) set time signal to hts of class SignalClass
            if ~hts.isdyadicHypercube
                hts.padd2dyadic(0,'post');
                warning('signal has been padded to a dyadic length');
            end
            set_signal@FrameTrafo(obj,hts);
        end
        
        function set_basisname(obj, wvn)
            % ~(wvn) set wavelet name and reset trafo
            obj.requireOnly(ischar(wvn),'local',' wvn is wavelet name');
            set_basisname@MultiResTransform2(obj, wvn);
            nr=obj.wvnumber;
            if obj.OW_mlab2wlab.isKey(obj.wvfamily)
                ow=obj.OW_mlab2wlab(obj.wvfamily);
                obj.qmf=MakeONFilter(ow,2*nr);  % associate filter with wavelet
            elseif obj.BOW_mlab2wlab.isKey(obj.wvfamily)
                biow=obj.BOW_mlab2wlab(obj.wvfamily);
                par=num2str(nr);  % assume par=='x.y'
                par1=str2double(par(1)); par2=str2double(par(3));
                [obj.qmf, obj.dqmf]=MakeBSFilter(biow,[par1,par2]);  % associate filter with wavelet
            else
                error('The specified key is not present in  the container of wavelet names');
            end
            obj.C=[];
            obj.ensureOnly(~obj.isTrafodone,'local', 'trafo reset');
        end
        
        function set_dwtmode(obj, modestr)
            % ~(modestr) border extension mode cf. help dwtmode
            % e.g. 'sym' (symmetric extension), 'zpd' (zero-padding), etc.
            % ??? obj.dwtmode_active=modestr;
        end
        
        function set_algo(obj)
            set_algo@MultiResTransform2(obj);
            obj.algo.name='DWT';   % decimated wavelet trafo
            obj.algo.toolbox='wavelab';
        end
        
        
        function dec(obj,lvl)
            % discrete multiresolution analysis
            % creates obj.C as a matrix of same shape as signal
            obj.requireOnly(~isempty(obj.img),'local','non-empty sample size');
            obj.requireOnly(~isempty(obj.qmf),'local','needs filter');
            if nargin <2
                lvl=obj.deepestlev;
            end
            obj.lev=lvl;
            if obj.OW_mlab2wlab.isKey(obj.wvfamily) % orhtogonal WT
                obj.C = FWT2_PO(obj.img,obj.level_bottomUp,obj.qmf);
            elseif  obj.BOW_mlab2wlab.isKey(obj.wvfamily)  % bi-orthogonal WT
                % symmetric border extension improves PSNR by ca. 0.5dB
                obj.C = FWT2_SBS(obj.img,obj.level_bottomUp,obj.qmf, obj.dqmf);
                %obj.C = FWT2_PB(obj.img,obj.level_bottomUp,obj.qmf, obj.dqmf);
            end
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wc)
            %yn=~(cC) % signal reconstruction from wavelet decomp wc
            obj.requireOnly(~isempty(obj.qmf),'local','needs filter');
            if isvector(wc)
                wc=reshape(wc,obj.ts.N);
            end
            if isempty(obj.lev) % rec might be called before dec
                obj.lev=obj.deepestlev;
            end
            if obj.OW_mlab2wlab.isKey(obj.wvfamily) % orhtogonal WT
                yn= IWT2_PO(wc,obj.level_bottomUp,obj.qmf);
            elseif obj.BOW_mlab2wlab.isKey(obj.wvfamily)  % bi-orthogonal WT
                % symmetric border extension improves PSNR by ca. 0.5dB
                yn= IWT2_SBS(wc,obj.level_bottomUp,obj.qmf, obj.dqmf);
                %yn= IWT2_PB(wc,obj.level_bottomUp,obj.qmf, obj.dqmf);
            end
            
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x
            obj.lev=obj.deepestlev;
            yn=FWT2_PO(reshape(x,obj.N),obj.level_bottomUp,obj.qmf);
        end
        
        function xn= synthesize(obj,y)
            % xn=~(x) synthesize
            if isempty(obj.lev) % synthesize might be called before analyze
                obj.lev=obj.deepestlev;
            end
            xn= IWT2_PO(reshape(y,obj.N),obj.level_bottomUp,obj.qmf);
        end
        
        
    end
    
    
    %% queries
    methods
        function ok=isTrafodone(obj)
            % ok=~(): is result of DWT available?
            ok=~isempty(obj.C);
        end
        
        function L=wmaxlev(obj)
            % maximal useful wavelet decomposition level
            % analog wavelet toolbox of matlab
            obj.requireOnly(~isempty(obj.qmf),'local','needs filter');
            
            lw = length(obj.qmf);
            lx=obj.N;
            lx=lx(1);
            L = fix(log(lx/(lw-1))/log(2));
            if L<1 , L = 0; end
            
        end
        
        function L=level(obj)
            L= obj.lev;
        end
        
        function n=wvsupport(obj)
            % length of support of wavelet
            n=length(obj.qmf);
        end
        
        function [img_dwt, levelshown]=C2graph(obj, lev)
            % img=~() transform dwt to image
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || lev<=obj.level,'local', 'smaller than max. #levels');
            if nargin < 2
                lev=obj.lev;
            end
            levelshown=obj.level;
            img_dwt=reshape(obj.C,obj.ts.N);
        end
        
        
        function coeff=detcoef(obj,level,type)
            % coeff=~(level[,type]) extract detail coefficients from wavelet transform
            % type ='d' for diagonal fluctutation;
            % if type is missing, it collects all filters of level lev.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(level>=1 && level<=obj.level,'local',...
                'level is below computed level');
            obj.requireOnly(nargin <3 || ischar(type) || type<=3,'local',...
                'type is number or "h" or "v" or "d".');
            d=round(log2(obj.N));
            ls=2.^(d-level+1);
            if nargin >=3
                if isnumeric(type)
                    type=obj.detail_idx2str(type);
                end
                switch type
                    case 'd'
                        coeff=obj.C(floor(ls(1)/2)+1:ls(1), floor(ls(2)/2)+1:ls(2));
                    case 'v'
                        coeff=obj.C(1:floor(ls(1)/2), floor(ls(2)/2)+1:ls(2));
                    case 'h'
                        coeff=obj.C(floor(ls(1)/2)+1:ls(1), 1:floor(ls(2)/2));
                    otherwise
                        error('undefined type');
                end
            else
                % recursion
                L=obj.detcoef_rangemax(level);
                coeff=cell(1,L);
                for j=1:L
                    coeff{j}=obj.detcoef(level,j);
                end
            end
        end
        
        function A= appcoef(obj)
            % coeff=~() extract approximate coefficients from wavelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            d=round(log2(obj.N));
            A=obj.C(1:2^(d(1)-obj.level),1:2^(d(2)-obj.level));
        end
        
        function cC=vec2C(obj,vC)
            % convert vector to form usable for reconstruction rec
            cC=reshape(vC,obj.ts.N);
        end
        
        
    end
    
    %% filters
    methods
        
        
        function hardthresh2_inplace(obj,thresh1, thresh2)
            % filter with 2 thresholds thresh1<thresh2
            obj.requireOnly(thresh1<=thresh2,'local', 'thresh1<=thresh2');
            coarsest=obj.wcoarsestlev;
            ll  = obj.C(1:coarsest,1:coarsest);
            obj.C(N/2+1:N,1:N/2)= (abs(obj.C(N/2+1:N,1:N/2))> thresh1).*obj.C(N/2+1:N,1:N/2);   % Finest vertical.
            obj.C(1:N/2,N/2+1:N)= (abs(obj.C(1:N/2,N/2+1:N))> thresh1).*obj.C(1:N/2,N/2+1:N);   % Finest horizontal.
            obj.C(N/2+1:N,N/2+1:N)=(abs(obj.C(N/2+1:N,N/2+1:N))> thresh1).*obj.C(N/2+1:N,N/2+1:N); % Finest diagonal.
            obj.C(1:N/2,1:N/2)= (abs(obj.C(1:N/2,1:N/2)) > thresh2).*obj.C(1:N/2,1:N/2);     % Other scales.
            obj.C(1:coarsest,1:coarsest) = ll;
            obj.C=wc;
        end
        
        function [nodes,mask]=nodesPDF(obj,c,params)
            % nodes=~(c) nodes sampled with compression rate c
            % trend part of signal is sampled denser than detail part;
            % for the matlab toolbox the trend signal is in the beginning
            % of the decomposition vector.
            obj.requireOnly(c>1,'local','compression rate in ]0,1[');
            obj.requireOnly(~obj.ts.isemptydata,'local','non-empty signal');
            if nargin<3 || ~isfield(params,'P')
                params.P=5; % Variable density polymonial degree
            end
            if ~isfield(params,'iter')
                params.iter=1;
            end
            obj.dec(1);
            trend=obj.appcoef();
            L=numel(trend);
            p=1/c;
            N=obj.numelC();
            r=0.85/0.15;
            % p*N==p1*L+p2*(N_L)
            % p1/p2=r;
            p2=p*N/(r*L+N-L);
            p1=r*p2;
            if p1>1  % we can choose all trend coefficients!
                p1=1;
                p2=(p*N-L)/(N-L);
            end
            
            % generate sampling mask: p2 for detail part (here for all in a
            % first step):
            pdf=p2*ones(obj.ts.N);
            lvl=1;
            d=round(log2(obj.ts.N));
            Ntrend=[2^(d(1)-lvl),2^(d(2)-lvl)];
            % trend signal gets p1:
            pdf(1:Ntrend(1),1:Ntrend(2))=p1*ones(Ntrend);
            
            assert(abs(sum(pdf(:))-p*N)<1,'definition of pdf');
            mask = genSampling(pdf,params);
            % compute nodes:
            nodes=SampledNodes(find(mask),obj.ts.N,obj);
        end
        
        function nodes=nodesDefault(obj,c,params)
            % nodes=~() default nodes
            if nargin <2 || isempty(c)
                c=4;
            end
            if ~exist('params','var')
                params=struct;
            end
            nodes=obj.nodesPDF(c,params);
        end
        
    end
    
    %% tests
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Wavelet2D_wlab();
        end
        
        function w=show_test(w)
            % w=~(ts) testing the wavelet trafo with signal ts
            assert(nargin <1 || isa(w,'Wavelet2D_wlab'),'w belongs to wavelet class');
            if nargin ==0
                hN=256; v=8;
                hts=Signal2D.make_star(v,hN);
                w=Wavelet2D_wlab(hts);
                w.figopt.pixsizeX=1000;
            end
            
            w.show_trafos_selection({'db1','db2','db3'});
            
        end
        
    end
    %% graphics
    methods
        function cC=dec2graph(obj,cC)
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.C=cC;
            cC=obj.C2graph;
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 2d';
            ok= isa(obj.ts,'Signal2D');
            if ok
                descr='signal is a dyadic square';
                ok=obj.ts.isdyadicHypercube();
            end
        end
    end
    
end



