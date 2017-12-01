classdef Wavelet2D_mlab <MultiResTransform2
    % 2-dim. decimated wavelet transformation using the matlab toolbox <<Wavelet_Toolbox>>
    % features:
    %    --- both discrete multi-resolution analysis (.wavedec)
    %    --- noise filter (.show_noisefilter)
    %
    %
    % Example:
    % ====================================
    % --- choose signal (image):
    %     signal1=Signal2D.make_exp_symm();
    %     signal2=Signal2D.make_fromImage('cameraman.bmp');
    %
    % --- discrete wavelet trafo (multi-resolution analysis):
    %     w=Wavelet2D_mlab();  w.set_basisname('db2'); w.show_test;
    %     res=w.test_framefeatures(); display(res);
    %
    % --- db2 (not db1) sensitive to discont. in 1st derivative of signal exp_symm
    %     w.set_signal(signal1);
    %     w.set_dwtmode('sym');
    %     -- db2,3 show edges at border and along image axis:
    %     w.show_trafos_selection({'db1','db2','db3'});
    %     w.fig=2;
    %     w.set_dwtmode('sp1');
    %     -- only edge along image axis remains:
    %     w.show_trafos_selection({'db1','db2','db3'});
    %
    % --- noise filter
    %     sig=0.1;
    %     w= Wavelet2D_mlab(signal2);
    %     w.graph_test_denoise(sig);
    %
    % --  photos have many edges where 1st derivative is discontinuous
    %     w= Wavelet2D_mlab(signal2);
    %     .... limit deepest level of decomposition: w.set_deepestlev(2); 
    %     .... compute decomposition:  w.dec;
    %     .... show results: w.show_resultd;
    %     w.figopt.pixsizeX=1000;
    %     w.graph_trafo;
    %     w.show_trafo_components(1);
    %
    % -- compute wavelet basis
    %     w= Wavelet2D_mlab();
    %     hN=16;  w.set_basisname('db1'); w.computeFrame(hN);
    %     w.set_signal(signal2);
    %     w.test_frame;
    %
    
    properties
        
        % computed values:
        S    %@<vector> wavelet bookkeeping vector S
    end
    
    %% constructor and commands
    methods
        function obj=Wavelet2D_mlab(signal)
            % constructor
            if nargin==0
                signal=Signal2D();
            end
            obj = obj@MultiResTransform2(signal);
            obj.requireOnly(license('test','Wavelet_Toolbox'),'local','needs wavelet toolbox');
            obj.set_basisname('db1');
            %obj.dwtmode_active=dwtmode('status','nodisp');  % default border extension
            % default border extension: sp1 (smooth padding) works well for
            % smooth signals, 'sym' works well for photos.
            % photos:
            obj.dwtmode_active=dwtmode('sym','nodisp');
            
            % detail signal:
            detail_labels={'v','h','d'};
            detail_idx=num2cell(1:3);
            obj.detail_str2idx=containers.Map(detail_labels,detail_idx);
            obj.detail_idx2str=containers.Map(detail_idx,detail_labels);
            
        end
        
        function set_dwtmode(obj, modestr)
            % ~(modestr) border extension mode cf. help dwtmode
            % e.g. 'sym' (symmetric extension), 'zpd' (zero-padding),
            % 'ppd' periodic, 'sp1 ' smooth padding
            obj.dwtmode_active=modestr;
            dwtmode(modestr,'nodisp');
        end
        
        function set_algo(obj)
            set_algo@MultiResTransform2(obj);
            obj.algo.name='DWT'; % decimated wavelet trafo
            obj.algo.toolbox='Wavelet\_Toolbox';
            
        end
        
        
        function dec(obj,lev)
            % discrete multiresolution analysis
            % creates obj.C as a vector.
            obj.require(~isempty(obj.img),'local','non-empty sample size');
            if nargin <2
                lev=obj.deepestlev();
            end
            % wavelet decomposition from matlab's wavelet toolbox            
            [obj.C,obj.S] = wavedec2(obj.img,lev,obj.basisname);
            
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wc)
            %yn=~(cC) % signal reconstruction from wavelet decomp wc
            % wavelet recon from matlab's wavelet toolbox
            obj.requireOnly(~isempty(obj.S),'local','needs recon info');
            yn = waverec2(wc,obj.S,obj.basisname) ;
        end
        
        
    end
    
    
    %% queries
    methods
        
        function ok=isTrafodone(obj)
            % ok=~(): is result of DWT available?
            ok=~isempty(obj.C) && ~isempty(obj.S);
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x and returns a vector of transform coeff.
            yn=wavedec2(reshape(x,obj.N),obj.deepestlev(),obj.basisname);
        end
        
        function xn= synthesize(obj,y)
            % xn=~(y) synthesize, i.e. returns signal content as matrix
            obj.requireOnly(~isempty(obj.S),'local','needs recon info, call dec.');
            xn=waverec2(y,obj.S,obj.basisname) ;
        end
        
        function LN=wmaxlev(obj)
            % wavelet decomposition level
            if ~isempty(obj.basisname)
                try
                    LN=wmaxlev(max(size(obj.img)),obj.basisname);
                catch
                    % unknown wavelets                    
                    LN=max(0,floor(log2(max(obj.size))));
                end
            else
                LN=max(0,min(floor(log2(min(size(obj.img))))));
            end
        end
        
        function L=level(obj)
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            L= size(obj.S,1)-2;
        end
        
        
        function [phi,psi,xval] = wavefun(obj,iter)
            % [phi,psi,xval]=~(iter) returns scaling and wavelet function
            [phi,psi,xval] = wavefun(obj.basisname,iter);
        end
        
        function coeff=detcoef(obj,level,type)
            % coeff=~(level[,type]) extract detail coefficients from wavelet transform
            % e.g. type ='d' for diagonal fluctutation;
            % if type is missing, it collects all filters of level lev.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(level>=1 && level<=obj.level,'local',...
                'level is below computed level');
            obj.requireOnly(nargin <3 || ischar(type) || type<=3,'local',...
                'type is number or "h" or "v" or "d".');            
            if nargin>=3
                if isnumeric(type)
                    type=obj.detail_idx2str(type);
                end
                coeff=detcoef2(type,obj.C,obj.S,level);
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
            % coeff=~() extract deepest level approximate coefficients
            % from wavelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            A= appcoef2(obj.C,obj.S,obj.basisname,obj.level);
        end
        
        function [img_dwt, levelshown]=C2graph(obj, lev)
            % img=~(lev) transform transform coefficients to image
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || lev<=obj.level,'local', 'smaller than max. #levels');
            if nargin < 2
                lev=obj.level;
            end
            % compose wavelet image from parts:            
            img_dwt=zeros(sum(obj.S(1:end-1,1)), sum(obj.S(1:end-1,2)));
            levelshown=size(obj.S,1)-2;
            rows=cumsum(obj.S(1:end-1,1));
            cols=cumsum(obj.S(1:end-1,2));
            img_dwt(1:rows(1),1:cols(1))=...
                appcoef2(obj.C,obj.S,obj.basisname,levelshown);
            for j=levelshown:-1:1
                k=levelshown-j+1;
                img_dwt(1:obj.S(k+1,1),cols(k)+1:cols(k+1))=...
                    detcoef2('h',obj.C,obj.S,j);
                img_dwt(rows(k)+1:rows(k+1),1:obj.S(k+1,2))=...
                    detcoef2('v',obj.C,obj.S,j);
                img_dwt(rows(k)+1:rows(k+1),cols(k)+1:cols(k+1))=...
                    detcoef2('d',obj.C,obj.S,j);
            end
            
        end
        
        function cC=dec2graph(obj,cC)
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.C=cC;
            cC=obj.C2graph;
        end
        
    end
    
    
    %% filters
    methods
        
        
        function hardthresh2_inplace(obj,thresh1, thresh2)
            % filter with 2 thresholds thresh1<thresh2
            obj.requireOnly(thresh1<=thresh2,'local', 'thresh1<=thresh2');
            % filter trend signal with thresh2
            obj.C(abs(obj.C)<thresh2)=0;
            
            % filter detail signal with thresh1
            cD=obj.C(obj.S(1,1)*obj.S(1,2)+1:end);
            cD(abs(cD)<thresh1)=0;
            % combine again with trend signal:
            obj.C(obj.S(1,1)*obj.S(1,2)+1:end)=cD;
        end
        
        function [nodes,mask]=nodesPDF(obj,c,params)
            % nodes=~(c) nodes sampled with compression rate c
            % trend part of signal is sampled denser than detail part;
            % for the matlab toolbox the trend signal is in the beginning
            % of the decomposition vector.
            obj.requireOnly(isscalar(c) && c>1,'local','compression rate in ]0,1[');
            obj.requireOnly(~obj.ts.isemptydata,'local','non-empty signal');            
            if nargin<3 || ~isfield(params,'P')
                params.P=5; % Variable density polymonial degree
            end
            if ~isfield(params,'iter')
                params.iter=1;
            end
            if ~isfield(params,'w1')
                params.w1=0.85;  % weight for trend signal
            end
            
            obj.dec(1);
            trend=obj.appcoef();
            L=numel(trend);
            p=1/c;
            N=obj.numelC();  
            r=params.w1/(1-params.w1);  % ratio of weights of trend and detail signals
            % p*N==p1*L+p2*(N_L)
            % p1/p2=r;
            p2=p*N/(r*L+N-L);
            p1=r*p2;
            if p1>1  % we can choose all trend coefficients!
                p1=1;
                p2=(p*N-L)/(N-L);
            end
           
            % generate sampling mask:
            pdf= [p1*ones(L,1); p2*ones(N-L,1)];   % trend signal is in the head  
            assert(abs(sum(pdf)-p*N)<1,'definition of pdf');
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
    
    %% statics
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Wavelet2D_mlab();
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
        
    end
    
end

