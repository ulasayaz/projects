classdef Curvelet2_clab <MultiResTransform2
    % 2-dim. curvelet transformation using the open source toolbox curvelab
    % features:
    %    --- equal norm parseval tight frame
    %    --- all atoms have the same nonunit l2-norm, namely
    %        frame_redundancy^(-1/2)
    %
    %
    % Example: (cf. methods test...()):
    %====================================
    % addpathSplittingSolvers(false); % path needed for curvelets
    % (Source: modified mfiles from chap7, Starck, Sparse Image and signal
    % Processing)
    % --- choose signal
    %     signal1=Signal2D.make_exp_symm();
    %     signal2=Signal2D.make_fromImage('cameraman.bmp');
    %
    % --- discrete curvelet trafo :
    %     w=Curvelet2_clab(); w.show_test(16);
    %     res=w.test_framefeatures(); display(res);
    %
    % --- Zernikes
    %     nr=3; nr=48; nr=150;
    %     sig=0; sig=0.01;
    %     signal1=Signal2D.make_zernike(nr,sig);
    %     w.set_signal(signal1);
    %     w.dec;
    %     w.show_resultd;
    %
    % --- noise filter
    %     sig=0.1;
    %     w= Curvelet2_clab(signal2);
    %     w.graph_test_denoise(sig);
    %
    % --- test image
    %     w= Curvelet2_clab(signal2);
    %     .... limit deepest level of decomposition: w.set_deepestlev(2);
    %     .... compute decomposition:  w.dec;
    %     .... show results: w.show_resultd;
    %     w.figopt.pixsizeX=1000;
    %     w.graph_trafo;
    %     w.show_trafo_components(1);
    %
    % --- compute curvelet frame
    %     w= Curvelet2_clab();
    %     hN=16;  w.computeFrame(hN);
    %     w.set_signal(signal2);
    %     w.test_frame;
    % -- test transform and frame:
    %    res=B.test_framefeatures(); display(res);
    %    [ok,err1,err2]= f.test_DecRec(); disp(ok);
    %    [ok,err1,err2]= f.test_AnalyzeSynthsize(); disp(ok);
    %
    
    properties
        allcurvelets   % @<logical> coeff. at finest level wavelet or curvelet
        % redundancy is smaller if false (ca. 4 instead of 7)
        % best quality for hardthresholding if true.
        is_real        %@<logical> chooses between real and complex valued curvelets
    end
    
    properties (SetAccess=protected)
        lev     %<integer> actual level of curvelet decomp
        
    end
    
    %% constructor and commands
    methods
        function obj=Curvelet2_clab(signal)
            if nargin==0
                signal=Signal2D();
            end
            obj = obj@MultiResTransform2(signal);
            
            obj.repfun=@(x) sign(x).*sqrt(abs(x));
            obj.set_basisname('dctg2');
            obj.frameSet={'dctg2'};
            obj.allcurvelets=true;
            obj.is_real=false;
            obj.frame_issparse=false;  % curvelet basis is not sparse expressed in standard ONB
            % default (only?)border extensionof wavelab is periodic:
            obj.dwtmode_active=[];
            obj.wcL=4;  % $$$ should be 3, but there is an offset by 1 which has to checked
            
        end
        
        
        function set_dwtmode(obj, modestr)
            % ~(modestr) border extension mode cf. help dwtmode
            % e.g. 'sym' (symmetric extension), 'zpd' (zero-padding), etc.
            obj.dwtmode_active=modestr;
        end
        
        function set_algo(obj)
            set_algo@MultiResTransform2(obj);
            obj.algo.name='Curvelet';
            obj.algo.toolbox='curvelab';
        end
        
        
        function dec(obj, lvl)
            % ~() curvelet transform, result obj.C is a cell
            obj.requireOnly(~isempty(obj.img),'local','non-empty sample size');
            obj.lev=obj.deepestlev; % default: floor(log2(min(m,n)))-2)
            if nargin <2 || lvl>obj.lev
                lvl=obj.lev;
            else
                obj.lev=lvl;
            end           
            obj.C=fdct_wrapping_GT(obj.img,obj.is_real, obj.allcurvelets, lvl+1);
            obj.ensureOnly(obj.isTrafodone,'local', 'curvelet trafo done');
        end
        
        function yn=rec(obj,wc)
            %yn=~(cC) % signal reconstruction from curvelet decomp wc;
            % argument wc is of same datatype as result of dec)
            % returns a matrix (image).            
            if isempty(obj.lev)  % rec might be called before dec
                obj.lev=obj.wmaxlev-obj.wcoarsestlev;
            end
            yn=real(ifdct_wrapping_GT(wc,obj.is_real,obj.allcurvelets,obj.lev+1));%, obj.wmaxlev-obj.wcoarsestlev));
        end
        
        
        function yn= analyze(obj,x)
            % yn=~(x) analyze operates like dec but x can be matrix or vector;
            % returns a vector instead of a cell
            obj.lev=obj.deepestlev; % default: floor(log2(min(m,n)))-2)            
            yn=fdct_wrapping_GT(reshape(x,obj.N),obj.is_real, obj.allcurvelets, obj.lev+1);
            yn=obj.cell2vec(yn);
        end
        
        function xn= synthesize(obj,y)
            % xn=~(y) synthesize operates like rec but accepts a vector (datatype
            % of analyze), instead of the original datatype of a decomposition
            % (here cell) and returns a matrix.
            obj.requireOnly(isvector(y),'local','y is vector being the result of analyze');            
            if isempty(obj.lev)  % synthesize might be called before analyze
                obj.lev=obj.wmaxlev-obj.wcoarsestlev;
            end
            xn=real(ifdct_wrapping_GT(obj.vec2cell(y),obj.is_real,obj.allcurvelets,obj.lev+1));
        end
        
        
    end
    
    
    %% queries
    methods
        function ok=isTrafodone(obj)
            % ok=~(): is result of DWT available?
            ok=~isempty(obj.C);
        end
        
        function L=level(obj)
            L= obj.lev;
        end
        
        function nv=frame_norm(obj)
            % l2-norm of 2nd gen. curvelets
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            nv=1/sqrt(obj.frame_redundancy);
        end
        
        function E=frame_norm_exact(obj)
            % E=~() norm of all atoms (exact) returns data compatible with obj.C
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            E = cell(size(obj.C));
            for s=1:length(obj.C)
                E{s} = cell(size(obj.C{s}));
                for w=1:length(obj.C{s})
                    A = obj.C{s}{w};
                    E{s}{w} = sqrt(sum(sum(A.*obj.Conj(A))) / numel(A));
                end
            end
        end
        
        function nv=frame_length(obj)
            % n=~() number of elements of transform
            if ~obj.isTrafodone
                % results from fdct_wrapping_range do not match!
                %  n=max(obj.ts.size);
                % [~,nv] = fdct_wrapping_range(n,nextpow2(n)-obj.wcoarsestlev+1);
                obj.dec;
            end
            
            try
                nv= numel(obj.C{1}{1});
                for s=2:length(obj.C)
                    for w=1:length(obj.C{s})
                        nv=nv+numel(obj.C{s}{w});
                    end
                end
            catch  % empty case
                nv=0;
            end
            
        end
        
        function [img_dwt, levelshown]=C2graph(obj, lvl)
            % img=~() transform curvelet trafo to an image
            % Remark: each corona is individually normalized.
            % 1. The low frequency (coarse scale) coefficients are stored at
            % the center of the display.')
            % 2. The Cartesian concentric coronae show the coefficients at different
            %  scales; the outer coronae correspond to higher frequencies.
            % 3. There are four strips associated to each corona, corresponding to
            %  the four cardinal points; these are further subdivided in angular panels.
            % 4. Each panel represent coefficients at a specified scale and along
            %   the orientation suggested by the position of the panel.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || lvl<=obj.level,'local', 'smaller than max. #levels');
            if nargin < 2 || isempty(lvl)
                lvl=obj.level;
            end
            levelshown=lvl;
            img_dwt= fdct_wrapping_dispcoef(obj.C, lvl+1);
            if ~obj.is_real
                img_dwt=abs(img_dwt);
            end
        end
        
        function vC=C2vec(obj,cC)
            % vec=~() convert transformation result to vector form
            % obj.C contains atoms curvelet frame stored as cell to reduce
            % redundancy
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin <2
                vC = obj.C{1}{1}(:);
                for j=2:length(obj.C)
                    for w=1:length(obj.C{j})
                        vC = [vC;obj.C{j}{w}(:)];
                    end
                end
            else
                vC = cC{1}{1}(:);
                for j=2:length(cC)
                    for w=1:length(cC{j})
                        vC = [vC;cC{j}{w}(:)];
                    end
                end
            end
        end
        
        function nc=numelC(obj,cC)
            % nc=~([cC]) number of elements in decomposition
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            
            if nargin <2
                nc = numel(obj.C{1}{1});
                for j=2:length(obj.C)
                    for w=1:length(obj.C{j})
                        nc = nc+ numel(obj.C{j}{w});
                    end
                end
            else
                nc = numel(cC{1}{1});
                for j=2:length(cC)
                    for w=1:length(cC{j})
                        nc = nc+numel(cC{j}{w});
                    end
                end
            end
        end
        
        function cC=vec2C(obj,vC)
            % cC=~(vC) convert frame coefficients to vector form
            cC=obj.vec2cell(vC);
        end
        
        function cC=vec2cell(obj,vC)
            % converts vector result to cell result
            % using obj.C as pattern
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            % obj.require(numel(vC)==obj.frame_length,'local','same length');
            ptr1=1;
            ptr2=numel(obj.C{1}{1});
            cC{1}{1}=reshape(vC(ptr1:ptr2),size(obj.C{1}{1}));
            cC{1}{2}=obj.C{1}{2};
            for j=2:length(obj.C)
                for w=1:length(obj.C{j})
                    ptr1=ptr2+1;
                    ptr2=ptr2+numel(obj.C{j}{w});
                    cC{j}{w}=reshape(vC(ptr1:ptr2),size(obj.C{j}{w}));
                end
            end
        end
        
        function vC=cell2vec(obj,cC)
            % vec=~() convert transformation result to vector form
            % obj.C contains atoms curvelet frame stored as cell to reduce
            % redundancy
            vC = cC{1}{1}(:);
            for j=2:length(cC)
                for w=1:length(cC{j})
                    vC = [vC;cC{j}{w}(:)];
                end
            end
        end
        
        
        function Ccell=C2cell(obj,lvl)
            % convert transformation result to cell of components
            % must be redefined for those subclasses, where C is not a cell
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || isempty(lvl) || (lvl>=1 && lvl<=obj.level),'local',...
                'level is below computed level');
            if nargin <2
                lvl=[];
            end
            Ccell{1}=obj.appcoef();
            k=1;
            if isempty(lvl)
                for j=1:obj.level
                    for w=1:length(obj.C{j+1})
                        k=k+1;
                        Ccell{k} = obj.detcoef(j,w);
                    end
                end
            else
                j=lvl;
                for w=1:length(obj.C{j+1})
                    k=k+1;
                    Ccell{k} = obj.detcoef(j,w);
                end
            end
        end
        
        function coeff=detcoef(obj,lvl,type)
            % coeff=~(lvl[,type]) extract detail coefficients from curvelet transform
            % type is orientation;
            % if type is missing, it collects all filters of level lev.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(lvl>=1 && lvl<=obj.level,'local',...
                'level is below computed level');
            obj.requireOnly(nargin <3 || isnumeric(type),'local',...
                'detail type is number');
            idx=length(obj.C)+1-lvl;
            if nargin >=3
                coeff=obj.C{idx}{type};
            else
                % recursion
                L=obj.detcoef_rangemax(lvl);
                coeff=cell(1,L);
                for j=1:L
                    coeff{j}=obj.detcoef(lvl,j);
                end
            end
            
        end
        
        function A= appcoef(obj)
            % coeff=~() extract approximate coefficients from curvelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            A=obj.C{1}{1};
        end
        
        
    end
    
    methods (Hidden)
        
        function t=detcoef_rangemax(obj,lvl)
            % t=~() max value of parameter type of the detail coefficients
            t=length(obj.C{end+1-lvl});
        end
        
    end
    
    %% filters
    methods
        
        function estimateNoise(obj)
            % estimate noise by measuring the standard deviation of the first
            % diagonal fluctutation of an orthogonal(!) wavelet transform of NoisyImage.
            % observe, correction factor obj.frame_norm~=1 for curvelets!
            estimateNoise@FrameTrafo(obj);
        end
        
        
        function hardthresh1_inplace(obj,thresh)
            % filter all coefficients obj.C
            % operation on place to save memory
            L=length(obj.C);
            for j=2:L
                for w = 1:length(obj.C{j})
                    obj.C{j}{w} = obj.C{j}{w}.* (abs(obj.C{j}{w})> thresh);
                end
            end
        end
        
        function softthresh1_inplace(obj,thresh)
            % simple soft filter: reduce all coeff. obj.C by size
            % operation on place to save memory
            L=length(obj.C);
            for j=2:L
                for w = 1:length(obj.C{j})
                    obj.C{j}{w} = obj.C{j}{w}.* (abs(obj.C{j}{w})> thresh);
                    obj.C{j}{w}=obj.C{j}{w}-sign(obj.C{j}{w})*thresh;
                end
            end
        end
        
        function hardthresh2_inplace(obj,thresh1, thresh2)
            % filter with 2 thresholds thresh1<thresh2
            % filter only finest details with thresh1
            obj.requireOnly(thresh1<=thresh2,'local', 'thresh1<=thresh2');
            L=length(obj.C);
            for j=2:L
                if j==L
                    thresh=thresh1;
                else
                    thresh=thresh2;
                end
                for w = 1:length(obj.C{j})
                    obj.C{j}{w} = obj.C{j}{w}.* (abs(obj.C{j}{w})> thresh);
                end
            end
        end
        
    end
    
    %% statics
    methods (Static)
        
        function ok= isvalid_framename(wvn)
            ok=ismember(wvn,{'dctg2','dctg2fft'});
        end
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Curvelet2_clab();
        end
        
    end
    %% graphics
    methods
        
        function tit2= title2(obj,levels,~)
            % second title line
            
            if nargin <2 || isempty(levels)
                levels=obj.level;
            end
            tit2=title2@MultiResTransform2(obj,levels,false);
            if obj.allcurvelets
                tit3='finest=curvelet';
            else
                tit3='finest=wavelet';
            end
            tit2=[tit2,', ',tit3,', each corona normalized'];
        end
    end
    
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 2d';
            ok= isa(obj.ts,'Signal2D');
            if ok && ~isempty(obj.ts.xn)
                descr='minimal size';
                ok=min(obj.ts.size(1:2))>16;
            end
        end
    end
    
end




