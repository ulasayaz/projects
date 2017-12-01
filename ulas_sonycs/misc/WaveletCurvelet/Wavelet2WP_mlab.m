classdef Wavelet2WP_mlab < MultiResTransform2
    % 2-dim. wavelet packet transformation using the matlab toolbox <<Wavelet_Toolbox>>   
    %
    % Example: 
    %====================================
    % --- choose signal (image):
    %     signal=Signal2D.make_exp_symm();
    %     signal=Signal2D.make_fromImage('cameraman.bmp');
    %
    % --- discrete wavelet trafo (multi-resolution analysis):
    %     w=Wavelet2WP_mlab();  w.set_basisname('db2'); 
    %     w.show_test;
    %     res=w.test_framefeatures(); display(res);
    % --- set signal: 
    %     w= Wavelet2WP_mlab(signal);
    % --- discrete wavelet packet trafo (multi-resolution analysis):
    %     w.dec;
    %     w.show_resultd;   
    %
   
    
    properties        
        entropy_name  %<string> method to chose optimal subtree
        entropy_param      %<double> method parameter
        chose_besttree %<logical> if true chose best tree with smallest entropy_param       
    end
   
    %% constructor and commands
    methods
        function obj=Wavelet2WP_mlab(signal)
            % constructor
            if nargin==0
                signal=Signal2D();
            end
            obj = obj@MultiResTransform2(signal);
            obj.C=[];
            obj.requireOnly(license('test','Wavelet_Toolbox'),'local','needs wavelet toolbox');
            obj.set_basisname('db1');
            %obj.dwtmode_active=dwtmode('status','nodisp');  % default border extension
            % default border extension: symmetric padding is best for
            % photos:
            obj.dwtmode_active=dwtmode('sym','nodisp');
            obj.entropy_name='norm'; % lp-norm
            obj.entropy_param=1; % p=1: l1-norm 
            obj.chose_besttree=false;   % must be false for compressive sensing
            
        end
        
        function set_dwtmode(obj, modestr)
            % ~(modestr) border extension mode cf. help dwtmode
            % e.g. 'sym' (symmetric extension), 'zpd' (zero-padding),
            % 'ppd' periodic
            obj.dwtmode_active=modestr;
            dwtmode(modestr,'nodisp');
        end
        
        function set_algo(obj)
            set_algo@MultiResTransform2(obj);
            obj.algo.name='WPT'; % wavelet packet trafo
            obj.algo.toolbox='Wavelet\_Toolbox';
            
        end
        
        
        function dec(obj,lvl)
            % ~([lvl]) Wavelet packet decomposition, returns a wptree (!)
            if nargin <2
                lvl=obj.deepestlev();
            end
            obj.require(~isempty(obj.img),'local','non-empty sample size');
            obj.C = wpdec2(obj.img,lvl,obj.basisname,...
                obj.entropy_name, obj.entropy_param);
            if obj.chose_besttree  % else use leaves of complete tree
                obj.C = besttree(obj.C);
            end
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,cC)
            %yn=~(cC) % signal reconstruction from wptree cC;            
            if ~isa(cC,'wptree')
                cC=obj.vec2C(cC);
            end
            yn=wprec2(cC);
        end
        
        
    end
    
    
    %% queries
    methods
        function ok=isTrafodone(obj)
            % ok=~(): is result of DWT available?
            ok=~isempty(obj.C);
        end
        
        function LN=wmaxlev(obj)
            % wavelet decomposition level
            try
                LN=wmaxlev(max(size(obj.img)),obj.basisname);
            catch % unknown wavelets
                LN=max(0,floor(log2(max(size(obj.img)))));
            end
        end
        
        function L=level(obj)
            obj.require(obj.isTrafodone,'local',obj.msgDoDWT);
            L=treedpth(obj.C);
        end
        
        
        function vC=C2vec(obj,cC)
            % vC= ~() convert frame coefficients to vector form
            obj.requireOnly(isa(obj.C,'wptree'),'local','column vector defined for input. Redefine method?');
            %vC= read(obj.C,'cfs',leaves(obj.C));
            if nargin<2
            	vC = read(obj.C,'allcfs');
            else
            	vC = read(cC,'allcfs');
            end
            obj.ensureOnly(isvector(vC),'local','result is column vector');
        end
        
        function cC=vec2C(obj,vC)
            % cC=~(vC) convert frame coefficients to vector form
            obj.require(obj.isTrafodone,'local',obj.msgDoDWT);
            order=4; % for images 4, for 1d-signals 2
            % cfs2wpt needs a row vector !!!
            cC=cfs2wpt(obj.basisname,obj.ts.size,tnodes(obj.C), order,...
                reshape(vC,1,[]));
        end
        
        function coeff=detcoef(obj,level,type)
            % coeff=~(type,level) extract detail coefficient from wavelet transform
            % type ='d' for diagonal fluctutation.
            obj.require(obj.isTrafodone,'local',obj.msgDoDWT);
            obj.requireOnly(level>=1 && level<=obj.level,'local',...
                'level is below computed level');           
            coeff=wpcoef(obj.C,[level,type]);
        end
        
        function A= appcoef(obj, level)
            % coeff=~(type,level) extract detail coefficient from wavelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoDWT);
            A=wpcoef(obj.C,[level,0]);
        end
        
        
        function yn=wpdencmp(obj,thresh)
            % De-noising or compression using wavelet packets
            yn=wpdencmp(obj.ts.xn,'h',obj.deepestlev(),obj.basisname,'threshold',thresh,0);
        end
        
        function [img_dwt, levelshown]=C2graph(obj, lev)
            % img=~(lev) transform transform coefficients to image
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || lev<=obj.level,'local', 'smaller than max. #levels');
            if nargin < 2
                lev=obj.level;
            end
            levelshown=lev;
            % compose wavelet image from parts:
            S=obj.appcoef(lev);
            H=obj.detcoef(lev,1);
            V=obj.detcoef(lev,2);
            D=obj.detcoef(lev,3);
            
            s1=min([size(S,1),size(H,1),size(V,1),size(D,1)]);
            s2=min([size(S,2),size(H,2),size(V,2),size(D,2)]);
            
            img_dwt=[S(1:s1,1:s2), H(1:s1,1:s2);
                     V(1:s1,1:s2),D(1:s1,1:s2)];
            
            
        end
        
        function estimateNoise(obj)
            % estimate noise by measuring the standard deviation of the first
            % diagonal fluctutation of an orthogonal(!) wavelet transform of NoisyImage.
            % observe, correction factor obj.frame_norm~=1 for curvelets!
            estimateNoise@FrameTrafo(obj);
        end
        
    end
    
    methods (Hidden)
        
        function t=detcoef_rangemax(obj,lvl)
            % t=~() max value of parameter type of the detail coefficients
            % is also valid for wave packets (each node can be
            % addressed by 2 indices (level and orientation (0 for trend,
            % 1:3 for detail signal)
            t=3;
        end
        
    end
    
    %% transforms
    methods
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x (like dec), but returns a vector of transform coeff.
            yn = wpdec2(reshape(x,obj.N),obj.deepestlev(),obj.basisname,...
                obj.entropy_name, obj.entropy_param);
            if obj.chose_besttree % else use leaves of complete tree
                yn = besttree(yn);
            end            
            yn=read(yn,'allcfs');
        end
        
        function xn= synthesize(obj,y)
            % xn=~(y) synthesize operates like rec but needs a vector input
            obj.requireOnly(isvector(y),'local','y is vector being the result of analyze');            
            xn=wprec2(obj.vec2C(y));
        end
        
        
    end
    
    %% filters
    methods
        
        function [signal2, thresh, fCoeff]=hardthresh1(obj,thresh)
            % [signal2,thresh.fCoeff]=~([thresh]) hard threshold transform coefficients  
            % redefined from class FrameTrafo.   
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin < 2 || isempty(thresh) || isnan(thresh)
                obj.estimateNoise();
                thresh=4*obj.frame_norm *obj.sigma;
            end
            
            oldC=obj.C; % must be clone
            
            cC=obj.C2vec; 
            cC(abs(cC)<thresh)=0;
            obj.C=obj.vec2C(cC);
            
            % reconstructed signal
            yn = obj.rec(obj.C) ;
            if nargin >2
                fCoeff=obj.C;
            end
            obj.C=oldC;
            
            signal2=obj.ts.clone();
            signal2.replace_signal(yn);
            signal2.signalname=[obj.ts.signalname,', hard thresh. (',obj.basisname,...
                ', ',num2str(thresh,'%3.1e'),')'];
            signal2.colormap_active=obj.ts.colormap_active;
            
        end
        
        function [signal2, thresh1, thresh2,fCoeff]=hardthresh2(obj,thresh1, thresh2)
            % [signal2,thresh1,thresh2]=~(thresh) denoise using coefficient thresholding
            % cf. hardthresh2 using 2 thresholds for different kinds of
            % coefficients.
            % $$$ tO BE IMPLEMENTED: cases for thresh1 missing
            
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
            
            cC=obj.C2vec;
            cC(abs(cC)<thresh2)=0;
            obj.C=obj.vec2C(cC);
            
            % reconstructed signal
            yn = obj.rec(obj.C) ;
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
            % [signal2,thresh]=~(thresh) denoise using coefficient thresholding
            % cf. hardthresh1,2
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin < 2 || isempty(thresh)
                obj.estimateNoise();
                thresh=4*obj.frame_norm *obj.sigma;
            end
            
            oldC=obj.C; % must be clone
            
            fCoeff=obj.C2vec;
            fCoeff(abs(fCoeff)<thresh)=0;
            % now also reduce other coefficients
            fCoeff=fCoeff-sign(fCoeff)*thresh;
            obj.C=obj.vec2C(cC);
                        
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
        
        
    end
    
    %% tests
    methods (Static)
        
        function str=msgDoDWT()
            str='Wavelet packet decomposition available (call dec).';
        end
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Wavelet2WP_mlab();
        end
        
    end
    
    %% graphics
    methods
        
    end
    
end

