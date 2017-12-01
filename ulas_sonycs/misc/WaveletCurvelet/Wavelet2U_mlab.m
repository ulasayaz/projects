classdef Wavelet2U_mlab <MultiResTransform2
    % 2-dim. undecimated decimated wavelet transformation using the matlab toolbox
    % <<Wavelet_Toolbox>>
    % features:
    %    --- both discrete multi-resolution analysis (.wavedec)
    %    --- noise filter (.show_noisefilter)
    %
    %
    % Example:
    %====================================
    % --- choose signal (image):
    %     signal1=Signal2D.make_exp_symm();
    %     signal2=Signal2D.make_fromImage('cameraman.bmp'):
    % --- discrete wavelet trafo (multi-resolution analysis):
    %     w=Wavelet2U_mlab(); w.set_basisname('db2'); w.show_test;
    % --- db2 (not db1) sensitive to discont. in 1st derivative of signal exp_symm
    %     w=Wavelet2U_mlab(signal1);
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
    %     w= Wavelet2U_mlab(signal2);
    %     w.graph_test_denoise(sig);
    %
    % --  photos have many edges where 1st derivative is discontinuous
    %     w= Wavelet2U_mlab(signal2);
    %     w.set_deepestlev(2); w.dec;
    %     w.show_resultd;
    %     w.figopt.pixsizeX=1000;
    %     w.graph_trafo;
    %     w.show_trafo_components(1);
    %
    
    properties
        
    end
    
    %% constructor and commands
    methods
        function obj=Wavelet2U_mlab(signal)
            if nargin==0
                signal=[];
            end
            obj = obj@MultiResTransform2(signal);
            obj.requireOnly(license('test','Wavelet_Toolbox'),'local','needs wavelet toolbox');
            obj.set_basisname('db1');
            %obj.dwtmode_active=dwtmode('status','nodisp');  % default border extension
            % default border extension: symmetric padding is best for
            % photos:
            obj.dwtmode_active=dwtmode('sym','nodisp');
            
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
            obj.algo.name='UWT'; % undecimated wavelet trafo
            obj.algo.toolbox='Wavelet\_Toolbox';
            
        end
        
        
        function dec(obj, lvl)
            % ~([lvl])discrete multiresolution analysis
            obj.require(~isempty(obj.img),'local','non-empty sample size');
            if nargin <2
                lvl=obj.deepestlev();
            end
            % wavelet decomposition from matlab's wavelet toolbox
            % obj.C.dec: order: 'a', 'v','h','d'
            obj.C =ndwt2(obj.img,lvl,obj.basisname);
            
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wc)
            %yn=~(cC) % signal reconstruction from wavelet decomp wc
            % wc must be of type obj.C
            % wavelet recon from matlab's wavelet toolbox
            yn =indwt2(wc);
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x 
            yn=ndwt2(x,obj.deepestlev(),obj.basisname);
        end
        
        function xn= synthesize(obj,y)
            % xn=~(x) synthesize
            xn =indwt2(y);
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
            if ~isempty(obj.basisname)
                try
                    LN=wmaxlev(min(size(obj.img)),obj.basisname);
                catch
                    % unknown wavelets
                    LN=max(0,min(floor(log2(min(size(obj.img))))));
                end
            else
                LN=max(0,min(floor(log2(min(size(obj.img))))));
            end
        end
        
        function L=level(obj)
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            L= obj.C.level;
        end
        
        
        function nv=frame_length(obj)
            % n=~() number of elements of transform
            nv=0;
            try
                for j=1:length(obj.C.dec )
                    nv= nv+numel(obj.C.dec{j});
                end
            catch  % empty case
                nv=0;
            end
        end
        
        function coeff=detcoef(obj,lvl,type)
            % coeff=~(level[,type]) extract detail coefficients from wavelet transform
            % e.g. type ='d' for diagonal fluctutation.
            % if type is missing, it collects all filters of level lev.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(lvl>=1 && lvl<=obj.level,'local',...
                'level is below computed level');
            obj.requireOnly(nargin <3 || ischar(type) || type<=3,'local',...
                'type is number or "h" or "v" or "d".');
            dlvl=obj.level;
            if nargin >=3
                if isnumeric(type)
                    type=obj.detail_idx2str(type);
                end
                switch type
                    case 'd'
                        idx=3*(dlvl-lvl)+1+3;
                    case 'v'
                        idx=3*(dlvl-lvl)+1+1;
                    case 'h'
                        idx=3*(dlvl-lvl)+1+2;
                    otherwise
                        error('undefined type');
                end
                coeff=obj.C.dec{idx};
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
            % coeff=~() extract approximate coefficients from wavelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            A = obj.C.dec{1,1};
        end
        
        function dw=copy2decimated(obj)
            % dw=~() make copy as decimated wavelet trafo
            dw=Wavelet2D_mlab(obj.ts);
            dw.set_basisname(obj.basisname);
            dw.dwtmode_active=obj.dwtmode_active;
            dw.wcL=obj.wcL;
        end
        
        function [img_dwt, levelshown]=C2graph(obj, lev)
            % transform dwt to image (for UWT show 4 parts)
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || lev<=obj.level,'local', 'smaller than max. #levels');
            if nargin < 2
                lev=obj.level;
            end
            levelshown=1;
            % compose wavelet image from parts:
            % obj.C.dec: order: 'a', 'v','h','d'
            % wkeep(d(:)',len);
            lvlup=obj.level-lev+1;
            A= wkeep(obj.C.dec{1,1},obj.N);
            V= wkeep(obj.C.dec{(lvlup-1)*3+2,1},obj.N);
            H= wkeep(obj.C.dec{(lvlup-1)*3+3,1},obj.N);
            D= wkeep(obj.C.dec{(lvlup-1)*3+4,1},obj.N);
            img_dwt=[A,V;
                H,D];
        end
        
        function Cvec=C2vec(obj,cC)
            % convert transformation result to vector form
            % default assumes obj.C is matrix
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            Cvec=[];
            if nargin <2
            	for j=1:length(obj.C.dec )
                   Cvec = [Cvec;obj.C.dec{j}(:)];
            	end
            else
            	for j=1:length(cC.dec )
	           Cvec = [Cvec;cC.dec{j}(:)];
                end
            end
        end
        
        
    end
    
    %% filters
    methods
        
        
        function [signal2, thresh, fCoeff]=hardthresh1(obj,thresh)
            % [signal2,thresh.fCoeff]=~([thresh]) hard threshold transform coefficients  
            % redefined from class FrameTrafo.            
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin < 2 || isempty(thresh)
                obj.estimateNoise();
                thresh=4*obj.sigma;
            end
            fCoeff=obj.C; % copy by value needed (assignment works if obj.C is array)
            
            % filter
            for j=1:length(obj.C.dec)
                fCoeff.dec{j}(abs(fCoeff.dec{j})<thresh)=0;
            end
            % reconstructed signal
            yn = obj.rec(fCoeff) ;
            % 			yn = yn-min(yn(:));
            % 			yn = yn/max(yn(:));
            signal2=obj.ts.clone();
            signal2.replace_signal(yn);
            signal2.signalname=[obj.ts.signalname,', hard thresh. (',obj.basisname,...
                ', ',num2str(thresh,'%3.1e'),')'];
            signal2.colormap_active=obj.ts.colormap_active;
            
        end
        
        function [signal2, thresh, fCoeff]=softthresh1(obj,thresh)
            % [yn,thresh]=~([thresh]) denoise using coefficient thresholding
            % apply e.g. to a test signal created by w= Signal2D.make_sinus(5,0.1).
            % cf. hardthresh2 using 2 thresholds for different kinds of
            % coefficients.
            % hardthresh1 is used by sparseApprox; subclasses redefine hardthresh
            % if necessary.
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin < 2 || isempty(thresh)
                obj.estimateNoise();
                thresh=4*obj.sigma;
            end
            fCoeff=obj.C; % copy by value needed (assignment works if obj.C is array)
            %
            for j=1:length(obj.C.dec)
                fCoeff.dec{j}(abs(fCoeff.dec{j})<thresh)=0;
                % simple soft filter: reduce all coeff. by size
                fCoeff.dec{j}=fCoeff.dec{j}-sign(fCoeff.dec{j})*thresh;
            end
            
            % reconstructed signal
            yn = obj.rec(fCoeff) ;
            
            % renormalisation prevents calculation of PSNR, hence commented
            % out:
            %           yn = yn-min(yn(:));
            % 			yn = yn/max(yn(:));
            signal2=obj.ts.make_like();
            signal2.replace_signal(yn);
            signal2.signalname=[obj.ts.signalname,', soft thresh. (',obj.basisname,...
                ', ',num2str(thresh,'%3.1e'),')'];
            signal2.colormap_active=obj.ts.colormap_active;
            
        end
        
        function [signal2, thresh1, thresh2,fCoeff]=hardthresh2(obj,thresh1, thresh2)
            % [yn,thresh1, thresh2]=~(thresh1, thresh2) 2 level wavelet thresholding
            % apply e.g. to a test signal created by w= Signal2D.make_sinus(5,0.1).
            % use ndwt2 and indwt2 (undecimated WT) to ensure shift invariance
            % [yn,thresh]=~(thresh) denoise using wavelet thresholding
            % apply e.g. to a test signal created by w= Signal2D.make_sinus(5,0.1).
            % use ndwt2 and indwt2 (undecimated WT) to ensure shift invariance
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if  nargin<2 || isempty(thresh1) || nargin <3
                obj.estimateNoise;
            end
            if nargin<2 || isempty(thresh1)
                thresh1=4*obj.sigma;
            end
            if nargin < 3
                thresh2=3*obj.sigma;
            end
            assert(thresh1>=thresh2,'thresh1 is for finest details, thresh2 for coarsest');
            
            fCoeff=obj.C; % copy by value needed (assignment works if obj.C is array)
            
            % filter trend signal with thresh2
            fCoeff.dec{1}(abs(fCoeff.dec{1})<thresh2)=0;
            
            % filter detail signals with thresh1
            for j=2:length(obj.C.dec)
                fCoeff.dec{j}(abs(fCoeff.dec{j})<thresh1)=0;
            end
            % reconstructed signal
            yn = obj.rec(fCoeff) ;
            % renormalisation prevents calculation of PSNR, hence commented
            % out:
            %           yn = yn-min(yn(:));
            % 			yn = yn/max(yn(:));
            signal2=obj.ts.make_like();
            signal2.replace_signal(yn);
            signal2.signalname=[obj.ts.signalname,', 2 hard thresh. (',obj.basisname,...
                ', ',vec2str([thresh1,thresh2],'%3.1e'),')'];
            signal2.colormap_active=obj.ts.colormap_active;
            
        end
        
        function [signal2,BSE,s,ortho_defect,fCoeff]=sparseApprox(obj,cs,p)
            % [signal2,BSE,s]=~(c_s,p) best s-sterm approx. at compression factor c_s
            % using curvelet decomposition.
            % BSE ... error of best s-term approximation (sparsity defect)
            % s ... sparsity of filtered curvelet trafo
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(cs>1,'local','compression rate >1');
            if nargin < 2 || isempty(cs)
                cs=16;
            end
            if nargin <3 || isempty(p)
                p=2;  % lp-norm
            end
            
            % one cannot just take the 1/cs largest coefficients of
            % the undecimated WT, because that would leave only the
            % coarsest levels. Instead compute threshold from
            % the decimated WT.
            % Of course this will not achieve the required compression
            % rate, it rather elucidates the distorsions of the
            % decimated wavelet trafo.
            % Alternatively, one could average the levels of the UWT
            % to construct an averaged DWT, threshold the averaged DWT
            % and reconstruct with the inverse of DWT.
            dw= obj.copy2decimated();
            dw.dec;
            cC=dw.C;
            sC=sort(abs(cC(:)),'descend');
            obj.check(sC,'all(local>=0)','sC is non-neg. vector');
            
            % threshold is chosen s.t. 1/cs of all elements are bigger
            thresh=sC(ceil(numel(cC)/cs));            
            % apply thershold and reconstruct:
            if nargout<5
                signal2=obj.hardthresh1(thresh);
            else % return also filtered coefficients
                [signal2,~,fCoeff]= obj.hardthresh1(thresh);
            end                               
            
            if nargout>=2   % ||W-F(W)||/||W||
                BSE=(sum(sC(sC<thresh)).^p)^(1/p);   % needs original values
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
        
    end
    
    %% statics
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Wavelet2U_mlab();
        end
        
    end
    
    %% graphics
    methods
        
        
    end
    
end



