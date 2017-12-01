classdef Wavelet2U_wlab <MultiResTransform2
    % 2-dim. unddecimated wavelet transformations using the open source toolbox wavelab
    % features:
    %    --- both discrete multi-resolution analysis (.wavedec)
    %    --- noise filter (.show_noisefilter)
    %    --- needs paths set by WavePath()
    %
    % Example:
    %====================================
    % --- choose signal (image):
    %     signal1=Signal2D.make_exp_symm();
    %     signal2=signal2:
    % --- discrete wavelet trafo (multi-resolution analysis):
    %     w=Wavelet2U_wlab(); w.show_test;
    % --- db2 (not db1) sensitive to discont. in 1st derivative of signal exp_symm
    %     w=Wavelet2U_wlab(signal1);
    %     w.set_dwtmode('sym');
    %     -- db2,3 show edges at border and along image axis:
    %     w.show_trafos_selection({'db1','db2','db3'});
    %     w.fig=2;
    %     w.set_dwtmode('sp1');
    %     -- only edge along image axis remains:
    %     w.show_trafos_selection({'db1','db2','db3'});
    %
    % --- noise filter
    %     omega0=5;sig=0.2;hN=512;
    %     w= Wavelet2U_wlab(Signal2D.make_sinus(omega0,sig,hN));
    %     w.dec;
    %     w.estimateNoise; disp(w.sigma);
    %     w.show_resultd;
    %     w.show_noisefilter;
    %
    % --  photos have many edges where 1st derivative is discontinuous
    %     w= Wavelet2U_wlab(signal2);
    %     w.set_deepestlev(2); w.dec;
    %     w.show_resultd;
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
        function obj=Wavelet2U_wlab(signal)
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
            obj.set_signal(signal);
            
            % default (only?)border extensionof wavelab is periodic:
            obj.dwtmode_active='ppd';
            obj.wcL=2;
            
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
            obj.algo.name='TIWT';   % decimated wavelet trafo
            obj.algo.toolbox='wavelab';
        end
        
        
        function dec(obj, lvl)
            % discrete multiresolution analysis
            obj.requireOnly(~isempty(obj.img),'local','non-empty sample size');
            obj.requireOnly(~isempty(obj.qmf),'local','needs filter');
            if nargin <2
                lvl=obj.deepestlev();
            end
            obj.lev=lvl;
            if obj.OW_mlab2wlab.isKey(obj.wvfamily) % orthogonal WT
                obj.C = FWT2_TI(obj.img,obj.level_bottomUp,obj.qmf);
            elseif  obj.BOW_mlab2wlab.isKey(obj.wvfamily)  % bi-orthogonal WT
                % symmetric border extension improves PSNR by ca. 0.5dB
                error('biorthogonal undecimated WT not yet implemented');
            end
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wc)
            %yn=~(cC) % signal reconstruction from wavelet decomp wc
            obj.requireOnly(~isempty(obj.qmf),'local','needs filter');
            if isempty(obj.lev) % rec might be called before dec
                obj.lev=obj.deepestlev;
            end
            if obj.OW_mlab2wlab.isKey(obj.wvfamily) % orhtogonal WT
                yn= IWT2_TI(wc,obj.level_bottomUp,obj.qmf);
            elseif obj.BOW_mlab2wlab.isKey(obj.wvfamily)  % bi-orthogonal WT
                % symmetric border extension improves PSNR by ca. 0.5dB
                error('biorthogonal undecimated WT not yet implemented');
            end
            
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
            % transform dwt to image (for UWT show 4 parts)
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || lev<=obj.level,'local', 'smaller than max. #levels');
            if nargin < 2
                lev=1;
            end
            levelshown=1;
            A= obj.appcoef();
            V= obj.detcoef(lev,'v');
            H= obj.detcoef(lev,'h');
            D= obj.detcoef(lev,'d');
            img_dwt=[A,V;
                H,D];
        end
        
        
        function coeff=detcoef(obj,lvl,type)
            % coeff=~(level[,type]) extract detail coefficients from wavelet transform
            % e.g. type ='d' for diagonal fluctutation;
            % if type is missing, it collects all filters of level lev.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(lvl>=1 && lvl<=obj.level,'local',...
                'level is number below computed level');
            obj.requireOnly(nargin <3 || ischar(type) || type<=3,'local',...
                'type is number or "h" or "v" or "d".');            
            if nargin>=3
                if isnumeric(type)
                    type=obj.detail_idx2str(type);
                end
                [N,J]=quadlength(obj.ts.xn);
                switch type
                    case 'v'
                        offset=0;
                    case 'h'
                        offset=1;
                    case 'd'
                        offset=2;
                end
                %coeff=obj.C(end-(3*lvl*N-1+offset),:);
                n=(3*(lvl-1)+offset)*N;
                coeff=obj.C(n+1:n+N(1),:);
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
            % coeff=~() extract approximate coefficients from wavelet
            % transform%
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            [N,J]=quadlength(obj.ts.xn);
            A=obj.C(end-N+1:end,:);
        end
        
        function dw=copy2decimated(obj)
            % dw=~() make copy as decimated wavelet trafo
            dw=Wavelet2D_wlab(obj.ts);
            dw.set_basisname(obj.basisname);
            dw.dwtmode_active=obj.dwtmode_active;
            dw.wcL=obj.wcL;
        end
        
    end
    
    %% filters
    methods
        
        function [yn, thresh1, thresh2]=hardthresh2(obj,thresh1, thresh2)
            % [yn,thresh1, thresh2]=~(thresh1, thresh2) 2 level wavelet thresholding
            % apply e.g. to a test signal created by w= Signal2D.make_sinus(5,0.1).
            % use ndwt2 and indwt2 (undecimated WT) to ensure shift invariance
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if  nargin<1 || isempty(thresh1) || nargin <3
                obj.estimateNoise;
            end
            if nargin<1 || isempty(thresh1)
                thresh1=4*obj.sigma;
            end
            if nargin < 3
                thresh2=3*obj.sigma;
            end
            obj.assert(thresh1>=thresh2,'thresh1 is for finest details, thresh2 for coarsest');
            
            uwc=obj.C;
            [N,J]=quadlength(obj.ts.xn);
            % filter all (also coarsest level?)
            uwc(1:3*N,:)       = HardThresh(uwc(1:3*N,:),thresh1);       % Finest.
            uwc(3*N+1:end-N,:) = HardThresh(uwc(3*N+1:end-N,:),thresh2); % Other scales.
            
            % reconstructed signal
            yn = obj.rec(uwc) ;
            
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
            % instead of the following lines we call hardthresh1
            % to faciliated redefintions in subclasses:
            % cC(abs(cC)<thresh)=0; % filter -> cC has changed!
            % yn = obj.rec(cC) ; % reconstruct
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
    
    %% tests
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Wavelet2U_wlab();
        end
        
        function w=show_test(w)
            % w=~(ts) testing the wavelet trafo with signal ts
            assert(nargin <1 || isa(w,'Wavelet2U_wlab'),'w belongs to wavelet class');
            if nargin ==0
                hN=256; v=8;
                hts=Signal2D.make_star(v,hN);
                w=Wavelet2U_wlab(hts);
                w.figopt.pixsizeX=1000;
            end
            w.set_deepestlev(2);
            w.show_trafos_selection({'db1','db2','db3'});
            
        end
        
    end
    %% graphics
    methods
        
        
    end
    
end



