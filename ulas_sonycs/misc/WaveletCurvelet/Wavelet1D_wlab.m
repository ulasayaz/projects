classdef Wavelet1D_wlab <MultiResTransform1
    % 2-dim. decimated wavelet transformations using the open source toolbox wavelab
    % set pathes by calling WavePath;
    %
    % features:
    %    --- both discrete multi-resolution analysis (.dec)
    %    --- noise filter (.show_noisefilter)
    %
    %
    % Example:
    %====================================
    % --- discrete wavelet trafo (multi-resolution analysis):
    %     w=Wavelet1D_wlab(); w.show_test;
    % --- choose signal (image):
    %     signal1=TimeSignal.make_sinus();
    % --- db2 (not db1) sensitive to discont. in 1st derivative of signal exp_symm
    %     w=Wavelet1D_wlab(signal1);
    %     -- db2,3 show edges at border and along image axis:
    %     w.show_trafos_selection({'db1','db2','db3'});
    % --- noise filter
    %     omega0=5;sig=0.2;hN=512;
    %     w= Wavelet1D_wlab(TimeSignal.make_sinus(omega0,sig,hN));
    %     w.set_basisname('db3');
    %     w.dec;
    %     w.estimateNoise; disp(w.sigma);
    %     w.show_resultd;
    %     w.show_noisefilter;
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
        function obj=Wavelet1D_wlab(signal)
            global WAVELABPATH
            if nargin==0
                signal=TimeSignal();
            end
            obj = obj@MultiResTransform1(); % cannot set signal before filter)
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
            set_basisname@MultiResTransform1(obj, wvn);
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
            set_algo@MultiResTransform1(obj);
            obj.algo.name='DWT';   % decimated wavelet trafo
            obj.algo.toolbox='wavelab';
        end
        
        
        function dec(obj)
            % discrete multiresolution analysis
            obj.requireOnly(obj.ts.numel>0,'local','non-empty sample size');
            obj.requireOnly(~isempty(obj.qmf),'local','needs filter');
            obj.lev=obj.deepestlev;
            if obj.OW_mlab2wlab.isKey(obj.wvfamily) % orhtogonal WT
                obj.C = FWT_PO(obj.sdata,obj.level_bottomUp,obj.qmf);
            elseif  obj.BOW_mlab2wlab.isKey(obj.wvfamily)  % bi-orthogonal WT
                % symmetric border extension improves PSNR by ca. 0.5dB
                obj.C = FWT_SBS(obj.sdata,obj.level_bottomUp,obj.qmf, obj.dqmf);
                %obj.C = FWT_PB(obj.sdata,obj.level_bottomUp,obj.qmf, obj.dqmf);
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
                yn= IWT_PO(wc,obj.level_bottomUp,obj.qmf);
            elseif obj.BOW_mlab2wlab.isKey(obj.wvfamily)  % bi-orthogonal WT
                % symmetric border extension improves PSNR by ca. 0.5dB
                yn= IWT_SBS(wc,obj.level_bottomUp,obj.qmf, obj.dqmf);
                %yn= IWT_PB(wc,obj.level_bottomUp,obj.qmf, obj.dqmf);
            end
            
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x
            obj.lev=obj.deepestlev;
            yn=FWT_PO(reshape(x,obj.N),obj.level_bottomUp,obj.qmf);
        end
        
        function xn= synthesize(obj,y)
            % xn=~(x) synthesize
            if isempty(obj.lev) % synthesize might be called before analyze
                obj.lev=obj.deepestlev;
            end
            xn= IWT_PO(reshape(y,obj.N),obj.level_bottomUp,obj.qmf);
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
        
        function coeff=detcoef(obj,level)
            % coeff=~(level,type) extract detail coefficients from wavelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(level>=1 && level<=obj.level,'local',...
                'level is below computed level');
            d=round(log2(obj.N));
            ls=2.^(d-level+1);
            coeff=obj.C(1:floor(ls(1)/2));
        end
        
        function A= appcoef(obj)
            % coeff=~() extract approximate coefficients from wavelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            d=round(log2(obj.N));
            A=obj.C(1:2^(d-obj.level));
        end
        
        
    end
    
    %% filters
    methods
        
    end
    
    %% tests
    methods
        
    end
    
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Wavelet1D_wlab();
        end
        
        
    end
    %% graphics
    methods
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 1d';
            ok= isa(obj.ts,'TimeSignal');
            if ok
                descr='forms dyadic cube';
                ok=obj.ts.isdyadicHypercube();
            end
        end
    end
end





