classdef Wavelet3D_wlab<MultiResTransform3
    % 3d wavelet transform using WaveLab
    %  **************  set pathes by calling WavePath; *********
    %
    % Example:
    % =========
    % Example:
    % --- 3d sparsification:
    % signal2=Signal2D.make_fromImage('cameraman.bmp');
    % signal2.resize(128);
    % L=128; signal3=Signal3D.make_CyclefromLowerDimSig(signal2,1,L);    
    % crop to dyadic and cubic:
    % signal3.crop();
    % w3= Wavelet3D_wlab(signal3);
    % w3.set_basisname('db2');    
    % w3.dec;
    % c=33; test3=w3.sparseApprox(c);
    % test3.graph_signal;
    %
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
    
    
    %% commands
    methods
        
        function obj=Wavelet3D_wlab(signal)
            % constructor
            global WAVELABPATH
            if nargin==0
                signal=Signal3D();
            end
            obj = obj@MultiResTransform3(signal);
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
            detail_labels={'HLL','LHL','HHL','LLH','HLH','LHH','HHH'};
            detail_idx=num2cell(1:7);
            obj.detail_str2idx=containers.Map(detail_labels,detail_idx);
            obj.detail_idx2str=containers.Map(detail_idx,detail_labels);
            
            obj.ensure(isa(obj.ts,'Signal3D'),'local','3D signal created');
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
            set_basisname@MultiResTransform3(obj, wvn);
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
            set_algo@MultiResTransform3(obj);
            obj.algo.name='DWT';   % decimated wavelet trafo
            obj.algo.toolbox='wavelab';
        end
        
        
        function dec(obj)
            % discrete multiresolution analysis
            obj.requireOnly(~isempty(obj.sdata),'local','non-empty sample size');
            obj.requireOnly(~isempty(obj.qmf),'local','needs filter');
            obj.lev=obj.deepestlev;
            obj.C = FWT3_PO(obj.sdata,obj.level_bottomUp,obj.qmf);
            
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
            yn= IWT3_PO(wc,obj.level_bottomUp,obj.qmf);
            
        end
        
    end
    
    
    %% queries
    
    methods
        
        function nv=frame_norm(obj)
            % l2-norm of transform coefficients
            % tests in this form only isometries (e.g. ONS)
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            nv=1;  % norm(obj.C,2)/obj.ts.norm(2) will not work in matlab toolbox for higher dbx
        end
        
        
        function L=level(obj)
            L= obj.lev;
        end
        
        function L=level_bottomUp(obj)
            L=min(obj.size_dyadic)-obj.level();
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x 
            obj.lev=obj.deepestlev;
            yn=FWT3_PO(reshape(x,obj.N),obj.level_bottomUp,obj.qmf);
        end
        
        function xn= synthesize(obj,y)
            % xn=~(x) synthesize
            if isempty(obj.lev) % synthesize might be called before analyze
                obj.lev=obj.deepestlev;
            end
            xn= IWT3_PO(reshape(y,obj.N),obj.level_bottomUp,obj.qmf);
        end
        
        function n=wvsupport(obj)
            % length of support of wavelet
            n=length(obj.qmf);
        end
        
        
        
    end
    
    %% filter
    
    methods
        
        
        
        
    end
    
    %% transforms
    methods                
        
        function coeff=detcoef(obj,lvl,type)
            % coeff=~(lvl,type) extract detail coefficients from wavelet transform
            % type ='d' for diagonal fluctutation.
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(lvl>=1 && lvl<=obj.level,'local',...
                'lvl is below computed level');
            obj.requireOnly(ischar(type) || (isnumeric(type) && type>=1 && type <=7),'local',...
                'detail type is number in 1..7 or type string');
            if ischar(type)
                type=obj.detail_str2idx(type);
            end
            d=round(log2(obj.N));
            ls=2.^(d-lvl+1);
            % convert type to a binary number with leading zeros
            b=arrayfun(@str2num,dec2bin(type));
            L=3-length(b);
            b=padarray(b(:),L,0,'pre');
            coeff=obj.C(1+b(1)*ls(1)/2:ls(1)/2+b(1)*ls(1)/2, ...
                        1+b(2)*ls(2)/2:ls(2)/2+b(2)*ls(2)/2, ...
                        1+b(3)*ls(3)/2:ls(3)/2+b(3)*ls(3)/2);        
        end
        
        function A= appcoef(obj)
            % coeff=~() extract approximate coefficients from wavelet transform
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            d=round(log2(obj.N));
            A=obj.C(1:2^(d(1)-obj.level),1:2^(d(2)-obj.level),1:2^(d(3)-obj.level));
        end
        
        
    end
    
    
    %% statics
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Wavelet3D_wlab();
        end
        
    end
    
    %% graphics
    
    methods
        
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='signal is 3d';
            ok= isa(obj.ts,'Signal3D');
            if ok
                descr='forms dyadic cube';
                ok=obj.ts.isdyadicHypercube();
            end
        end
    end
    
    
    
    
end

