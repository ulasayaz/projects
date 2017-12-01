classdef Wavelet3D_mlab<MultiResTransform3 
    % 3d wavelet transform using the wavelet-toolbox of matlab
    %
    % Example:
    % ==============
    % signal2=Signal2D.make_fromImage('cameraman.bmp');
    % L=32; signal3=Signal3D.make_CyclefromLowerDimSig(signal2,1,L);
    % w3= Wavelet3D_mlab(signal3);
    % w3.set_basisname('db1'); 
    % w3.set_deepestlevel(1);
    % w3.dec;
    % w3.graph_trafo;
    % w3.show_resultd;
    %
    % --- test quality of sparse approximation of 3d signal:
    % c=33; test3=w3.sparseApprox(c);
    % test3.play_signal;
    %
    
    properties
       
    end
    
    properties (SetAccess=protected)
       
    end
    
    
    %% commands
    methods
        
        function obj=Wavelet3D_mlab(signal)
            % constructor            
            if nargin==0
                signal=Signal3D();
            end
			obj = obj@MultiResTransform3(signal);
            obj.requireOnly(license('test','Wavelet_Toolbox'),'local','needs wavelet toolbox');
            
            obj.set_basisname('db1');
            obj.dwtmode_active=dwtmode('sym','nodisp');
            
            % detail signal:
            detail_labels={'HLL','LHL','HHL','LLH','HLH','LHH','HHH'};
            detail_idx=num2cell(1:7);
            obj.detail_str2idx=containers.Map(detail_labels,detail_idx);
            obj.detail_idx2str=containers.Map(detail_idx,detail_labels);
            
            obj.ensure(isa(obj.ts,'Signal3D'),'local','3D signal created');
        end
        
        function set_dwtmode(obj, modestr)
            % ~(modestr) border extension mode cf. help dwtmode
            % e.g. 'sym' (symmetric extension), 'zpd' (zero-padding),
            % 'ppd' periodic, 'sp1 ' smooth padding
            obj.dwtmode_active=modestr;
            dwtmode(modestr,'nodisp');
        end
        
        
        function set_algo(obj)
            set_algo@MultiResTransform3(obj);
            obj.algo.name='DWT'; % decimated wavelet trafo
            obj.algo.toolbox='Wavelet\_Toolbox';
            
        end
                            
        
        function dec(obj,lev)
            % discrete multiresolution analysis
            obj.require(obj.ts.numel>0,'local','non-empty sample size');
            if nargin <2
                lev=obj.deepestlevAnisotropic();
            end
            % wavelet decomposition from matlab's wavelet toolbox
            obj.C = wavedec3(obj.sdata,lev,obj.basisname);
            
            obj.ensureOnly(obj.isTrafodone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wc)
            %yn=~(cC) % signal reconstruction from wavelet decomp wc (same
            %type as result of dec)
            % wavelet recon from matlab's wavelet toolbox
            yn = waverec3(wc) ;
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
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            L= obj.C.level;
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x, which may be matrix or vector
            yn=obj.C2vec(wavedec3(reshape(x,obj.N),obj.deepestlevAnisotropic(),obj.basisname));
        end
        
        function xn= synthesize(obj,y)
            % xn=~(y) synthesize from decomp y, which is a matrix or vector
            xn=waverec3(obj.vec2C(y(:))) ;
        end
        
        
    end
    
    %% filter
    
    methods
        
        function hardthresh1_inplace(obj,thresh)
            % filter all coefficients obj.C
            % operation on place to save memory
            for j=1:length(obj.C.dec )
                obj.C.dec{j}(abs(obj.C.dec{j})<thresh)=0;
            end
        end
        
        function hardthresh2_inplace(obj,thresh1, thresh2)
            % filter with 2 thresholds
            obj.requireOnly(thresh1<=thresh2,'local', 'thresh1<=thresh2');
            % filter all with thresh2
            L=length(obj.C.dec );
            for j=1:L
                obj.C.dec{j}(abs(obj.C.dec{j})<thresh2)=0;
            end
            % filter finest detail signal with thresh1
            obj.C.dec{L}(abs(obj.C.dec{L})<thresh1)=0;
        end
        
        function softthresh1_inplace(obj,thresh)
            % simple soft filter: reduce all coeff. obj.C by size
            % operation on place to save memory
            for j=1:length(obj.C.dec )
                obj.C.dec{j}(abs(obj.C.dec{j})<thresh)=0;
                % now also reduce other coefficients
                obj.C.dec{j}=obj.C.dec{j}-sign(obj.C.dec{j})*thresh;
            end
        end
        
        
    end
    
    %% transforms
    methods
        
        function coeff=detcoef(obj,lvl,type)
            % coeff=~(level[,type]) extract detail coefficients from wavelet transform
            % e.g. type ='HHH' for high freq. filter in all dim;
            % if type is missing, it collects all filters of level lvl.
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(lvl>=1 && lvl<=obj.level,'local',...
                'level is below computed level');
            obj.requireOnly(nargin <3 || ischar(type) || type<=7,'local',...
                'type is number or label.');            
            if ischar(type)
                type=obj.detail_str2idx(type);
            end
            dlvl=obj.level;
            if nargin >=3
                idx=7*(dlvl-lvl)+1+type;
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
        
        function Cvec=C2vec(obj,cC)
            % convert transformation result to vector form
            % default assumes obj.C is matrix
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            Cvec=[];
            if nargin<2
                for j=1:length(obj.C.dec )
                    Cvec = [Cvec;obj.C.dec{j}(:)];
                end
            else
                for j=1:length(cC.dec )
                    Cvec = [Cvec;cC.dec{j}(:)];
                end                
            end
        end
        
        function nc=numelC(obj,cC)
            % nc=~([cC]) number of elements in decomposition
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            nc=0;
            if nargin<2
                for j=1:length(obj.C.dec )
                    nc = nc +numel(obj.C.dec{j});
                end
            else
                for j=1:length(cC.dec )
                    nc = nc+numel(cC.dec{j});
                end
            end
        end
        
        function cC=vec2C(obj,vC)
            % convert vector to form usable for reconstruction rec
            cC=obj.C;
            offset=0;
            for j=1:length(obj.C.dec )
                L=numel(obj.C.dec{j});
                ss=size(obj.C.dec{j});
                cC.dec{j}=reshape(vC(offset+1:offset+L),ss);
                offset=offset+L;
            end
        end
        
        
    end
    
    
    %% statics
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Wavelet3D_mlab();
        end
        
    end
    
    %% graphics
    
    methods
               
        
    end
    
    %% invariant
    
    
    
end

