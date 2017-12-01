classdef FrameTrafoPair < FrameTrafo
    % pair of Parseval frame transformations
    % GT Stuttgart, 2014
    %
    properties
        
        ftrafo1   %@<FrameTrafo>  parseval frame transform 1
        ftrafo2   %@<FrameTrafo>  parseval frame transform 2                
        
    end
    
    properties (SetAccess=protected)
        synweight1  %<real> weight for synthesis of ftrafo1, synweight2=1-synweight1
        anweight1 %<real> weight for analysis coeff. of ftrafo1, controlling contribution to 1-norm
    end
    
    %% commands and constructors
    methods
        
        function obj=FrameTrafoPair(signal)
            % constructor obj=~(signal)
            assert(nargin<1 || isa(signal,'SignalClass'),'needs SignalClass');
            if nargin<1
                signal=[];
            end
            obj = obj@FrameTrafo(signal);
            % each frame gets the same weight:
            obj.synweight1=0.5;  
            % no rel. amplification or attenuation of either type of decomposition coefficients
            obj.anweight1=0.5;
        end
        
        function set_algo(obj)
            set_algo@FrameTrafo(obj);
            obj.algo.name='frame pair';
        end
                
        function set_transforms(obj, hftrafo1, hftrafo2)
            % ~(hftrafo1, hftrafo2) set 2 decoder transformations
            obj.requireOnly(isempty(hftrafo1) || isa(hftrafo1,'FrameTrafo'),...
                'local','needs 2 FrameTrafos');
            obj.requireOnly(isempty(hftrafo2) || isa(hftrafo2,'FrameTrafo'),...
                'local','needs 2 FrameTrafos');            
            obj.ftrafo1=hftrafo1;
            obj.ftrafo2=hftrafo2;
            obj.reset_Trafo();
            
            obj.colormap_active='default';            
        end
        
        function set_signal(obj,hts)
            % ~(hts) set time signal to hts of class SignalClass
            obj.ftrafo1.set_signal(hts);
            obj.ftrafo2.set_signal(hts);
            obj.ts=hts;
            obj.reset_Trafo;
            % do not check invariant because caller may be constructor
            obj.ensure(~obj.isTrafodone,'local', ...
                'decomposition is reset');
        end
        
        function set_synweight1(obj,fw1)
           % ~(fw1) set weight for frame 1, e.g. 0.5 (frame 2 gets weight 1-fw1). 
            obj.requireOnly(fw1>=0 && fw1<=1,'local', 'weight is >=0 and <=1');
            obj.synweight1 =fw1;            
        end
        
        function set_anweight1(obj,fw1)
           % ~(fw1) set weight for frame 1, e.g. 0.5 (frame 2 gets weight 1-fw1). 
            obj.requireOnly(fw1>0 && fw1<1,'local', 'weight is > 0 and < 1');
            obj.anweight1 =fw1;      
            obj.ts=[];  % to detect change of anweight1
        end
        
        function set_C(obj,cC)
            %~(cC) set transform coefficients ( to save or test recon)
            obj.requireOnly(isnumeric(cC),'local','matrix or vector needed');
            L1=obj.ftrafo1.frame_length;
            obj.C=cC;
            obj.ftrafo1.set_C(cC(1:L1));
            obj.ftrafo2.set_C(cC(L1+1:end));
        end
        
        function reset_Trafo(obj)
            % reset trafo (e.g. after signal has changed)
            if ~isempty(obj.ftrafo1)
                obj.ftrafo1.C=[];
            end
            if ~isempty(obj.ftrafo2)
                obj.ftrafo2.C=[];
            end
        end
        
        function dec(obj)
            % decomposition (analysis) yields frame coeffs. obj.C
            obj.requireOnly(~isempty(obj.ftrafo1)&&~isempty(obj.ftrafo2),'local', 'ftrafos set');
            obj.ftrafo1.dec;
            obj.ftrafo2.dec;
            obj.C=[obj.anweight1*obj.ftrafo1.C2vec;(1-obj.anweight1)*obj.ftrafo2.C2vec];
            
        end
        
        
    end
    
    %% queries
    methods
        
        function ok=isTrafodone(obj)
            ok=obj.ftrafo1.isTrafodone &&  obj.ftrafo2.isTrafodone;
        end
        
        function yn=rec(obj,cC)
            % signal reconstruction from frame coefficients C
            % here we assume that both frames are Parseval,
            % hence their union is a tight frame with frame constant 2;
            % a weighted average will again produce a Parseval frame.
            % x=sum(<x,x_i>x_i, x= sum(<x,y_i>y_i) =>
            % x=w*sum(<x,x_i>x_i+(1-w)*sum(<x,y_i>y_i)
            obj.requireOnly(~isempty(obj.ftrafo1)&&~isempty(obj.ftrafo2),'local', 'ftrafos set');
            L1=obj.ftrafo1.frame_length;
            if obj.synweight1~=1  % needs both transforms:
                yn=obj.synweight1/obj.anweight1*obj.ftrafo1.synthesize(cC(1:L1)) + ...
                    (1-obj.synweight1)/(1-obj.anweight1)*obj.ftrafo2.synthesize(cC(L1+1:end));
            else % needs only 1st transform:
                yn=1/obj.anweight1*obj.ftrafo1.synthesize(cC(1:L1));
            end
                
        end
        
        function xn= analyze(obj,cC)
            % yn=~(x) decompose x 
            % same as dec but redefinitions must return a matrix of
            % which can be used as input of synthesize.
            xn=[reshape(obj.anweight1 * obj.ftrafo1.analyze(cC),[],1);...
                reshape((1-obj.anweight1) * obj.ftrafo2.analyze(cC),[],1)];
        end
        
        function yn= synthesize(obj,x)
            % yn=~(x) decompose x 
            % here we assume that both frames are Parseval,
            % hence their union is a tight frame with frame constant 2.
            L1=obj.ftrafo1.frame_length;
            if obj.synweight1~=1  % needs both transforms:
                yn=obj.synweight1/obj.anweight1*obj.ftrafo1.synthesize(x(1:L1))+...
                    (1-obj.synweight1)/(1-obj.anweight1)*obj.ftrafo2.synthesize(x(L1+1:end));
            else % needs only 1st transform:
                yn=1/obj.anweight1*obj.ftrafo1.synthesize(x(1:L1));
            end
        end
        
        function d=frame_dim(obj)
            d=obj.ftrafo1.frame_dim();
        end
        
        function d=frame_length(obj)
            d=obj.ftrafo1.frame_length()+obj.ftrafo2.frame_length();
        end
        
        function bn= basisname(obj)
            % bn=~() basis name
            add1=[]; add2=[];
            if obj.isTrafodone
                %                 add1=['[',num2str(obj.ftrafo1.frame_length()),']'];
                %                 add2=['[',num2str(obj.ftrafo2.frame_length()),']'];
                add1=['[',num2str(obj.anweight1),']'];
                add2=['[',num2str(1-obj.anweight1),']'];
            end
            bn=[obj.ftrafo1.basisname,add1,' & ',obj.ftrafo2.basisname,add2];
        end
        
        
        function vC=C2vec(obj,cC)
            % convert frame coefficients to vector form
            obj.requireOnly(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin<2
            	vC = [obj.anweight1*obj.ftrafo1.C2vec;(1-obj.anweight1)*obj.ftrafo2.C2vec];
            else
            	vC = [obj.anweight1*obj.ftrafo1.C2vec(cC);(1-obj.anweight1)*obj.ftrafo2.C2vec(cC)];
            end
            obj.ensureOnly(isvector(vC),'local','result is column vector');
        end
        
        function nc=numelC(obj,cC)
            % nc=~([cC]) number of elements in decomposition
            % must be redefined for those subclasses, where C is not a
            % matrix but e.g. a cell array.
            if nargin <2
                nc= obj.ftrafo1.numelC + obj.ftrafo2numelC;
            else
                nc=numel(obj.C2vec(cC));
            end
        end
        
        function mat=C2graph(obj)
            % mat=~() convert coefficient vector to graphical output form
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            mat1=obj.ftrafo1.C2graph();
            mat2=obj.ftrafo2.C2graph();
            if isvector(mat1)
                mat1=reshape(mat1,1,[]);
            end
            if isvector(mat2)
                mat2=reshape(mat2,1,[]);
            end
            sdiff=size(mat1,1)-size(mat2,1);
            if sdiff>0
                mat=[mat1, padarray(mat2,[sdiff,0],NaN,'post')];
            elseif sdiff <0
                mat=[padarray(mat1,[-sdiff,0],NaN,'post'),mat2];
            else
                mat=[mat1,mat2];
            end
        end                
        
    end
    
    %% filters
    methods
        
         function [signal2,BSE,s,ortho_defect]=sparseApprox(obj,cs,p)
            % [signal2,BSE,s]=~(c_s,p) best p-norm s-term approx. at compression c_s
            % signal2 ... approximation of signal obj.ts by best s-terms.
            % BSE ... lp-error of best s-term approximation (sparsity defect)
            % s ... sparsity of filtered transform
            
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            obj.requireOnly(cs>1,'local','compression rate >1');
            if nargin < 2 || isempty(cs)
                cs=16;
            end
            if nargin <3 || isempty(p)
                p=2;  % lp-norm
            end
            
            % obj.C is already vector with analysis weights applied:            
            sC=sort(abs(obj.C),'descend');
            obj.check(sC,'all(local>=0)','sC is non-neg. vector');
            % threshold is chosen s.t. 1/cs of all elements are bigger
            thresh=sC(ceil(numel(sC)/cs));
            % filter coefficientw
            cC=obj.C;
            cC(abs(cC)<thresh)=0;
            
            signal2=obj.ts.make_like();
            signal2.set_signal(reshape(obj.rec(cC),obj.ts.N));
            
            signal2.signalname=[obj.ts.signalname,', hard thresh. (',obj.basisname,...
                ', c=',num2str(cs,'%3.1f'),')'];
            
            if nargout>=2   % ||W-F(W)||/||W||
                BSE=(sum((sC(sC<thresh)).^p))^(1/p);   % needs original values
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
    
    %% graphics
    methods
        
        function graph_trafo(obj, open_new, whichpart)
            % ~([open_new,whichpart) show transform coefficients
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin <2
                open_new=true;
            end
            if nargin<3
                whichpart=[];
            end
            if open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            if whichpart==1
                obj.ftrafo1.graph_trafo(false);
            elseif whichpart==2
                obj.ftrafo2.graph_trafo(false);
            else
                s2=obj.ts.make();
                s2.set_signal(obj.C2graph());
                s2.signalname={['coeff. in frame ',obj.basisname], obj.ts.signalname};
                s2.graph_signal(false);
            end
            
        end
        
        function graph_distribution(obj, open_new, whichpart)
            % ~([open_new,whichpart) show distribution of transform coefficients in open window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            if nargin <2
                open_new=true;
            end
            if nargin<3
                whichpart=1;
            end
            if open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            if whichpart==1
                obj.ftrafo1.graph_distribution(false);
            else
                obj.ftrafo2.graph_distribution(false);
            end
            
        end
        
        function show_recon(obj)
            % ~() show signal reconstructed from transform in new window
            obj.require(obj.isTrafodone,'local',obj.msgDoTrafo);
            prepfigure(obj.fig,obj.algo,obj.figopt);
            
            subplot(2,2,1);
            obj.ts.graph_signal(false);
            
            subplot(2,2,2);
            obj.graph_recon(false)
            
            subplot(2,2,3);
            obj.ftrafo1.graph_recon(false);
            
            subplot(2,2,4);
            obj.ftrafo2.graph_recon(false);
            
        end
        
        
    end
    
end

