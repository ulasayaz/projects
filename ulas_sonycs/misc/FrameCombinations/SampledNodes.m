classdef SampledNodes < DC
    % nodes to be sampled
    % GT Stuttgart, 2014
    %
    % example:
    % -- set signal:
    % signal=Signal2D.make_fromImage('cameraman.bmp');
    % -- set frame transform:
    % B=WalshHadamard2(signal); param=struct;
    % -- or alternatively:
    % B=Wavelet2D_mlab(signal); param.w1=0.95; B.set_deepestlev(1);
    % -- compute nodes
    % c=8; nodes=B.nodesPDF(c,param);
    % -- show graph of nodes:
    % nodes.graph_mask;
    %
    
    properties (SetAccess=protected)
        data        % nodes as vector or cell
        numnodes
    end
    
    properties
        signalSize  %@<vector> size of signal where nodes are sampled
        Csize       %@<vector> size of transformed signal (using ftrafo)
        ftrafo      %@<FrameTrafo> transform applied to signal before samling
        fL          %@<integer> frame number if same nodes are used for each frame
    end
    
    
    %% constructor and commands
    methods
        
        function obj= SampledNodes(hdata, ssize, hftrafo)
            % constructor obj=~(hdata, ssize,[hftrafo])
            obj.requireOnly(nargin<3 || isa(hftrafo,'FrameTrafo'),'local',...
                'hftrafo is FrameTrafo');
            if nargin <1
                hdata=[];
            end
            if nargin <2
                ssize=[0,0];
            end
            if nargin <3
                hftrafo=[];
            end
            obj.data=hdata;
            obj.signalSize=ssize;
            obj.ftrafo=hftrafo;
            obj.Csize=[];
            obj.set_numnodes();
            obj.fL=1;
        end
        
        function set_data(obj,hdata)
            % ~(hdata) set node data
            obj.data=hdata;
            obj.set_numnodes();
        end
        
        function set_cell(obj,j,hdata)
            % ~(hdata) set node data
            obj.data{j}=hdata;
            obj.numnodes=[];
        end
        
        function crop(obj)
            % ~() crop  so that all cells have the same number of nodes
            L=obj.frameNumber;
            if L>1
                nmin=min(cellfun(@numel, obj.data));
                for j=1:L
                    obj.data{j}=obj.data{j}(1:nmin);
                end
            end
        end
        
        
    end
    
    methods (Hidden)
        
        function set_numnodes(obj)
            if iscell(obj.data)
                n=0;
                obj.fL=1;
                for j=1:length(obj.data)
                    n=n+numel(obj.data{j});
                end
                obj.numnodes=n;
            else
                obj.numnodes=numel(obj.data);
            end
        end
        
    end
    
    
    %% queries
    methods
        
        function nd=get(obj,n)
            % ~([n]) get node data
            if iscell(obj.data)
                nd=obj.data{n};
            else
                nd=obj.data;
            end
        end
        
        function n=numel(obj,n)
            if nargin <2
                if isempty(obj.numnodes)
                    obj.set_numnodes;
                end
                n=obj.fL*obj.numnodes;
            elseif iscell(obj.data)
                n=numel(obj.data{n});
            else
                n=numel(obj.data);
            end
        end
        
        function L=frameNumber(obj)
            % L=~() number of node frames (different or identical)
            if iscell(obj.data)
                L=length(obj.data);
            else
                L=obj.fL;
            end
        end
        
        function L=frameNumberUnique(obj)
            % L=~() number of unique node frames
            if iscell(obj.data)
                L=length(obj.data);
            else
                L=1;
            end
        end
        
        function c=compression(obj)
            % c=~() compression rate card(signal)/card(nodes)
            if ~obj.ftrafo.isTrafodone
                obj.ftrafo.dec;
            end
            c= obj.frameNumber*obj.signalSize(1)*obj.signalSize(2)/obj.numel;
        end
        
        function mask=embedNodesInTrafo(obj,n)
            % ~() computes size of transformed signal
            obj.requireOnly(isa(obj.ftrafo,'FrameTrafo'),'local',...
                'needs transformation');
            if nargin <2
                n=1;
            end
            if ~isempty(obj.ftrafo) && isequal(obj.signalSize(1:2),obj.ftrafo.ts.N)
                if ~obj.ftrafo.isTrafodone
                    obj.ftrafo.dec;
                end
                % embed nodes in transformation vector C:
                cC=false(size(obj.ftrafo.C));
                cC(obj.get(n))=true;
                % transform C to matrix:
                mask=obj.ftrafo.dec2graph(cC);
            else
                % embed nodes directly into signal:
                mask=false(obj.signalSize);
                mask(obj.get(n))=true;
            end
            
        end
        
        function mask= nodes2mask(obj,n)
            % mask=~([n]) returns boolean mask of position of nodes
            % if nodes data are a cell array, the n-th cell is converted.
            obj.requireOnly(~isempty(obj.signalSize) && min(obj.signalSize)>0,...
                'true','needs signal size');
            obj.requireOnly(nargin<2 || n<=length(obj.data),...
                'true','needs signal size');
            if nargin <2
                n=[];
            end
            L=obj.frameNumber;
            
            if L>1
                if ~isempty(n)
                    mask=obj.embedNodesInTrafo(n);
                else
                    mask1=obj.embedNodesInTrafo(1);
                    mask=zeros([size(mask1),L]);
                    mask(:,:,1)=mask1;
                    for j=2:L
                        mask(:,:,j)=obj.embedNodesInTrafo(j);
                    end
                end
            else
                mask=obj.embedNodesInTrafo();
            end
        end
        
        function s=nodes2signal(obj)
            % s=~() returns nodes mask as signal
            mask=obj.nodes2mask;
            L=obj.frameNumber;
            ss=obj.signalSize;
            if L>1 && length(ss)<3
                ss=[ss(:)',L];
            end
            sL=sum(ss>1);
            if sL==1
                s=TimeSignal(squeeze(mask));
                s.marker='.';
                s.linestyle='';
                s.markercolor='b';
            elseif sL==2
                s=Signal2D(mask);
            elseif sL==3
                s=Signal3D(mask);
            end
            if ~isempty(obj.ftrafo)
                tn=obj.ftrafo.algo.name;
            else
                tn='Id';
            end
            sn=['node mask (transform=',tn,')'];
            if L>1
                sn=[sn,', L=',num2str(L)];
            end
            s.signalname=sn;
            obj.ensure(isa(s,'SignalClass'),'local','returns obkject of SignalClass');
        end
        
    end
    
    %% graphics
    methods
        
        function graph_mask(obj, open_new)
            % show signal in the open window
            if nargin <2
                open_new=true;
            end
            if open_new
                prepfigure(1);
            end
            ss=obj.nodes2signal;
            ss.graph_signal(false);
        end
        
    end
    
end

