classdef HMatrixSparse < HMatrix
    % sparse matrix implemented as handle (i.e. assignment does not double
    % required memory).
    % Example real
    % sm= SparseMatrix(0.1*sprandn(1000,1000,0.05));    
    % sm.item(2,3)
    % sm.put(2,3,1.1);
    % sm.item(2,3)
    % sm.show;
    % -- complex case:
    % sm= SparseMatrix(0.1*sprandn(1000,1000,0.025)+1i*0.2*sprandn(1000,1000,0.025));
    %
    
    properties               
             
    end
    
    
    % constructor & commands
    methods
        
        function obj=SparseMatrix(hdata,u)
            % ~([v],[u]) constructor:hdata can be a sparse matrix or size,
            % u (string) are units
            obj.requireOnly(nargin<1 || issparse(hdata) || length(hdata)==2,...
                'local','if specified data must be sparse matrix or size information');
            obj.requireOnly(nargin<2 || ischar(u),'local','if specified units must be a string');
            if nargin==0
                % needed to allow empty calls from subclasses
                % and prevent their premature testing of invariant
                ssize=[100,100];
                hdata=[];
            elseif ~issparse(hdata)
                ssize=hdata;
                hdata=[];
            end
            
            %set units
            if nargin>=2
                obj.units=u;
            end
            
            if isempty(hdata)
                obj.data=sparse(ssize(1),ssize(2));
            else
                obj.data=hdata;
            end
            
            obj.fig=1;
            obj.figopt=struct;
            obj.descr='sparse matrix';           
            
            if nargin>0
                % prevent premature testing of invariant by subclasses
                obj.ensure([],'true','verify invariant');
            end
        end        
        
        function sparse(obj,i,j,n)
             % ~(i,j,n) sets obj.data to a sparse matrix created by vectors
             % i,j,n, s.t. obj.data(i(k),j(k)) = s(k)
             obj.requireOnly(nargin==4 && length(i)==length(j) && length(j)==length(n),...
                 'local','3 vectors of the same length');
             obj.data=sparse(i,j,n);                        
        end
        
    end
    
    %% queries
    methods
        
        function c=mat2vec(obj)
            % c=~() retrieves data matrix as colum vector ignoring zeros
            % important redefinition of superclass method
            c=nonzeros(obj.data);
        end
        
        function sm2=copy(obj,sm)
            %sm2=~(sm) copy sm by copying data into a new sparse matrix
            obj.requireOnly(isa(sm,'SparseMatrix'),'local','same class');
            sm2=SparseMatrix(sm.data);
            try  % case of non-matching versions
                sm2.units=sm.units;
                sm2.descr=sm.descr;                
            catch
                
            end
        end               
                
               
    end
    
    %% graphics
    
    methods
        
           
        function graph_data(obj,open_new)
            % show spy results in the open window (default)
            if nargin <2 || open_new
                prepfigure(1);
            end
            
            if obj.density <0.3 || numel(obj.data)>1e7
               spy(obj.data);
            else
               imagesc(obj.data); colorbar; 
            end
            
            occupation= nnz(obj.data)/numel(obj.data);
            title({strtrim(obj.descr),...
                ['occupation ',num2tex(100*occupation,...
                '%3.2g','none'), '%']},...
                'fontsize',12);
        end
        
        function graph_distribution(obj, open_new)
            % ~() show distribution density of non-zero elements in open win
            if nargin <2
                open_new=true;
            end
            graph_distribution@HMatrix(obj, open_new);
            
            addstr=[];
            if ~isreal(obj.data)
                addstr='real part of ';
            end
            
            xlabel(strtrim([addstr,'values of non-zeros',' ',obj.units_str]),'fontsize',12);
            title({strtrim(obj.descr),['distribution of non-zeros \mu=',...
                num2str(obj.mean,'%3.2g'),' \pm ', num2str(3*obj.std,'%3.2g'),...
                ' (3\sigma)']},...
                'fontsize',12);
        end
        
        
        
    end
    
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='is sparse matrix';
            ok= issparse(obj.data);
        end
        
    end
    
end



