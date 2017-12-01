classdef FrameSparse < Frame
    % frame using a sparse matrix
    
    properties
    end
    
    methods
        
        function obj=FrameSparse(hdata,u)
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
            obj.descr='sparse frame';
            
            if nargin>0
                % prevent premature testing of invariant by subclasses
                obj.ensure([],'true','verify invariant');
            end
        end
        
        function c=mat2vec(obj)
            % c=~() retrieves data matrix as colum vector ignoring zeros
            % important redefinition of superclass method
            c=nonzeros(obj.data);
        end
        
        function fm2=copy(obj,fm)
            %sm2=~(sm) copy sm by copying data into a new sparse matrix
            obj.requireOnly(isa(sm,'SparseFrame'),'local','same class');
            fm2=FrameSparse(fm.data);
            try  % case of non-matching versions
                fm2.units=fm.units;
                fm2.descr=fm.descr;
            catch
                
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
    
    methods (Static)
        
        function F=test_operators(M)
            % test operations with Parseval frame:
            if nargin <1 || isempty(M)
                M=sqrt(2/3)*[1, -1/2,  -1/2; 0,  sqrt(3)/2,  -sqrt(3)/2];
            end
            if ~isa(M,'Frame');
                F=FrameSparse(M);
            else
                F=M;
            end
            test_operators@Frame(F);
        end
        
    end
    
end

