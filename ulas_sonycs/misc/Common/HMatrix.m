classdef HMatrix < DC
    % matrix implemented as handle
    % Example
    % hm= HMatrix(randn(2,3));
    % hm.item(2,3)
    % hm.put(2,3,1.1);
    % hm.item(2,3);
    % hm.show;
    %
    
    properties (SetAccess=?DC)
        data                    %@ (matrix) data matrix
        units=[]                %@ (string) physical units of matrix entries
        descr                   %@ (string) description
    end
    
    
    % constructor & commands
    methods
        
        function obj=HMatrix(hdata,u)
            % ~([hdata],[u]) constructs object from matrix hdata with optional unit u.
            obj.requireOnly(nargin<1 || isnumeric(hdata),'local','if specified values must be numeric');
            obj.requireOnly(nargin<2 || ischar(u),'local','if specified units must be a string');
            if nargin==0
                % needed to allow empty calls from subclasses
                % and prevent their premature testing of invariant
                hdata=[];
            end
            
            %set units
            if nargin>=2
                obj.units=u;
            end
            obj.data=hdata;
            obj.descr='matrix data';
            if nargin>0
                % prevent premature testing of invariant by subclasses
                obj.ensure([],'true','verify invariant');
            end
        end
        
        
        function put(obj,i,j,s)
            % ~(i,j,s) set data at index i,j to s
            obj.require(i,'obj.isnatural(local) && local<=obj.size(1)','valid index');
            obj.requireOnly(j,'obj.isnatural(local) && local<=obj.size(2)','valid index');
            obj.requireOnly(s,'obj.isnumber(local)','scalar number');
            obj.data(i,j)=s;
            obj.ensure({i,j,s},'local{3}==obj.data(local{1},local{2})', 'data set');
        end
        
        function set_data(obj,hdata)
            % ~(hdata) set data matrix to hdata
            obj.requireOnly(isnumeric(hdata),'local','hdata is numeric');
            obj.data=hdata;
            obj.ensure(hdata,'all(local==obj.data)', ' data has been set');
        end
        
        function set_units(obj, u)
            % ~(u) sets units to u
            obj.require(u,'ischar(local)',' is string');
            obj.units=u;
            obj.ensure(u,'strcmp(local,obj.units)', 'units set');
        end
        
        
    end
    
    %% queries
    methods
        
        function c=item(obj,i,j)
            % c=~(i,j) retrieve data element at index i,j
            obj.require(i,'obj.isnatural(local) && local<=obj.size(1)','valid index');
            obj.requireOnly(j,'obj.isnatural(local) && local<=obj.size(2)','valid index');
            c=obj.data(i,j);
        end
        
        function cget(obj)
            % c=~() retrieves data matrix
            c=obj.data;
        end
        
        function c=mat2vec(obj)
            % c=~() retrieves data as column vector
            % can be redefined in subclasses, e.g. SparseMatrix
            c=obj.data(:);
        end
        
        function sm2=copy(obj,sm)
            %sm2=~(sm) copy sm by copying data into a new matrix
            obj.requireOnly(isa(sm,'SparseMatrix'),'local','same class');
            sm2=HMatrix(sm.data);
            try  % case of non-matching versions
                sm2.units=sm.units;
                sm2.descr=sm.descr;
            catch
                
            end
        end
        
        function sm2=clone(obj)
            % sm2=~() clone object
            sm2=obj.copy(obj);
        end
        
        function d = size(obj, dim)
            % d=~() yields size of data matrix
            if nargin<2
                d=size(obj.data);
            else
                d=size(obj.data,dim);
            end
        end
        
        function v=nnz(obj)
            % =~() number of non-zeros
            v=nnz(obj.data);
        end
        
        function v=nonzeros(obj)
            % v=~() returns all non-zero elements
            v=nonzeros(obj.data);
        end
        
        function d = length(obj)
            % d=~() yields length of data matrix
            d=length(obj.data);
        end
        
        function ok =isvector(obj)
            % ok= ~() tests if data matrix is a vector
            ok=isvector(obj.data);
        end
        
        function ok=isnumeric(obj)
            % ok=~() tests property <<data>>
            ok=isnumeric(obj.data);
        end
        
        function ok=isempty(obj)
            % ok=~() tests property <<data>>
            ok=isempty(obj.data);
        end
        
        function ok=issparse(obj)
            % ok=~() tests property <<data>>
            ok=issparse(obj.data);
        end
        
        function v=density(obj)
            % density of sparse matrix
            v=nnz(obj.data)/numel(obj.data);
        end
        
        function mm=max(obj)
            % mm=~() row max of matrix
            mm=max(obj.mat2vec); % use redefinable mat2vec
        end
        
        function mm=min(obj)
            % mm=~() row min of matrix
            mm=min(obj.mat2vec);  % use redefinable mat2vec
        end
        
        function mm=mean(obj)
            % mm=~() row min of matrix
            mm=mean(obj.mat2vec);  % use redefinable mat2vec
        end
        
        function mm=std(obj)
            mm=std(obj.mat2vec);  % use redefinable mat2vec
        end
        
    end
    
    
    %% conversion queries
    
    methods
        
        function mm=make(obj,data)
            mm=HMatrix(data);
        end
        
        function mm=real(obj)
            % mm=~() real part of data
            mm=obj.make(real(obj.data));
        end
        
        function mm=imag(obj)
            % mm=~() imaginary part of data
            mm=HMatrix(imag(obj.data));
        end
        
        function mm=abs(obj)
            % mm=~() absolute value of data
            mm=HMatrix(abs(obj.data));
        end
        
        function mm=phase(obj)
            % mm=~() phase of data (units rad)
            mm=AngleMatrix(atan2(imag(obj.data),real(obj.data)));
        end
        
    end
    
    %% graphics
    methods
        
        function s=units_str(obj)
            if isempty(obj.units)
                s=[];
            else
                s=['[',obj.units,']'];
            end
        end
        
        function prepfigure(obj)
            % made redefinable
            fig=1;
            opt.pixsizeX=1000;
            prepfigure(fig,[],opt);
        end
        
        function graph_data(obj,open_new, maxsize)
            % show spy results in the open window (default)
            if nargin <2
                open_new=true;
            end
            if nargin <3
                maxsize=1e3;
            end
            if open_new
                obj.prepfigure();
            end
            occupation= nnz(obj.data)/numel(obj.data);
            occ_str=[];
            if occupation<0.9
                occ_str=['density ',num2tex(100*occupation,...
                    '%3.2g','none'), '%'];
            end
            if obj.density <0.01
                spy(obj.data);
            else
                colormap('default');
                % resize to save memory
                imgsize=numel(obj.data);
                resize_factor=min(1,maxsize/sqrt(imgsize));
                if resize_factor<1
                    % choose method 'nearest' to preserve qualitative
                    % differences of value distributions             
                    imagesc([1,size(obj.data,2)],[1,size(obj.data,1)],...
                        imresize(real(obj.data),resize_factor,'nearest'));
                    xlabel(['resized by ',num2str(resize_factor,'%3.1g')]);
                else
                    imagesc(real(obj.data));
                end
                colorbar;
                cbfreeze;
            end
            freezeColors;
            
            title({strtrim(obj.descr),occ_str},...
                'fontsize',12);
        end
        
        
        function graph_distribution(obj, open_new)
            % ~() show distribution density in open win
            if nargin <2 || open_new
                obj.prepfigure;
            end
            
            v=obj.mat2vec;
            isdatareal=isreal(v);
            if ~isdatareal
                v=real(v);
            end
            
            N=numel(obj.data);
            sample_idx=randi(N,[min(10,N),1]);
            count_values=length(unique([obj.max; obj.min; obj.data(sample_idx)]));
            
            if  count_values<=3 ... % ksdensity gives false impression in this case
                    || numel(v)>1e6   % ksdensity too slow
                hist(v);
                ylabel('frequency');
            else
                try
                    spread=max(v)-min(v);
                    ksdensity(v,'width',spread/10);
                    ylabel('density','fontsize',12);
                catch
                    hist(v);
                    ylabel('frequency');
                end
            end
            addstr=[];
            if ~isdatareal
                addstr='real part of ';
            end
            xlabel(strtrim([addstr,'data values ',obj.units_str]),'fontsize',12);
            title({['data distribution'],...  % first line should be replacable by clients
                ['\mu=',num2str(obj.mean,'%3.2g'),' \pm ', num2str(3*obj.std,'%3.2g'),...
                ' (3\sigma)']},...
                'fontsize',12);
        end
        
        function show(obj)
            % ~() contents and stats of matrix
            obj.prepfigure;
            p1=subplot(1,2,1);
            obj.graph_data(false);
            sp1Axs = gca;
            sp1AxsRatio = get(sp1Axs,'PlotBoxAspectRatio');
            p2=subplot(1,2,2);
            obj.graph_distribution(false);
            sp2Axs = gca;
            sp2AxsRatio = get(sp2Axs,'PlotBoxAspectRatio');
            
            pos1 = get(p1, 'Position');
            pos2 = get(p2, 'Position');
            %
            % Set width of second axes equal to first
            pos2(2) = pos1(2);
            set(p2,'Position',pos2)
            % set same aspect ratio as for first subplot
            if ~isequal(sp1AxsRatio,sp2AxsRatio)
                set(sp2Axs,'PlotBoxAspectRatio', sp1AxsRatio);
            end
            
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='is numeric matrix';
            ok= isnumeric(obj.data);
        end
        
    end
    
end

