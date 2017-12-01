classdef HartleyTrafo1 < Fourier1
    % 2-dim. Hartley transformation.
    %
    % Example:
    % Example
    % 1.) define a signal, e.g.
    %     signal=TimeSignal.make_sinus(3);
    %
    % 2.) create object
    %     f=HartleyTrafo1(signal);
    %
    % 3.) apply transform
    %     f.dec;
    % 4.) show result
    %     f.graph_trafo;
    %
    % Test frame features:
    % res=f.test_framefeatures(); display(res);
    %
    % GT Stuttgart, 2014
    %
    %% commands and constructor
    methods
        function obj=HartleyTrafo1(signal)
            % constructor
            if nargin==0
                signal=TimeSignal();
            end
            obj= obj@Fourier1(signal);
            obj.repfun=@(F) F;
           
        end
        
        
        function set_algo(obj)
            obj.algo.version=0.9;
            obj.algo.versiondate='14.3.2014';
            obj.algo.name='Hartley1';
            obj.algo.toolbox='matlab';
            obj.algo.mfilename=class(obj);
            obj.algo.author='G. Troll';
        end               
        
        function dec(obj)
            % ~() decomposition whose inverse is rec, i.e. rec(dec)=Id
            % yielding an ONB (check with test_framefeatures() );
            % i.e. operation computeFrame by applying the linear dec operation
            % on the canonical basis (e_i) will yield an ONB (but is no
            % longer involution).
            obj.require(obj.ts.numel>0,'local','non-empty sample size');
            % both dec and rec must either do without or with fftshift
            % to be inverses of each other! Here we chose not to shift:
            obj.C= fft(obj.ts.xn);
            % normalize to ONB, otherwise involution
            obj.C=obj.C/sqrt(numel(obj.C));  
            obj.C=real(obj.C)-imag(obj.C);
            
        end
        
        function yn=rec(obj,cC)
            % ~(cC) signal reconstruction from frame coefficients C
            % should satisfy rec(dec)=Id            
            
            yn= fft(cC);
            % normalize to ONB otherwise it is an isinvolution
            yn=yn/sqrt(numel(yn));  
            yn=real(yn)-imag(yn);            
        end
        
        function yn= analyze(obj,x)
            % yn=~(x) decompose x using an ONB
            yn= fft(x);
            yn=yn/sqrt(numel(yn));  % normalize to ONB
            yn=real(yn)-imag(yn);     
        end
        
        function xn= synthesize(obj,y)
            % yn=~(x) synthesize as inverse operation to analyze
            
            xn= fft(y);
            xn=xn/sqrt(numel(xn));  % normalize to ONB
            xn=real(xn)-imag(xn);     
        end
        
    end
    
    %% queries
    methods
    
        function wvn= basisname(obj)
            % multi-resolution short name
            wvn='Hartley';
        end     
        
        function mat=C2graph(obj)
            % mat=~() convert coefficient vector to graphical output form   
           mat=fftshift(obj.C); 
        end
       
    end
    
     %% tests
    methods (Static)
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=HartleyTrafo1();
        end
        
    end
    
    %% graphics
    
    methods 
        
        function graph_trafo(obj,open_new)
            % show |fft2| without additional data in open figure window
            % used together with obj.dec
            if nargin <2 || open_new
               prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            plot(obj.dec2graph(obj.repfun(obj.C)));          
            tittext={['Transformation: |',obj.basisname,'|'],obj.ts.signalname};
            title(tittext,'fontsize',12);
        end
        
        function graph_distribution(obj, open_new)
            % ~() show distribution of frame coefficients in open window
            obj.require(~isempty(obj.C),'local',obj.msgDoTrafo);
            obj.requireOnly(nargin<2 || islogical(open_new),'local', 'boolean');
            
            if nargin <2 || open_new
                prepfigure(obj.fig,obj.algo,obj.figopt);
            end
            
            img_dwt=obj.C2vec();
            HM=HMatrix(img_dwt);
            HM.graph_distribution(false);
            
            h=get(gca,'Title');
            present_titstr=get(h,'String');
            addtit=['distribution of coefficients (',obj.basisname,')'];
            new_titstr=present_titstr;
            new_titstr{1}=addtit;
            if ~iscell(present_titstr)
                new_titstr{2}=present_titstr;
            end
            title(new_titstr,'fontsize',12);
            
        end
        
    end
    
    
end


