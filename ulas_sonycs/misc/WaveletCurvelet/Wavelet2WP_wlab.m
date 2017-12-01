classdef Wavelet2WP_wlab < Wavelet2D_wlab
    % 2-dim. wavelet packet transformation using the toolbox <<wavelab>>
    %
    % Example: 
    %====================================
    % --- choose signal (image):
    %     signal=Signal2D.make_exp_symm();
    %     signal=Signal2D.make_fromImage('cameraman.bmp');
    %
    % --- discrete wavelet trafo (multi-resolution analysis):
    %     w=Wavelet2WP_wlab();  w.set_basisname('db2'); 
    %     w.show_test;
    %     res=w.test_framefeatures(); display(res);
    % --- set signal: 
    %     w= Wavelet2WP_wlab(signal);
    % --- discrete wavelet packet trafo (multi-resolution analysis):
    %     w.dec;
    %     w.show_resultd;   
    %
    
    properties
        stats  %  stat tree of wavelet packet
        bob    %  Best2dBasis
        entropy_type_name  %<string> method to chose optimal subtree
        entropy_param      %<double> method parameter
    end
    
    %% constructor and commands
    methods
        function obj=Wavelet2WP_wlab(signal)
            % constructor
            global WAVELABPATH
            global WLVERBOSE            
            if nargin==0
                 signal=Signal2D([]);
            end
            obj = obj@Wavelet2D_wlab(signal);
            WLVERBOSE='No';  % 'No' or 'Yes'
            
            obj.entropy_name='l^p';
            obj.entropy_param=1; % l1-norm
            
        end
        
        function set_algo(obj)
            set_algo@MultiResTransform2(obj);
            obj.algo.name='WPT'; % wavelet packet trafo
            obj.algo.toolbox='wavelab';
        end
        
        
        function dec(obj,lvl)
            % ~([lvl]) Wavelet packet decomposition
            obj.requireOnly(~isempty(obj.img),'local','non-empty sample size');
            obj.requireOnly(~isempty(obj.qmf),'local','needs filter');
            if nargin <2
                lvl=obj.deepestlev();
            end                        
            obj.lev=lvl;
            obj.stats=Calc2dStatTree('WP',obj.img-mean(reshape(obj.img,[],1)),...
                obj.lev,obj.qmf,obj.entropy_name,obj.entropy_param);
            obj.bob=Best2dBasis(obj.stats,obj.lev);
            obj.C=FPT2_WP(obj.bob,obj.img,obj.qmf);
            
            obj.ensureOnly(obj.isDWTdone,'local', 'wavelet trafo done');
        end
        
        function yn=rec(obj,wpc)
            %yn=~(wpc) % signal reconstruction from wavelet packet decomp wpc
            yn=IPT2_WP(obj.bob,wpc,obj.qmf);
        end
        
        
    end
    
    
    %% queries
    methods
        
        
        
    end
    
    %% filters
    methods
        
        
    end
    
    %% tests
    methods (Static)
        
        function str=msgDoDWT()
            str='Wavelet packet decomposition available (call dec).';
        end
        
        function w=make()
            % contructor with universal name (callable by parent)
            w=Wavelet2WP_wlab();
        end
        
    end
    
    %% graphics
    methods
        
    end
    
end

