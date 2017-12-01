classdef SessionConst < Singleton1
    % SessionConst is a singleton class containing session constants
    % i.e. constants which are specific to a particular session o
    %
    % Example
    % sc= SessionConst.instance;    
    %
    
    properties (Constant)
        
    end
    
    properties (SetAccess=private)
        warningsOff={}           % switched off warnings
        hash                     %@ <Hashtable> for session objects
    end
    
    methods(Access=private, Hidden = true)
        
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.
        function obj = SessionConst()
            % private constructor
            % one needs the initialisations in the constructor
            % instead of in the property list in order to avoid
            % problems when reloading the workspace.
            % These problems appear when copying objects after reloading
            % (tested under Matlab 2011b).
            obj.hash = Hashtable(250);
        end
    end
    
    methods(Static)
        
        function singleObj = instance(c)
            % ~([c]) get unique instance of class (normal call without parameter)
            % without parameter: call private constructor
            % with parameter of this class: assign persistant variable to
            % that parameter (use this option when loading a workspace)
            persistent uniqueInstance
            if nargin>0
                uniqueInstance=c;
            end
            if isempty(uniqueInstance) || ~isvalid(uniqueInstance)
                uniqueInstance = SessionConst();
            end
            singleObj = uniqueInstance;
        end
        
        function obj = loadobj(a)
            % ~(a) is called automatically when loading a workspace
            % where a is the saved object which has to be assigned to obj.
            % passing thru instance ensures the identity of the singleton
            % by assigning also the persistant variable to a.
            obj= SessionConst.instance(a);
        end
        
        
    end
    
   
end