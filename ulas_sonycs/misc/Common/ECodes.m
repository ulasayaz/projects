classdef ECodes
    % ECodes is an enumeration of error and warning codes
    % each code has an integer level, the lower levels have higher priority
    % warnings whose code has a level higher than the corresponding limit level
    % in class NConst will be suppressed.
    %
    % Example
    % ECodes.error_none.description
    % ECodes.error_none.level
    %
    % Note
    % not realisable as singleton, because the constructor cannot be called
    % outside enumeration block and is automatically called even before static methods
    
    properties
        level
        description
    end
    
    enumeration        
        error_none ('no error',9)
        error_fileNotFound ('File not found')
        error_overflow ('overflow error');
        error_unknown_source_type('unknown source file');
        error_not_implemented_yet('not implemented yet',1);
        error_OS('operating system error');
        error_AP('application program error');
        error_matlab('matlab error',4);
        error_constructor('constructor failed',1);

        warning_none ('no warnings',9);
        
         
    end
    
    methods(Access=private)
        function obj= ECodes(str,lv)
            if nargin==1
                lv=4;
            end
            obj.level=lv;
            obj.description=str;
        end
    end
    
    
    
end

