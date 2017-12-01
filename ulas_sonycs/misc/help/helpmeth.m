function helpmeth( className, methodName )
% ~(className, [methodName]) shows short form of all public methods of a class
% or optionally the short form of the method <<methodName>> in the stated class.
% If methodName is missing or a regular expression only public methods will
% be shown.
%
% The short form of a methods contains the function signature, the first comment line, 
% its preconditions and its postconditions. The optional function name funcname 
% may be a regular expression.
%
% Example
% helpmeth NConst
% helpmeth NConst fig.*
%

%
% author : Guenter Troll
% date   : November 2011 - March 2012
%
    
    nconst=NConst.instance;
    if nargin ==1 && isvarname(className)
        nconst.shortForm(className);
    elseif nargin==2 && isvarname(className) && ...
            ischar(methodName) && ~isempty(methodName)        
        nconst.shortForm(className, methodName);
    else
        warning([nconst.idtoolbox,':',mfilename],...
            'please specify valid class name and/or non-empty method name');
    end

end

