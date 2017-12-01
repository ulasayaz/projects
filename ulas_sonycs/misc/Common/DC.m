classdef DC < handle
    % DesignByContract emulates assertions & handles errors and warnings
    % all classes using design by contract are subclasses of this class.
    %
    % Construction:
    % dc=DesignByContract();
    %
    % switch contract verification on/off
    % -----------------------------------------
    % switch off: dc.set_level(0)
    % switch on preconditions only: dc.set_level(dc.nconst.requireLevel)
    % switch on invariants too: dc.set_level(dc.nconst.invariantLevel)
    % switch on all contracts: dc.set_level(dc.nconst.ensureLevel)
    %
    % consequence of a contract violation
    % -----------------------------------------
    % will throw an exception if nconst.stopOnError or otherwise write the
    % violation into the log file (cf. class NConst).
    % use <<dbstop if error>> to stop execution in block throwing the
    % exception.
    %
    % emulated contract types
    %------------------------------------------
    % Precondition & Invariant  : require(local, condition, name)
    % Postcondition & Invariant : ensure(local, condition, name)
    % ~ without invariant       : requireOnly(...), ensureOnly(...)
    % General(without Invariant): check(local, condition, name)
    % class invariant           : [ok,descr]=invariant(obj)
    %
    % where
    %          local     .... (cell) list of local variables
    %          condition ...  (string) logical proposition to be verified
    %          name      ...  (string) message/name of assertion
    %
    %
    % Example
    % testDC;
    %
    %
    % Sources: 111007_GTroll_ContractDesign_Principles&Matlab_Emulation.ppt
    %          111007_GTroll_KontraktDesign_Principien&Matlab_Emulation.ppt
    %          120305_GTroll_EmulationOfDesignByContractInMatlab.ppt
    %
    
    %
    % author : Guenter Troll
    % date   : November 2011 - March 2012
    %
    
    properties
        nconst                     %@<NConst> single instance of constants
    end
    
    properties (SetAccess=protected, Hidden)
        recursion_flag=false       %@<logical> used to detect recursion e.g. dcopy, update
    end
    
    
    
    %% commands
    methods
        
        function obj=DC()
            % ~() constructor
            % one needs the initialisations in the constructor
            % instead of in the property list in order to avoid
            % problems when reloading the workspace.
            % These problems appear when copying objects after reloading
            % (tested under Matlab 2011b).
            obj.nconst	=   NConst.instance;
            obj.recursion_flag=false;
        end
        
        function c=copyConstructor(obj)
            % c=~() is the generic, i.e. empty constructor
            % must be redefined if a non-empty constructor call is needed
            c=eval(class(obj));
        end
        
        function set_level(obj,lv)
            % ~(lv) set test of contracts (contract level) to level lv
            obj.require(lv,'obj.isnatural0(local)', 'lv is non-negative integer');
            obj.nconst.set_contractLevel(lv)
            obj.ensure(lv,'obj.get==local','contract level set');
        end
        
        function  errorReport(obj,stackpos,additionalMessage)
            % ~([stackpos],[additionalMessage]) error handling analogous to contract design methods
            % ~(stackpos),
            % ~([],additionaMessage)
            % throws an exception if nconst.stopOnError, otherwise write
            % error to log file.
            if nargin>1
                obj.require(stackpos,'isempty(local) || (isnumeric(local) && local>0)',...
                    'if entered, stackpos must be positive');
                if nargin>2
                    obj.require(additionalMessage,'ischar(local)','is string');
                end
            end
            
            if nargin==1 || isempty(stackpos)
                stackpos=2;
            end
            if nargin<3
                additionalMessage='';
            else
                additionalMessage=[': ', additionalMessage];
            end
            
            stack=dbstack();
            stacklevel=stack(min(stackpos, length(stack)));
            
            reportString=([stacklevel.name, ' >> line number ', num2str(stacklevel.line) ,...
                additionalMessage]);
            if obj.nconst.stopOnError
                %  error(reportString);
                err = MException([class(obj),':error'],reportString);
                throw(err);
            else
                obj.writeToLog(['ERROR: ', reportString]);
            end
        end
        
        function warningReport(obj,id,code,addMsg)
            % ~(id, code, addMsg) shows/logs warning, code belongs to class ECodes
            % id --> messageID analogue to Matlab warning('MSGID','MESSAGE')
            % example: id='Mechatronics:classname'; code=ECodes.error_fileNotFound;
            % suppresses warnings if code.level> warningLevel
            % otherwise, writes into the log file and also to console
            % if nconst.sendWarningsToConsole
            if nargin<4
                addMsg=[];
            end
            obj.nconst.warningReport(id,code,addMsg,3);
        end
        
        function writeToLog(obj,str)
            % ~(str) writes string str into log file
            obj.nconst.writeToLog(str);
        end
        
    end
    
    methods (Hidden)
        
        function clone(obj, cc)
            % ~(cc) generic clone method to copy properties of cc onto own properties
            % present class needs get access to all attributes of cc
            mc=metaclass(obj);
            ls=mc.PropertyList;
            fun= @(j) ls(j).Name;
            ls=arrayfun(fun,1:length(ls),'UniformOutput',false);
            for j=1:length(ls)
                obj.(ls{j})=cc.(ls{j});
            end
        end
        
        function deep_clone(obj,c)
            % c=~()  recursive generic deep clone of object
            % present class needs get access to all attributes of cc
            % Cloned are only attributes which also have a deep_clone
            % method. All others are only assigned.
            mc=metaclass(obj);
            ls=mc.PropertyList;
            fun= @(j) ls(j).Name;
            ls=arrayfun(fun,1:length(ls),'UniformOutput',false);
            for j=1:length(ls)
                % call dcopy of each attribute to clone
                try
                    % attribute c.(ls{j}) must have method dcopy
                    obj.(ls{j})=c.(ls{j}).dcopy; % recursive call
                catch err_deepclone
                    % try clause fails if
                    % a) c.(ls{j}) does not have such a method (simple type
                    %    or external class)
                    % b) copyConstructor within dcopy call failed
                    % c) there is no read access to c.(ls{j}).
                    %
                    % case a) will occur without calling dcopy
                    % case c) will be forwarded and caught in dcopy
                    % case b) would be forwarded to this catch clause but is
                    % actually not thrown but treated as a warning within dcopy.
                    %
                    % Thus only cases a) and b) have to be treated here:
                    % use as alternative copy by assignment without throwing an exception.
                    % This is definitely the right thing to do for simple
                    % types like double and arguably the only way to
                    % proceed for external classes.
                    obj.(ls{j})=c.(ls{j});
                end
            end
        end
    end
    
    %% Assertions
    methods (Access = protected, Hidden)
        
        function ok = require(obj,local, condition, name)
            % ~(local, c, name) verifies precondition c (logical)
            % do nothing, if contractLevel <requireLevel
            % not ok implies throwing an exception (if stopOnError)
            % or reporting error in log file
            % local ... cell of local variables
            % condition ... string evaluating to a boolean
            % name ... string describing condition
            %
            if obj.nconst.contractLevel >= obj.nconst.requireLevel
                % test require condition:
                conditionException=[];
                try
                    ok = eval(condition);  % might be empty!
                catch conditionException
                    ok=false;
                end
                stackpos=2;
                if ok
                    if obj.nconst.contractLevel >= obj.nconst.invariantLevel
                        %test invariant:
                        try % try is faster than test ismethod !
                            [ok, contract.descr]= obj.invariant;
                            if ~ok
                                stackpos=stackpos+1;
                            end
                            contract.type=obj.nconst.strInvariant;
                            contract.loc=class(obj);
                        catch conditionException
                            % ok=true if no method <<invariant>>
                            ok=isequal(conditionException.identifier,'MATLAB:noSuchMethodOrField');
                            contract.type=obj.nconst.strInvariant;
                            contract.loc=class(obj);
                            contract.descr=conditionException.message;
                        end
                    end
                else
                    ok=false;   % ok might be empty which counts as false
                    contract.type=obj.nconst.strRequire;
                    contract.loc=[];
                    contract.descr=name;
                end
                if ~ok
                    stack=dbstack();
                    stacklevel=stack(min(length(stack),stackpos));
                    if isempty(contract.loc)
                        contract.loc=stacklevel.name;
                    end
                    contract.line=stacklevel.line;
                    
                    reportString=([contract.type,' ',obj.nconst.strViolated,contract.loc, ' >> line number ', num2str(contract.line) ,...
                        ': ' ,contract.descr]);
                    if obj.nconst.stopOnError
                        err = MException([class(obj),':',contract.type],reportString);
                        % use <<throwAsCaller>> instead of <<throw>> to throw
                        % exception from the calling function instead
                        % of from here.
                        if ~isempty(conditionException)
                            throwAsCaller(addCause(err,conditionException));
                        else
                            throwAsCaller(err);
                        end
                    else
                        obj.writeToLog(reportString);
                    end
                end
                
            end
        end
        
        function ok = ensure(obj,local, condition, name)
            % ~(local, c, name) verifies postcondition c  (logical)
            % do nothing, if contractLevel < ensureLevel
            % not ok implies throwing an exception (if stopOnError)
            % or reporting error in log file
            % local ... cell of local variables
            % condition ... string evaluating to a boolean
            % name ... string describing condition
            %
            if obj.nconst.contractLevel >= obj.nconst.ensureLevel
                % test ensure condition:
                conditionException=[];
                try
                    ok = eval(condition);
                catch conditionException
                    ok=false;
                end
                stackpos=2;
                if ok
                    %test invariant:
                    try % try is faster than ismethod !
                        [ok, contract.descr]= obj.invariant;
                        contract.type=obj.nconst.strInvariant;
                        contract.loc=class(obj);
                        if ~ok
                            stackpos=stackpos+1;
                        end
                    catch conditionException
                        % ok=true if no method <<invariant>>
                        ok=isequal(conditionException.identifier,'MATLAB:noSuchMethodOrField');
                        contract.type=obj.nconst.strInvariant;
                        contract.loc=class(obj);
                        contract.descr=conditionException.message;
                    end
                else
                    ok=false;   % ok might be empty which counts as false
                    contract.type=obj.nconst.strEnsure;
                    contract.descr=name;
                    contract.loc=[];
                end
                
                if ~ok
                    stack=dbstack();
                    stacklevel=stack(min(length(stack),stackpos));
                    if isempty(contract.loc)
                        contract.loc=stacklevel.name;
                    end
                    contract.line=stacklevel.line;
                    
                    reportString=([contract.type,' ',obj.nconst.strViolated,contract.loc, ' >> line number ', num2str(contract.line) ,...
                        ': ' ,contract.descr]);
                    if obj.nconst.stopOnError
                        err = MException([class(obj),':',contract.type],reportString);
                        % use <<throwAsCaller>> instead of <<throw>> to throw
                        % exception from the calling function instead
                        % of from here:
                        if ~isempty(conditionException)
                            throwAsCaller(addCause(err,conditionException));
                        else
                            throwAsCaller(err);
                        end
                    else
                        obj.writeToLog(reportString);
                    end
                end
                
            end
        end
        
        function ok = requireOnly(obj,local, condition, name)
            % ~(local, c, name) verifies precondition c (logical) w/o invariant.
            % Code copied from require (without invariant checking)
            % in order to make caller top of stack when
            % throwing an exception through throwAsCaller.
            % Details cf. method require.
            if obj.nconst.contractLevel >= obj.nconst.requireLevel
                % test require condition:
                conditionException=[];
                try
                    ok = eval(condition);
                catch conditionException
                    ok=false;
                end
                if ok
                else   % ok might be empty
                    stackpos=2;
                    contract.type=obj.nconst.strRequire;
                    contract.loc=[];
                    contract.descr=name;
                    
                    stack=dbstack();
                    stacklevel=stack(min(length(stack),stackpos));
                    if isempty(contract.loc)
                        contract.loc=stacklevel.name;
                    end
                    contract.line=stacklevel.line;
                    
                    reportString=([contract.type,' ',obj.nconst.strViolated,contract.loc, ' >> line number ', num2str(contract.line) ,...
                        ': ' ,contract.descr]);
                    if obj.nconst.stopOnError
                        err = MException([class(obj),':',contract.type],reportString);
                        if ~isempty(conditionException)
                            throwAsCaller(addCause(err,conditionException));
                        else
                            throwAsCaller(err);
                        end
                    else
                        obj.writeToLog(reportString);
                    end
                end
                
            end
        end
        
        function ok = ensureOnly(obj,local, condition, name)
            % ~(local, c, name) verifies postcondition c (logical) w/o invariant
            % Code copied from ensure (without invariant checking)
            % in order to make caller top of stack when
            % throwing an exception through throwAsCaller.
            % Details cf. method require.
            if obj.nconst.contractLevel >= obj.nconst.ensureLevel
                % test ensure condition:
                conditionException=[];
                try
                    ok = eval(condition);
                catch conditionException
                    ok=false;
                end
                stackpos=2;
                if ok
                    %no test of invariant:
                else
                    ok=false;   % ok might be empty which counts as false
                    contract.type=obj.nconst.strEnsure;
                    contract.descr=name;
                    contract.loc=[];
                end
                
                if ~ok
                    stack=dbstack();
                    stacklevel=stack(min(length(stack),stackpos));
                    if isempty(contract.loc)
                        contract.loc=stacklevel.name;
                    end
                    contract.line=stacklevel.line;
                    
                    reportString=([contract.type,' ',obj.nconst.strViolated,contract.loc, ' >> line number ', num2str(contract.line) ,...
                        ': ' ,contract.descr]);
                    if obj.nconst.stopOnError
                        err = MException([class(obj),':',contract.type],reportString);
                        if ~isempty(conditionException)
                            throwAsCaller(addCause(err,conditionException));
                        else
                            throwAsCaller(err);
                        end
                    else
                        obj.writeToLog(reportString);
                    end
                end
                
            end
        end
        
        
        function ok = check(obj,local, condition, name, ctype, stackpos)
            % ~(local, c, name) verifies local condition c  (logical)
            % do nothing, if contractLevel < ensureLevel
            % not ok implies throwing an exception (if stopOnError)
            % or reporting error in log file
            % local ... cell of local variables
            % condition ... string evaluating to a boolean
            % name ... string describing condition
            %
            if obj.nconst.contractLevel >= obj.nconst.ensureLevel
                conditionException=[];
                try
                    ok = eval(condition);
                catch conditionException
                    ok=false;
                end
                if nargin>4 && ~isempty(ctype)
                    contract.type=ctype;
                else
                    contract.type=obj.nconst.strCheck;
                end
                if nargin <6 || isempty(stackpos)
                    stackpos=2;
                end
                if ok
                else % ok might be empty!
                    stack=dbstack();
                    stacklevel=stack(min(stackpos, length(stack)));
                    contract.loc=stacklevel.name;
                    contract.line=stacklevel.line;
                    
                    reportString=([contract.type,' ',obj.nconst.strViolated,contract.loc, ' >> line number ', num2str(contract.line) ,...
                        ': ' ,name]);
                    if obj.nconst.stopOnError
                        
                        err = MException([class(obj),':',contract.type],reportString);
                        if ~isempty(conditionException)
                            throwAsCaller(addCause(err,conditionException));
                        else
                            throwAsCaller(err);
                        end
                    else
                        obj.writeToLog(reportString);
                    end
                end
            end
        end
        
        function a = condEval(obj,local, command)
            % ~(local, command) conditional evaluation of command if contractLevel > ensureLevel
            % local ... cell of local variables
            % command ... string evaluating a Matlab command
            %
            a=[];
            if obj.nconst.contractLevel >= obj.nconst.ensureLevel
                try
                    a = eval(command);
                catch commandException
                    a=NaN;
                end
            end
        end
        
    end
    
    
    %% visible queries
    methods
        function lv=get_level(obj)
            % ~() get current contract level
            lv= obj.nconst.contractLevel;
        end
        
        
        function shortForm(obj, funcName)
            % ~(funcName) shows the short form of a function
            % calling shortForm of class NConst.
            obj.require(nargin<2 || obj.isstring(funcName),'local', ...
                ' funcName is string if stated');
            if nargin<2
                funcName=[];
            end
            obj.nconst.shortForm(class(obj),funcName);
        end
        
        function umlForm(obj)
            % ~() shows the UML diagram of this class
            % calling umlForm of class NConst.
            obj.nconst.isUML_nocomments=false;
            obj.nconst.umlForm(class(obj));
        end
        
        function res=compare(obj,c)
            % ~(c) compare all properties
            obj.require(c,'isa(local,class(obj))','c must be based on object');
            mc=metaclass(obj);
            ls=mc.PropertyList;
            fun= @(j) ls(j).Name;
            ls=arrayfun(fun,1:length(ls),'UniformOutput',false);
            res=struct;
            for j=1:length(ls)
                % by precondition c has all the properties of obj
                res.(ls{j})=isequal(c.(ls{j}),obj.(ls{j}));
            end
        end
        
        function c= copy(obj)
            % c=~()  generic shallow copy of object
            % when a class redefines this method, it can follow this
            % pattern:
            %
            %    c=<<call constructor with no matter what arguments>>;
            %    c.clone(obj);
            %    obj.ensure(c,'obj.iscopy(local)','is copy');
            %
            try
                % needs a copy constructor, which can be redefined in
                % subclasses
                c=obj.copyConstructor;
            catch
                try
                    % try empty constructor explicitly (needed for MATLAB core
                    % classes)
                    c=eval(class(obj));
                catch err_copy
                    errCause = MException('DC:copy', ...
                        'provide empty constructor or redefine generic copy');
                    err_copy = addCause(err_copy, errCause);
                    rethrow(err_copy);
                end
            end
            try
                % needs get access to all attributes of obj
                c.clone(obj);
            catch err_clone
                errCause = MException('Mechatronics_Toolbox:TimeStamp', ...
                    'make all private/protected features accessible to DC or redefine generic clone');
                err_clone = addCause(err_clone, errCause);
                rethrow(err_clone);
            end
            obj.ensure(c,'obj.iscopy(local)','is copy');
        end
        
        function c=dcopy(obj)
            % c=~()  generic deep copy of object (takes longer than shallow copy)
            % when a class redefines this method, it can follow this
            % pattern (cf. class Edge1)
            %
            %    c=<<call constructor with no matter what arguments>>;
            %    c.deep_clone(obj);
            %    obj.ensure(c,'obj.iscopy(local)','is copy');
            %
            if obj.recursion_flag  % avoid infinite loops
                % Use direct assignment, which will be overwritten when returning to
                % the first access of this object. Therefore
                % do not reset the recursion flag here!
                c=obj;
                return;
            end
            obj.recursion_flag=true;
            try
                % needs a copy constructor, which can be redefined in
                % subclasses
                c=obj.copyConstructor;
            catch
                try
                    % try empty constructor explicitly (needed for MATLAB core
                    % classes)
                    c=eval(class(obj));
                catch err_copy
                    state='DC:dcopy';
                    msg=['copy constructor of class ',class(obj),' failed'];
                    obj.warningReport(state,ECodes.error_constructor,msg);
                end
            end
            try
                % needs get access to all attributes of obj
                c.deep_clone(obj);  % call deep clone recursively
            catch err_clone
                state='DC:deep_clone';
                errCause = MException(state, ...
                    'make all private/protected features accessible to DC or redefine generic clone');
                err_clone = addCause(err_clone, errCause);
                rethrow(err_clone);
            end
            
            c.recursion_flag=false;
            obj.recursion_flag=false;
            obj.ensure(c,'~obj.recursion_flag && obj.iscopy(local)','is copy');
        end
        
    end
    
    %% hidden or static queries
    methods (Static, Hidden)
        function ok= isnatural(x)
            % ~(x) true if x is a positive integer
            ok= isnumeric(x) && isscalar(x) && x>0 && floor(x)==x;
        end
        function ok= isnatural0(x)
            % ~(x) true if x is a non-negative integer
            ok= isnumeric(x) && isscalar(x) && x>=0 && floor(x)==x;
        end
        function ok= isnaturalArray(x)
            % ~(x) true if x is an array of positive integers
            ok= isnumeric(x) && all(x>0 & floor(x)==x);
        end
        function ok= ispositive(x)
            % ~(x) true if x is a non-negative integer
            ok= isnumeric(x) && isscalar(x) && x>0;
        end
        function ok= isnonnegative(x)
            % ~(x) true if x is a non-negative integer
            ok= isnumeric(x) && isscalar(x) && x>=0;
        end
        function ok= isnumber(x)
            % ~(x) true if x is a scalar number
            ok= isnumeric(x) && isscalar(x);
        end
        function ok=isdate(x)
            % ~(x) true if x is a numeric date
            ok=isnumeric(x) && isscalar(x) && x>= NConst.AD1900;
        end
        function ok=isstring(x)
            % ~(x) true if x is a nonempty string
            ok=~isempty(x) && ischar(x);
        end
        function ok=isvector3D(x)
            % ~(x) true if x is an row vector in R^3
            ok=isnumeric(x) && all(size(x)==[3,1]);
        end
        function ok=ismatrix3D(x)
            % ~(x) true if x is an matrix in R^3xR^3
            ok=isnumeric(x) && all(size(x)==[3,3]);
        end
        function ok=issingle(x,classname)
            % ~(x) is a single and not an array of objects of class <<classname>>
            ok= isa(x,classname) && isscalar(x)==1;
        end
        function ok= isvectorpair(a,b)
            %  ~(a,b) is a pair of comaptible numeric vectors
            ok=isnumeric(a) && isnumeric(b) && min(size(a))<=1 &&...
                all(size(a)==size(b));
        end
        
        function ok= implies(a, b)
            ok= ~a || b;
        end
        function ok= equiv(a, b)
            ok= (a && b) || (~a && ~b);
        end
        
        function ok = isCellOfStrings(c)
            ok = iscellstr(c(~cellfun(@(x) isempty(x),c)));
        end
        
        function value= iif(varargin)
            % value= iif(condition_1,value_1,...,true,value_final)
            assert(~rem(length(varargin),2),'even number of inputs');
            foo  = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();            
            value=foo(varargin{:});            
        end
        
        function value=iifempty(val1,val2)
            if ~isempty(val1)
                value=val1;
            else
                value=val2;
            end
        end
        
        
    end
    
    methods (Hidden)
        
        function ok=iscopy(obj,c)
            % ~(c) is c a copy by reference of obj?
            % only applicable for classes based on the handle class
            % where eq has not been redefined and tests handle equality.
            % Cannot be static in contrast to the other isXY methods
            ok=isequal(obj,c) && ~eq(obj,c);
        end
        
    end
    
    
    
end