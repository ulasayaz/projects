classdef NConst < Singleton1
    % NConst is a singleton containing constants (DP SINGLETON)
    % i.e. no matter where NConst is instantiated or reinstantiated
    % the same object with the same properties whether declared <<Constant>>
    % or not is used everywhere.
    %
    %
    % Example
    % -- accessing constants
    % nconst=NConst.instance;
    % nconst.idproject
    % -- or alternatively:
    % NConst.idproject
    % -- accessing non-constant properties:
    % nconst.contractLevel
    % -- short form of methods
    % nconst.shortForm('DC','.*')
    % -- uml form of classes (or use helpuml script with more features)
    % nconst.umlForm('DC');
    %
    %
    % Note
    % --- prevented assignment: nconst.pi2=1;
    %     (Observe that NConst.pi2=1 would overwrite by creating a structure
    %     NConst with field pi2)
    % --- NConst contains its own assertion methods (require and ensure) because it
    %     cannot be based on DC.
    %
    
    %
    % author  : Guenter Troll
    % date    : November 2011 - March 2012
    %
    
    %% fixed properties
    properties (Constant)
        
        % theses constant cannot be changed
        idproject='Project 1';   % identifier of this project
        requireLevel=1;   % level of preconditions < level of postconditions
        invariantLevel=2; % level of invariants is in between pre- and postconditions
        ensureLevel=3;    % level of postconditions > level of preconditions
        pi2=2*pi;         % 2*pi
        AD1900=693962     % datenum('01-01-1900') used to check if a number can be date
        
        
        %Units
        unit_lenght= 'm';
        unit_time= 's';
        
        %plot types
        plot_linear='plot';
        plot_semilogx='semilogx';
        plot_semilogy='semilogy';
        plot_loglog='loglog';
        plot_types='$plot$semilogx$semilogy$loglog$';
        
        
    end
    
    properties (Constant, Hidden = true)
        strRequire='REQUIRE'       % keyword for preconditions
        strEnsure='ENSURE'         % keyword for postconditions
        strInvariant='INVARIANT'   % keyword for invariants
        strCheck='CHECK'           % keyword for checks
        strViolated='violated at ' % violation string
        umlJava_keywords={'import'}
        umlJava_primitives=...        % primitive datatypes of Java used in UML methods
            {'boolean_','char_','int_','long_','float_','double_'};
        umlIgnoreClasses=... % classes excluded from UML graphs
            {'string','integer','cellarray','any','array','struct','logical','real','handle','function_handle'};
        umlIgnoreOptionalDefault={'DC','CellContainer','CellContainerTS',...
            'KeyClass','HVector','HMatrix'}  % ignored if not isUML_showAll
    end
    
    %% changeable properties
    properties
        warningLevel=4       % warnings of a higher level will be suppressed
        testLevel = 1        % test procedures of a higher level will be ignored
        contractLevel = 3    % 3 means all; contracts with higher level will not be verified
        logFileHandle        % handle for actual log file
        stopOnError = true   % flag to determine behaviour if an error occurs
        sendWarningsToConsole=true; % flag to determine output of warnings
        maxCountFigs_default=-1;    % maximal number of open figure for one builder, infinite if -1
        figLeftCorner=[]       % <pair(double)> position of left corner of figure win
        figSize=[]             % <pair(double)> width and height of figure win
        isUML_nocomments=true  % do not include class comment line in uml diagram
        isUML_short=false      % no attributes and methods shown in UML diagram
        isLayout_horizontal=false % graphics layout, e.g. umlgraph or graphviz
        umlIgnoreOptional       % hidden classes in UML
        isUML_showAll=false;   % also show hidden classes
    end
    
    properties (Access=private)
        logFileNameStandard = 'Project.log' % standard log file name
        logPathStandard = userpath;        % standard path of log file
        warningsOff={'MATLAB:DELETE:Permission',...
            'MATLAB:hg:set_chk'} % switched off warnings
        isUML_large=false        % abbreviate function list in uml graphs
        cellstr={}               % <cellarray(string)>
    end
    
    %% constructor methods
    methods(Access=private, Hidden = true)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.
        function obj = NConst()
            logpath=obj.logPathStandard;
            if isempty(logpath)
                logpath='/Users/ulasayaz/Desktop/temp';
            end
            if ~isstrprop(logpath(end),'alphanum')  
                logpath=logpath(1:end-1);
            end
            obj.set_logFile(obj.logFileNameStandard,logpath);
            for j=1:length(obj.warningsOff)
                warning('off',obj.warningsOff{j});
            end
            obj.umlIgnoreOptional=obj.umlIgnoreOptionalDefault;
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
                uniqueInstance = NConst();
            end
            singleObj = uniqueInstance;
        end
        
        function obj = loadobj(a)
            % ~(a) is called automatically when loading a workspace
            % where a is the saved object which has to be assigned to obj.
            % passing thru instance ensures the identity of the singleton
            % by assigning also the persistant variable to a.
            obj= NConst.instance(a);
        end
        
    end
    
    %% commands
    methods
        
        function set_contractLevel(obj,lv)
            % ~(lv) sets the maximal contractLevel to be verified to lv
            obj.require(lv,'DC.isnatural0(local)','contractLevel is non-negative integer');
            obj.contractLevel=lv;
        end
        function set_stopOnError(obj,clv)
            % ~(clv) if clv then a contract violation will throw an exception
            % otherwise the error will be written into the log file
            obj.stopOnError=clv;
        end
        
        function set_testLevel(obj,lv)
            % ~(lv) sets maximal testLevel to be executed to lv
            obj.require(lv,'DC.isnatural0(local)','testLevel is non-negative integer');
            obj.testLevel=lv;
        end
        
        function set_warningLevel(obj,lv)
            % ~(lv) sets the maximal warningLevel to be reported to lv
            obj.require(lv,'DC.isnatural0(local)','testLevel is non-negative integer');
            obj.warningLevel=lv;
        end
        
        function set_warningState(obj, flag)
            % ~(flag) set warning of toolbox routines on (flag=true) or off
            if flag
                warning('on', obj.idproject);
            else
                warning('off', obj.idproject);
            end
        end
        
        function set_figLeftCorner(obj,lc)
            % ~(lc) set left corner of figure windows to lc
            obj.require(lc,'isempty(local) || (obj.isnaturalArray(local) && length(local)==2)',...
                ' is empty or a pair of positive integers');
            obj.figLeftCorner=reshape(lc,1,[]);
        end
        
        function set_figSize(obj,s)
            % ~(lc) set width and height of figure windows to vector s
            obj.require(s,'isempty(local) || (obj.isnaturalArray(local) && length(local)==2)',...
                ' is empty or a pair of positive integers');
            obj.figSize=reshape(s,1,[]);
        end
        
        function set_figToUpperRight(obj)
            % ~() set figure position and size to upper right corner of screen
            screenSize = get(0,'ScreenSize');
            obj.figSize=floor([screenSize(3)/3, screenSize(4)/2]);
            obj.figLeftCorner=floor([2/3*screenSize(3), screenSize(4)/2]);
        end
        
        function set_logFile(obj,filename,pathname)
            % ~(filename, pathname) creates handle for log file
            obj.check({filename,pathname}, 'ischar(local{1}) && ischar(local{2})',...
                'filename and pathname are strings.');
            obj.logFileHandle=...
                fopen(fullfile(pathname,filename),'a');
            obj.ensure([],'true',' verify invariant');
        end
        
        function fn=getlogFileName(obj)
            % ~() actual log file name used
            fn=fopen(obj.logFileHandle);
        end
        
        function addToUMLIgnoreList(obj, classname)
            % ~(classname) add to cell umlIgnoreClasses
            obj.require(classname,'iscellstr(local) || isvarname(local)', ' is possible classname');
            obj.umlIgnoreOptional=[obj.umlIgnoreOptional, classname];
        end
        
        function set_umlIgnoreClasses(obj, classnames)
            % ~() sets umlIgnoreOptional to classnames, or if missing to default.
            % to delete umlIgnoreOptional, choose empty list of classnames {}.
            % to reset to default, class without argument.
            obj.require(nargin<2 || iscellstr(classnames),'local',...
                'if stated, classnames are a cellstring');
            if nargin <2
                obj.umlIgnoreOptional=obj.umlIgnoreOptionalDefault;
            else
                obj.umlIgnoreOptional=classnames;
            end
        end
        
        function writeToLog(obj,textstring)
            % ~(textstring) write into log file
            obj.require([],'~isempty(obj.logFileHandle)','LogFile not set');
            warning('off', 'MATLAB:printf:BadEscapeSequenceInFormat');
            fprintf(obj.logFileHandle , [datestr(clock),': ', textstring]);
            fprintf(obj.logFileHandle,'\n'); % 2nd command in case the first was badly formatted
        end
        
        
        
    end
    
    %% queries
    methods
        
        function s=get_warningState(obj)
            % ~() get warning state of this toolbox
            s=warning('query',obj.toolbox);
        end
        
        function fileList = getAllFiles(obj,dirName)
            % fileList=~([dirName]) file list of the directory tree under node dirName
            if nargin <1
                dirName=pwd;
            end
            obj.require(dirname,'~isempty(local) && ischar(local)', ' dirname is string');
            dirData = dir(dirName);    %# Get the data for the directory dirName
            dirIndex = [dirData.isdir];  %# Find the index for directories
            fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
            if ~isempty(fileList)
                fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
                    fileList,'UniformOutput',false);
            end
            subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
            validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
            %#   that are not '.' or '..'
            for iDir = find(validIndex)                  %# Loop over valid subdirectories
                nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
                fileList = [fileList; NConst.getAllFiles(nextDir)];  %# Recursively call getAllFiles
            end
            
        end
        
        function fullname= findFileOnPath(obj,filename)
            % fullname=~(filename) search for filename on matlab path
            % if found the full name including path will be returned
            % if not found an empty set will be returned
            obj.require(filename,'~isempty(local) && ischar(local) && isempty(strfind(local,filesep))',...
                ' dirname is string without file separator (path)');
            fullname=[pwd,filesep,filename];
            sdir=dir(fullname);
            found= ~isempty(sdir);
            
            if ~found
                % get all paths on matlab path
                pn=path;
                matlabpath=textscan(pn,'%s','Delimiter',pathsep,'bufsize',2*length(pn));
                matlbapath=matlabpath{1};
                pl=length(matlbapath);
                % locate m-file
                found=false;
                
                j=0;
                while ~found && j<pl
                    j=j+1;
                    fullname=[matlbapath{j},filesep,filename];
                    sdir=dir(fullname);
                    found=~isempty(sdir);
                end
            end
            if ~found
                fullname=which(filename);
                found=~isempty(fullname);
            end
            if ~found
                fullname=[];
            end
            obj.ensure(fullname,'isempty(local)  || exist(local)==2', 'if found file exists');
            
        end
        
        
        function scl=superclasses(obj, classname)
            % ~(classname) retrieves all superclasses of classname
            % removes package names
            obj.require(classname,'isvarname(local)','classname is valid name');
            scl={classname};
            try
                mc = meta.class.fromName(classname);
                superclasses=[];
                if ~isempty(mc)
                    superclasses=mc.SuperclassList;
                end
                for j=1:length(superclasses)
                    sc=superclasses(j).Name;
                    % remove package name
                    lp=strfind(sc,'.');
                    ispackform=~isempty(lp);
                    if ispackform  % reduce package form to standard form
                        sc=sc(lp(end)+1:end);
                    end
                    if ~isempty(sc)
                        scl=[scl; obj.superclasses(sc)] ;
                    end
                end
                scl=unique(scl);
            catch
            end
        end
        
        
        function scl= subclasses(obj,class0)
            % ~(classname) retrieves all subclasses of classname
            % found in any sister folder or sister subfolder of the folder
            % of classname
            obj.require(class0,'isvarname(local)','class0 is valid variable name');
            scl={};
            startpath=fileparts(which(class0));
            if ~isempty(startpath)
                startpath= fullfile(startpath,'..');
                dirs=dir(startpath);
                for j=1: length(dirs)
                    dn=dirs(j).name;
                    if dn(1) ~='.'
                        classnames=obj.lsclassnames(dn);
                        for k=1:length(classnames)
                            if ismember(class0,obj.superclasses(classnames{k}))
                                scl=[scl;classnames{k}];
                            end
                        end
                    end
                end
                scl=unique(scl);
            end
        end
        
        
        function scl= clientclasses(obj,supplierlist)
            % ~(supplierlist) retrieves all client classes of supplierlist
            obj.require(supplierlist,'~isempty(local) && iscellstr(local)',...
                ' supplierlist is non-empty cellstring');
            scl={};
            startpath=fileparts(which(supplierlist{1}));
            if ~isempty(startpath)
                wait_handle = waitbar(0,'Following client classes ...');
                startpath= fullfile(startpath,'..');
                dirs=dir(startpath);
                steps=length(dirs);
                for j=1: steps
                    waitbar(j / steps)
                    dn=dirs(j).name;
                    if dn(1) ~='.'
                        classnames=obj.lsclassnames(dn);
                        for k=1:length(classnames)
                            if obj.isclient(classnames{k},supplierlist)
                                scl=[scl;classnames{k}];
                            end
                        end
                    end
                end
                scl=unique(scl);
                close(wait_handle );
            end
        end
        
        
        
        
        function [status, msg]= umlForm(obj, classname)
            % ~(classname, [funcname]) creates a UML class diagram of <<classname>>
            % by transforming the Matlab class to a Java class, calling
            % UMLGraph to create the UML diagram as a png file in the temporary folder and
            % displaying it. OBSERVE that umlForm requires the installation of umlgraph.
            % Usage Note for developpers:
            % -- A matlab property a within a properties block is augmented by declaration tags, e.g.
            % a @ <TimeStamp>.
            % -- A declaration tag does not need a property
            % preceding it, to allow to include a method parameter type as a
            % supplier, e.g.  @ <Hashtable>
            % -- For a container type the element type can be entered e.g.
            % container  <CellContainer(SingleGraphics)>
            
            % --- usage:
            % nconst=NConst.instance; nconst.umlForm({'DC','NConst','ECodes'});
            % --- alternatively, call the matlab script helpuml which also accepts a package
            % instead of a class name.
            % --- another possibility is to call obj.umlForm for any object of a class
            % inheriting from DC.
            
            % @<umlgraph>  needs umlgraph
            
            obj.require(nargin<2 || (~isempty(classname) && (iscellstr(classname) || isvarname(classname))),'local', ...
                ' classname is valid variable name or cellstring of class names');
            if nargin<2
                classname=class(obj);
            end
            
            javafn='umlinput';
            umlsource=fullfile(tempdir,javafn);
            umlpng=[umlsource,'.png'];
            if exist(umlpng,'file')
                delete(umlpng);
            end
            
            funcname='.*';
            res=struct;
            res.classname={};res.sf={}; res.superclass={};
            res.errormsg=[];
            
            wait_handle = waitbar(0,'Analysing classes ...');
            
            if iscellstr(classname)  % analyse list of classes
                obj.isUML_large=true;
                steps=length(classname);
                for j=1:steps
                    if  steps==1 || (~ismember(classname{j},obj.umlIgnoreClasses) &&...
                            (obj.isUML_showAll || ~ismember(classname{j},obj.umlIgnoreOptional)))
                        res2=obj.analyseClass(strtrim(classname{j}), funcname, false);
                        res.classname=[res.classname,res2.classname];
                        res.sf=[res.sf; res2.sf];
                        waitbar(j / steps)
                    end
                end
            else  % analyse class and direct superclasses, if any
                obj.isUML_large=false;
                res=obj.analyseClass(classname, funcname, false);
                if ~isempty(res.superclass)
                    res.classname={classname};
                    steps=length(res.superclass);
                    for j=1:steps
                        waitbar(j/steps,wait_handle,'Analysing superclasses ...');
                        res2=obj.analyseClass(res.superclass{j}, funcname, false);
                        res.classname=[res.classname,res2.classname];
                        res.sf=[res.sf; res2.sf];
                    end
                end
            end
            close(wait_handle );
            
            % add genera options and hide certain classes
            res.sf=[obj.umloptions; res.sf; obj.java_ignore(setdiff(classname,obj.umlIgnoreOptional))];
            
            fid=fopen(fullfile(tempdir,[javafn,'.java']),'w');
            for j=1:length(res.sf)
                outtxt=strrep(res.sf{j},'%','%%');
                fprintf(fid,[outtxt,'\n'],'%s');
            end
            fclose(fid);
            
            if ~isempty(res.errormsg)
                obj.warningReport('NConst:umlForm',ECodes.error_fileNotFound,res.errormsg);
            else
                
                ec=[];
                [status, msg]=system(['umlgraph ',umlsource,' png']);
                if status==0  && exist(umlpng,'file') % show uml graph
                    os = system_dependent('getos');
                    if strcmp(os, 'Microsoft Windows XP')
                        [status, msg]=system(['mspaint ',umlpng, ' &']);
                    else
                        [status, msg]=system(umlpng);
                    end
                    
                elseif status==0
                    ec=ECodes.error_AP;
                end
                if status ~=0
                    ec=ECodes.error_OS;
                end
                if ~isempty(ec)
                    obj.warningReport('OS:umlgraph',ec,msg);
                end
                
            end
            
            
        end
        
        function shortForm(obj, classname, funcname)
            % ~(classname, [funcname]) shows the short form of all public functions
            % or of function <<funcname>> including signature, pre- and postconditions.
            % example: nc=NConst.instance; nc.shortForm('shortForm','NConst')
            % funcname can be a regular expression like 'get_.*' in which
            % case the output is restricted to public methods.
            % If funcname is missing or = '.*' the short form of all public methods
            % of a class will be shown.
            obj.require(nargin<2 || isvarname(classname),'local', ' classname is valid variable name');
            obj.require(nargin<3 || isempty(funcname) || ischar(funcname),...
                'local',' funcname need not be a variable name but must be at least a string');
            if nargin<2
                classname=class(obj);
            end
            if nargin<3 || isempty(funcname)
                funcname='.*';
            end
            res=obj.analyseClass(classname, funcname,true);
            
            for j=1:length(res.superclass)
                res2=obj.analyseClass(res.superclass{j}, funcname);
                res.sf=[res.sf; {'SUPERCLASS'}; {' '}; res2.sf];
            end
            
            % show result in the help window
            help.class=classname;
            help.func=funcname;
            if isempty(res.errormsg)
                obj.show(res.sf,help);
            else
                obj.show(res.errormsg,help);
            end
        end
        
        function warningReport(obj,id,code,addMsg, stackpos)
            % ~(id, code, addMsg) shows/logs warning, code belongs to class ECodes
            % id --> messageID analogue to Matlab warning('MSGID','MESSAGE')
            % example: id='Mechatronics:classname';
            % suppresses warnings if code.level> warningLevel
            % otherwise, writes into the log file and also to console
            % if sendWarningsToConsole
            obj.require(id,'ischar(local)','id is string');
            obj.require(code,'isa(local, ''ECodes'')','class of code is ECodes');
            obj.require(nargin<4 || isempty(addMsg) || ischar(addMsg),'local','if stated addMsg is string');
            if nargin<4 || isempty(addMsg)
                addMsg='';
            end
            if nargin <5
                stackpos=2;
            end
            if code.level<=obj.warningLevel
                stack=dbstack();
                stacklevel=stack(min(stackpos, length(stack)));
                reportString=([stacklevel.name, ' >> line number ', ...
                    num2str(stacklevel.line) , ': ',code.description]);
                if obj.sendWarningsToConsole
                    disp([reportString,10,addMsg]);
                end
                obj.writeToLog(['Warning: ', reportString,'\n',addMsg]);
            end
        end
        
        
        
        
    end
    
    
    %% local assertions
    % public and visible (so that they are available for non-object
    % scripts!
    methods
        
        %COPIED ensure, require and check " from DC
        %to avoid infinite recursion when inheriting from DC
        
        function ok = require(obj,local, condition, name)
            % local routine to verify precondition
            if obj.contractLevel >= obj.requireLevel
                conditionException=[];
                try
                    ok = eval(condition);
                catch conditionException
                    ok=false;
                end
                %eval invariants
                if ok && obj.contractLevel >= obj.invariantLevel && ismethod( obj, 'invariant')
                    [ok, contract.descr]= obj.invariant;
                    contract.type=obj.strInvariant;
                    contract.loc=class(obj);
                else
                    contract.type=obj.strRequire;
                    contract.descr=name;
                    contract.loc=[];
                end
                if ok
                else % ok might be empty!
                    stack=dbstack(-1);
                    contract.line=stack(2).line;
                    reportString=([contract.type,' ',obj.strViolated,contract.loc, ' >> line number ', num2str(contract.line) ,...
                        ': ' ,contract.descr]);
                    if obj.stopOnError
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
            % local routine to verify postcondition
            if obj.contractLevel >= obj.ensureLevel
                conditionException=[];
                try
                    ok = eval(condition);
                catch conditionException
                    ok=false;
                end
                if ok && ismethod( obj, 'invariant')
                    [ok, contract.descr]= obj.invariant;
                    contract.type=obj.strInvariant;
                    contract.loc=class(obj);
                else
                    contract.type=obj.strRequire;
                    contract.descr=name;
                    contract.loc=[];
                end
                if ok
                else % might be empty
                    stack=dbstack(-1);
                    contract.line=stack(2).line;
                    reportString=([contract.type,' ',obj.strViolated,contract.loc, ' >> line number ', num2str(contract.line) ,...
                        ': ' ,contract.descr]);
                    if obj.stopOnError
                        err = MException([class(obj),':E',contract.type],reportString);
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
        
        function ok = check(obj,local, condition, name)
            % local routine to check local condition
            if obj.contractLevel >= obj.ensureLevel
                conditionException=[];
                try
                    ok = eval(condition);
                catch conditionException
                    ok=false;
                end
                if ok
                else % ok might be empty
                    stack=dbstack(-1);
                    contract.type=obj.strCheck;
                    contract.line=stack(2).line;
                    contract.loc=[];
                    contract.descr=name;
                    reportString=([contract.type,' ',obj.strViolated,contract.loc, ' >> line number ', num2str(contract.line) ,...
                        ': ' ,contract.descr]);
                    if obj.stopOnError
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
        
    end
    
    
    methods (Hidden)
        % order relation  the same as with handle class but hidden
        function TF = eq(H1,H2)
            TF= eq@handle(H1,H2);
        end
        function TF = ne(H1,H2)
            TF= ne@handle(H1,H2);
        end
        function TF = lt(H1,H2)
            TF= lt@handle(H1,H2);
        end
        function TF = le(H1,H2)
            TF= le@handle(H1,H2);
        end
        function TF = gt(H1,H2)
            TF= gt@handle(H1,H2);
        end
        function TF = ge(H1,H2)
            TF= ge@handle(H1,H2);
        end
        
    end
    
    methods (Static, Hidden)
        
        function cn=lsclassnames(sourcename)
            % ~(sourcename) retrieves all classnames in folder sourcename
            c=what(sourcename);
            cn={};
            % extract possible class names
            % if what has found several folders
            % take the first one which contains classes
            for j=1:length(c)
                [~,cn,~]= cellfun(@fileparts,c(j).m,'UniformOutput',false);
                % extract actual class names
                isclass= @(x) exist(x,'class');
                cn=cn(cellfun(isclass,cn)>0);
                if ~isempty(cn)
                    break;
                end
            end
        end
        
        
    end
    
    methods (Access=private)
        
        function show(obj, celltxt, help)
            % ~(celltxt) show celltxt in help window
            % cf. http://undocumentedmatlab.com/blog/customizing-help-popup-contents/
            obj.require(celltxt,'iscell(local)', 'celltxt is cell');
            obj.require(help,'isstruct(local)',' help is struct');
            
            delete(fullfile(tempdir,'___*___.m'));
            
            if isfield(help, 'class') && ~isempty(help.class)
                help.title=['___',help.class,'___'];
            else
                help.title='SHORTFORM';
            end
            help.fn=fullfile(tempdir,[help.title,'.m']);
            fid=fopen(help.fn,'w');
            fprintf(fid,[help.title,'\n'],'%s');
            if isfield(help, 'class') && isfield(help,'func')
                fprintf(fid,['%%SHORTFORM: ',help.class,'.(',help.func,')\n'],'%s');
                fprintf(fid,('%%\n'),'%s');
            end
            for j=1:length(celltxt)
                if ~isempty(celltxt{j})
                    outtxt=strrep(celltxt{j},'%','%%');
                    fprintf(fid,['%%',outtxt,'\n'],'%s');
                end
            end
            fclose(fid);
            
            vers=version;
            try
                ok=str2double(vers(1:4))  >= 7.13;
            catch
                ok=false;
            end
            if ok
                % for Matlab Version from 7.13.0.564 (R2011b)
                helpUtils.errorDocCallback(help.fn);
            else
                jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
                jTextArea = jDesktop.getMainFrame.getFocusOwner;
                jClassName = 'com.mathworks.mlwidgets.help.HelpPopup';
                jPosition = java.awt.Rectangle(0,0,700,500);
                helpTopic = help.fn;
                javaMethodEDT('showHelp',jClassName,jTextArea,[],jPosition,helpTopic);
            end
        end
        
        function [res,decl]=analyseClass(obj, classname, funcname, isMatlab)
            % res=~(classname,[funcname] returns cell of strings
            % representing the short forms of functions
            % res.superclass   ... (string) name of superclass of class classname
            % res.classname    ... (string) name of class
            % res.sf           ... (cellstr) short form of class classname
            % res.errormsg     ... (string) error message
            
            obj.require(classname,'isvarname(local)',...
                ' classname is valid variable name');
            obj.require(nargin<3 || isempty(funcname) || ischar(funcname),...
                'local',' funcname need not be a variable name but must be at least a string');
            obj.require(nargin<4 || islogical(isMatlab),'local',' is logical');
            
            if nargin<3 || isempty(funcname)
                funcname='.*';
            end
            if nargin<4 || isempty(isMatlab)
                isMatlab=true;
            end
            
            res=struct;
            [funcNames,methodfilterOK]=obj.identifyPossibleFuncNames(funcname,classname);
            
            % search local path
            mfilename=[classname,'.m'];
            fullfilename=obj.findFileOnPath(mfilename);
            found= ~isempty(fullfilename); % && exist(classname,'class');
            res.superclass={};
            res.classname=classname;
            res.sf=[];
            res.errormsg=[];
            res.isAbstract=false;
            
            head=[];
            classcomment='';
            
            if ~found
                res.errormsg={'class m-file not found (name possibly misspelled)'};
                classtxt= ['classdef ', classname];
                if isMatlab
                    res.sf{1}=classtxt;
                else
                    decl=struct;
                    decl.vars={}; decl.types={}; decl.elem_types={};
                    res.sf{1}= obj.java_classhead(res, '',decl);
                end
            else
                res.sf=cell(length(funcNames)+1,1);
                sdir=dir(fullfilename);
                fid=fopen(fullfilename,'r');
                classtxt = textscan(fid,'%s','Delimiter','\n','whitespace','','bufsize',max(sdir.bytes,4095));
                classtxt=classtxt{1};
                fclose(fid);
                
                % locate classdef and first comment
                keyexpr='^\s*classdef\s+'; % keywords preceded only by whitespace
                fstart= @(x) regexp(x,keyexpr, 'once');
                idxtable=cell2mat(cellfun(fstart,classtxt,'UniformOutput',false));
                startline=find(idxtable,1);
                endline=startline;
                commentline=startline;
                if ~isempty(startline)
                    if length(classtxt)>endline && ~isempty(regexp(classtxt{startline},'\.{3}\s*$','once'))
                        % line continuation marked by <<...>>
                        endline=endline+1;
                    end
                    if length(classtxt)>endline && ~isempty(regexp(classtxt{endline+1},'^\s*%','once'))
                        % next line is a comment
                        commentline=endline+1;
                        classcomment= classtxt(commentline);
                        classcomment=strtrim(strrep(classcomment{1},'%',''));
                    end
                    head= classtxt(startline:max(endline,commentline));
                else
                    head={['classdef ', classname]};
                    startline=1; endline=1;
                end
                
                try
                    % identify classname as written in m-file and its superclasses
                    if ~isempty(head)
                        firstline=strcat(head{1:endline-startline+1});
                        % classname as written in m-file:
                        res.classname=obj.identify_classname(firstline);
                        res.superclass=obj.identify_superclasses(firstline);
                    end
                catch errNConst
                    if obj.warningLevel<=4
                        warning([obj.idproject,':NConst'],'identification of superclasses failed');
                        warning(errNConst.message);
                    end
                end
                
                interface_locator= @(x) regexp(x,'methods\s*(\s*Abstract','once');
                idxtable=cell2mat(cellfun(interface_locator,classtxt,'UniformOutput',false));
                res.isAbstract= any(idxtable);
                
                if ~isempty(head)
                    if isMatlab
                        res.sf{1}=head;
                    else
                        % extract attributes (properties) for java output
                        
                        try % matlab properties might fail if class has syntax errors
                            possattr=properties(classname);
                            attrfilterOK=true;
                        catch properties_exception
                            attrfilterOK=false;
                            possattr=[];
                            obj.warningReport('MATLAB:properties',ECodes.error_matlab,' calling <<properties>>');
                        end
                        
                        % prep step: find tagged lines:
                        % declaration tags are identified by %@<...>
                        decl_mark='%\s*@\s*<.*>';
                        decl_fun= @(x) regexp(x,decl_mark, 'once');
                        % idxtable contains the declaration tags, if any
                        [~,~,~,idxtable]=cellfun(decl_fun,classtxt,'UniformOutput',false);
                        locs= ~cellfun(@isempty,idxtable);
                        lines= classtxt(locs);
                        
                        % step 1: the variable name(s) (separated by commas) precede(s)
                        % the declaration tag in the lines containing such a tag:
                        var_fun= @(x) regexp(x, '^\s*(\w+)|(\w+\,?\w+)(?=.*%)','once');
                        [~,~,~,vars]=cellfun(var_fun,lines,'UniformOutput',false);
                        decl=struct;
                        decl.vars=strtrim(vars);
                        
                        % step 2: visibility of variables is determined by the
                        % result of the properties function: possattr
                        % appearance in this list means attribute is
                        % visible
                        if attrfilterOK
                            mfun= @(x) ismember(x,possattr);
                            decl.visible=cellfun(mfun, decl.vars);
                        else
                            decl.visible=ones(size(decl.vars));
                        end
                        
                        % step 3: try to determine type of variable
                        types=idxtable(locs);
                        % the leading alphanumeric string is the type:
                        type_fun= @(x) regexp(x, '\w+','once');
                        [~,~,~,decl.types]=cellfun(type_fun,types,'UniformOutput',false);
                        % the element type, if any, is contained in paranthesis within
                        % the declaration tag,  e.g.
                        % <CellContainer(SingleGraphics)> has type
                        % CellContainer and element type SingleGraphic
                        element_fun= @(x) regexp(x, '(?<=\()\w+(?=\))','once');
                        [~,~,~,decl.elem_types]=cellfun(element_fun,types,'UniformOutput',false);
                        res.sf{1}= obj.java_classhead(res, classcomment,decl);
                    end
                end
            end
            
            % Is classname entered case-sensitively identical with
            % the classname in the m-file ?
            isclassnameOk=isempty(res.classname) || strcmp(res.classname,classname);
            isclassidentified= found && isclassnameOk;
            
            if ~isclassnameOk
                res.errormsg={['classname in m-file is : ', res.classname];...
                    ['classname entered differs (case-sensitive): ',classname]};
                res.superclass={};
            end
            
            if ~isclassidentified
                classtxt={ ['classdef ', classname]};
                res.classname=classname;
            end
            
            
            % find all function identifiers in classtxt
            
            % require backwards: key word function possibly followed by a word and an equal
            alter='((\[.*\])|(\w+))';  % include case of several output variables [.,.]
            prefix=['^\s*function\s+((',alter,'\s*=\s*)|(\s*))'];
            root='\w+';
            % exclude ahead: alphanumerics followed possibly by blanks and definitively by an equal sign
            postfix='(?!(\w*\s*=))';
            func_fun= @(x) regexp(x, ['(?<=',prefix,')',root,postfix],'once');
            [~,~,~,foundFunc]=cellfun(func_fun,classtxt,'UniformOutput',false);
            idx=~cellfun(@isempty,foundFunc);
            % the functions start at lines startlinesFunc in classtxt
            startlinesFunc=find(idx);
            % sort function names
            [foundFunc,sidx]=sort(foundFunc(idx));
            startlinesFunc=startlinesFunc(sidx);
            
            % look for key words indicating the end of a function
            keyexpr=['(^\s*properties\s*)|',...   % keywords preceded only by whitespace
                '(^\s*methods\s*)|',...
                '(^\s*function\s+)'];
            fend=@(x) ~isempty(regexp(x,keyexpr, 'once'));
            idx=cell2mat(cellfun(fend,classtxt,'UniformOutput',false));
            endlinesCand=find(idx);  % lines where a function can end
            
            % alternative to get methods
            % shows all methods instead only the public and visible ones
            if ~methodfilterOK
                funcNames=foundFunc;
            else
                % select those functions that appear both in foundFunc and
                % funcNames, i.e. that are public and visible
                [funcNames, idxff] = intersect(foundFunc, funcNames);
                startlinesFunc=startlinesFunc(idxff);
            end
            
            
            
            for j=1:length(startlinesFunc)
                if isMatlab
                    % locate essentials of all functions in funcNames
                    % i.e. first comment line, pre- and postcondition for normal
                    % functions and ok and descr for invariants.
                    startline= startlinesFunc(j);
                    fn=funcNames{j};
                    if ~isempty(startline)
                        % find the end line of the present function fn
                        endlineidx=find(endlinesCand>startline,1);
                        endline=endlinesCand(endlineidx);
                        if isempty(endline)
                            endline=length(classtxt);
                        end
                        
                        % look for key words for signature, precondition, postcondition
                        if strcmp(fn,'invariant')
                            % keywords preceded only by white space
                            keyexpr='(^\s*(descr|ok)\s*=)|(invariant@)';
                        else
                            % keywords preceded by white space, alphabetic, numeric, or underscore character
                            % and one dot to detect the OO-call
                            % <<obj.require>>.
                            % Extended to also detect require without dot
                            % notation.
                            keyexpr='^((\s*\w*\.)|(\s*))(require|ensure)';
                        end
                        
                        % limit search to function block:
                        functxt=classtxt(startline:endline-1);
                        fpick=@(x) ~isempty(regexp(x,keyexpr, 'once'));
                        % find lines idxtable containig keyexpr
                        idxtable=cell2mat(cellfun(fpick,functxt,'UniformOutput',false));
                        % first line contains function signature
                        idxtable(1)=true;
                        
                        %
                        if length(idxtable)>1
                            % add second line if it is a comment line or has
                            % already been retained:
                            idxtable(2)=idxtable(2) || ~isempty(regexp(functxt{2},'^\s*\%.*', 'once'));
                            
                            % continuation lines are lines ending with <<...>>
                            fpick=@(x) ~isempty(regexp(x,'.*\.{3}\s*$', 'once'));
                            % all continuation line in block:
                            idxtable2=cell2mat(cellfun(fpick,functxt,'UniformOutput',false));
                            % add line following continuation line of a keyword line:
                            idxtable3=circshift(idxtable2 & idxtable,1);
                            % add line following the next iterated contination line if
                            % any
                            idxtable3=idxtable3 | circshift(idxtable2 & idxtable3,1);
                        else
                            idxtable3=false;
                        end
                        res.sf{j+1}=functxt(idxtable| idxtable3);
                    end
                    
                else
                    % java output: just function name
                    res.sf{j+1}= obj.java_function(res,funcNames,j);
                end
            end
            idx_close=min(length(res.sf),length(funcNames)+1);
            if ~isMatlab
                % end marker of java class
                res.sf{idx_close}=[res.sf{idx_close}; {'  }'}];
            end
            
            % collapse subcells of sf
            if ~all(cellfun(@isempty,res.sf))
                offs=0;
                % preallocate cell <<coll>>
                Lcells=cellfun(@length,res.sf);
                addlines=Lcells>0;
                LL= sum(Lcells+addlines);
                coll=cell(LL,1);
                
                for j=1:length(res.sf)
                    kk=0;
                    if iscell(res.sf{j})
                        for k=1:length(res.sf{j})
                            if ~isempty(res.sf{j}{k})
                                kk=kk+1;
                                coll{offs+kk,1}=res.sf{j}{k};
                            end
                        end
                    else
                        kk=1;
                        coll{offs+kk,1}=res.sf{j};
                    end
                    if ~isempty(res.sf{j})
                        offs=offs+kk+1;
                        coll{offs,1}=' ';
                    end
                end
                coll(cellfun(@isempty,coll)) = [];
                res.sf=coll;
            end
            
            obj.ensure(isstruct(res),'local','return value is structure');
            
        end % function analyseClass
        
        
        function c=umloptions(obj)
            % c= ~() options for umlgraph
            if obj.isUML_short
                alltag='*';
            else
                alltag='* @opt all';   % set umlgraph options
            end
            if obj.isLayout_horizontal
                hortag= '* @opt horizontal';
            else
                hortag='*';
            end
            c={'/**';alltag;hortag;'* @hidden';'*/';'class UMLOptions {}'};
            obj.ensure(iscellstr(c),'local','returns cellstr');
        end
        
        function c=java_ignore(obj, exceptions)
            % c=~() list of classes to ignore in UML diagrams if not in exceptions
            if nargin <2
                exceptions={};
            end
            if iscell(exceptions)
                ignore=setdiff(union(obj.umlIgnoreClasses,obj.umlIgnoreOptional),exceptions);
            else
                ignore=obj.umlIgnoreClasses;
            end
            ignore=[ignore,obj.umlJava_primitives];
            li=length(ignore);
            c=cell(2*li,1);
            for j=1:li
                c{2*j-1}='/** @hidden */';
                c{2*j}=['class ',ignore{j}, ' {}'];
            end
        end
        
        function c=java_classhead(obj,res,classcomment, decl)
            % c=~(res,classcomment,decl) java form of class head
            % e.g.   class BankAccount extends Asset {}
            %        /** @extends InsurableItem */
            obj.require(length(decl.vars)==length(decl.types),'local','same lengths');
            obj.require(length(decl.types)==length(decl.elem_types),'local','same lengths');
            lsc= length(res.superclass);
            lv=length(decl.vars);
            le=length(decl.elem_types)-sum(cellfun(@isempty,decl.elem_types));
            lh=max(1,lsc)+1;
            obj.cellstr=cell(lh+lv+le+5,1);
            p=1;
            % build obj.cellstr using obj.add
            
            % first the javadoc tags
            if obj.isUML_nocomments
                p=obj.add(p,'/** ');
            else
                lccom=length(classcomment);
                if lccom >20
                    p=obj.add(p,['/** @note ',classcomment(1:floor(lccom/2)),'-']);
                    p=obj.add(p,classcomment(floor(lccom/2)+1:end));
                elseif lccom > 0
                    p=obj.add(p,['/** @note ',classcomment]);
                else
                    p=obj.add(p,'/** ');
                end
            end
            
            p=obj.add(p,'* @opt nodefillcolor LemonChiffon');
            
            for j=2:lsc
                p=obj.add(p,['* @extends ',res.superclass{j}]);
            end
            
            if lv==0
                txt=' */';
            else
                txt='*';
            end
            p=obj.add(p,txt);
            
            
            for j=1:lv
                % client-provider relation
                if ~isempty(decl.types{j}) && (j==1 || ~isequal(decl.types{j-1},decl.types{j}))
                    txt=['* @navassoc - - - ',decl.types{j}];
                else
                    txt='* ';
                end
                % composition relation
                if ~isempty(decl.elem_types{j})  && (j==1 || ~isequal(decl.elem_types{j-1},decl.elem_types{j}))
                    p=obj.add(p,['* @has 1 -  1..* ',decl.elem_types{j}]);
                end
                if j==lv
                    txt= [txt, ' */'];
                end
                p=obj.add(p,txt);
            end
            
            % second the class declaration
            if res.isAbstract
                txt=['abstract class ', res.classname];
            else
                txt=['class ', res.classname];
            end
            if lsc>0
                txt=[txt, ' extends ', res.superclass{1}, '{'];
            else
                txt=[txt, ' {'];
            end
            p=obj.add(p,txt);
            
            % class attributes
            for j=1:lv
                visible='public ';
                if ~decl.visible(j)
                    visible='protected ';
                end
                if ~isempty(decl.vars{j}) && ~isempty(decl.types{j})
                    p=obj.add(p,[visible,decl.types{j},' ',decl.vars{j},';']);
                else % allow for the possibilty to declare a type without attribute
                    p=obj.add(p,' ');
                end
            end
            
            c=obj.cellstr;
            
            % open '{' will be close later
            
            obj.ensure([],'~isempty(strfind(obj.cellstr{1},''/**''))','comment opened in 1st line');
        end
        
        function p=add(obj,p,txt)
            % p=~(p,txt) add txt to the cellstr obj.cellstr
            obj.cellstr{p}=txt;
            p=p+1;
        end
        
        function c=java_function(obj,res,funcnames, j)
            % c=~(res,funcnames,j) add j-th function name to java class declaraton
            obj.require(j>=1 && j<=length(funcnames),'local',' j within bounds');
            obj.require(iscellstr(funcnames),'local','is cellstr');
            c='';
            funcname=funcnames{j};
            if ismember(funcname,obj.umlJava_keywords)
                funcname=strcat(funcname,'_');
            end
            if obj.isUML_large
                % show only 1 set-function, same for insert and remove
                [ok,funcname]=obj.funcname_collapse(j,funcnames,'set');
                if ~ok
                    [ok,funcname]=obj.funcname_collapse(j,funcnames,'insert');
                end
                if ~ok
                    [~,funcname]=obj.funcname_collapse(j,funcnames,'remove');
                end
                
                idx=strfind(funcname,'set');
                if ~isempty(idx) && idx(1)==1
                    idx2=[];
                    if j>1
                        idx2=strfind(funcnames{j-1},'set');
                    end
                    if ~isempty(idx2) && idx2(1)==1
                        funcname=[];
                    else
                        funcname='set_';
                    end
                end
            end
            
            if ~isempty(funcname) && ~strcmp(res.classname, funcname)
                % As constructor will be added automatically by umlgraph
                % it should not appear in function list
                c={['public void ', funcname, ' () {}']};
            end
        end
        
        function [ok,funcname]=funcname_collapse(obj,j,funcnames,head)
            % reduce certain function names with many variants
            funcname=funcnames{j};
            idx=strfind(funcname,head);
            ok=~isempty(idx) && idx(1)==1;
            if ok
                idx2=[];
                if j>1
                    idx2=strfind(funcnames{j-1},head);
                end
                if ~isempty(idx2) && idx2(1)==1
                    funcname=[];
                else
                    funcname=[head,'_'];
                end
            end
        end
        
        function cn=identify_classname(obj,str)
            % cn=~(str) identify class name cn from class header str
            obj.require(str,'~isempty(local) && ischar(local)',' is non-empty string');
            class_tag= '(?<=classdef\s*(\(.*\))?\s*)\w+';
            [~,~,~,cn]=regexp(str,class_tag,'once');
            obj.ensure(cn,'ischar(local)',' is (possibly empty) string');
        end
        
        function s= identify_superclasses(obj,str)
            % s= ~(str) identify superclasses fromm class header str
            obj.require(str,'~isempty(local) && ischar(local)',' is non-empty string');
            startpos= strfind(str,'<');
            if ~isempty(startpos)
                startpos=startpos(1);
                sc= str(min(length(str),startpos+1):end);
                ampersand_pos=strfind(sc,'&');
                ampersand_pos=[ampersand_pos,length(sc)+1];
                countclasses=length(ampersand_pos);
                s=cell(1,countclasses);
                startpos=1;
                for j=1:countclasses
                    s{j}=strtrim(sc(startpos:ampersand_pos(j)-1));
                    startpos=min(length(sc),ampersand_pos(j)+1);
                    lp=strfind(s{j},'.');
                    ispackform=~isempty(lp);
                    if ispackform  % reduce package form to standard form
                        s{j}=s{j}(lp(end)+1:end);
                    end
                end
            else
                s={};
            end
            obj.ensure(s,'iscellstr(local)',' is (possibly empty) cell of strings');
        end
        
        function [funcNames,methodfilterOK]=identifyPossibleFuncNames(obj,funcname, classname)
            % ~(funcname) restrict search space to visible names if possible
            % uses matlab function <<methods>>#
            % methods fails, if m-file has syntx errors
            obj.require(funcname,'~isempty(local) && ischar(local)',' is non-empty string');
            methodfilterOK=logical(exist(classname,'class'));
            
            if methodfilterOK && ~isvarname(funcname)
                % if funcname is no valid variable name, assume that
                % it is a regular expression.
                % filter VISIBLE methods by this regular expression
                try % matlab methods might fail if class has syntax errors
                    possmethods=methods(classname);
                catch methods_exception
                    methodfilterOK=false;
                    possmethods=[];
                    obj.warningReport('MATLAB:methods',ECodes.error_matlab,' calling <<methods>>');
                end
                possnames=[];
                if ~isempty(possmethods)
                    % filter possible methods with the regular expression
                    % funcname
                    f=@(x) ~isempty(regexpi(x,funcname,'once'));
                    possnames=find(cell2mat(cellfun(f,possmethods,'UniformOutput',false)));
                end
                if ~isempty(possnames)
                    funcNames=possmethods(possnames);
                else
                    funcNames={funcname};
                end
            else
                funcNames={funcname};
            end
            obj.ensure(methodfilterOK,'islogical(local)', 'is boolean');
        end
        
        function b= isclient(obj, checkclass, supplierlist)
            % ~(class, suppplierlist) has <<class>> one of the listed supppliers?
            obj.require(checkclass,'isvarname(local)',...
                ' client is valid variable name');
            obj.require(supplierlist,'~isempty(local) && iscellstr(local)',...
                ' supplierlist is non-empty cellstring');
            b=NaN;
            % search local path
            mfilename=[checkclass,'.m'];
            fullfilename=obj.findFileOnPath(mfilename);
            found= ~isempty(fullfilename);
            if found
                fid=fopen(fullfilename,'r');
                sdir=dir(fullfilename);
                classtxt = textscan(fid,'%s','Delimiter','\n','whitespace','','bufsize',max(sdir.bytes,4095));
                classtxt=classtxt{1};
                fclose(fid);
                suppliertag=supplierlist{1};
                for j=2:length(supplierlist)
                    suppliertag=[suppliertag,'|',supplierlist{j}];
                end
                % look for start tag "%+@+<" and end tag ")" or "(" or ">"
                supp_expr=['%\s*@\s*<.*(',suppliertag,')\s*(\)|\(|>)'];
                supp_fun= @(x) ~isempty(regexp(x,supp_expr, 'once'));
                s=cellfun(supp_fun,classtxt);
                b=max(s)>0;
            end
        end
        
    end
    
    %% invariant
    methods
        function [ok,descr]=invariant(obj)
            % class invariant
            descr='log file handle exists';
            ok=obj.logFileHandle>0;
        end
    end
end