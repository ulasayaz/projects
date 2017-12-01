function ok = requireOnly(local, condition, name)
% ~(local, c, name) verifies precondition <<condition>> (logical) w/o invariant.
% This script can be used instead of class DesignByContract
% to verify constructor preconditions that precede a call
% of a constructor of a superclass, because such calls are only allowed
% in Matlab (as of version 2011b) if there are no prior references to
% the object to be constructed.
% Details cf. method requireOnly in class DesignByContract
% In contrast to that method the present function has only access
% to class properties via the argument cell <<local>>.

%
% author  : Guenter Troll
% date    : November 2011 - March 2012
%
    

global nconst            % use the global nconst=NConst.instance if defined

try                      % test if nconst is defined
    docheck=nconst.contractLevel >= nconst.requireLevel;
catch exceptionNConst    % if nconst not defined verify condition 
    docheck=true;
end

if docheck
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
        contract.type='REQUIRE';
        contract.loc=[];
        contract.descr=name;
        
        stack=dbstack();
        stacklevel=stack(min(length(stack),stackpos));
        if isempty(contract.loc)
            contract.loc=strrep(stacklevel.name,'.',':');
        end
        contract.line=stacklevel.line;
        
        reportString=([contract.type,' ','violated at ',contract.loc, ' >> line number ', num2str(contract.line) ,...
            ': ' ,contract.descr]);
        if ~isempty(exceptionNConst) || nconst.stopOnError
            err = MException([contract.loc,':',contract.type],reportString);
            if ~isempty(conditionException)
                throwAsCaller(addCause(err,conditionException));
            else
                throwAsCaller(err);
            end
        else
            nconst.writeToLog(reportString);
        end
    end
    
end


end

