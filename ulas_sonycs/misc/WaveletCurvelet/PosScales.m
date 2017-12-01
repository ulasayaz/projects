classdef PosScales
    % position scales for signal axes
    
    properties        
        symbol         %@<char>
        label          %@<char>
        description    %@<char>
    end
    
    enumeration
       n({'n','m','k'},'index n', 'integer index')
       x({'x','y','z'},'spatial scale','length scale')               
       t('t','time t', 'time scale')    
       quantile('q','quantile','quantile of distribution')
    end
        
    methods(Access=private)
        function obj= PosScales(sym,lbl, txt)
            obj.symbol=sym;
            obj.label=lbl;
            obj.description=txt;
        end
    end
    
end


