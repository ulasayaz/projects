function sf= num2tex(h,format,interpreter)
% str=~(h,format,[interpreter]) ... interpreter in ['none','latex'], default latex
assert(isnumeric(h) && isscalar(h),'first parameter must be a scalar number');

if nargin <2 || isempty(format)
    format='%3.2e';
end
if nargin <3
    interpreter='latex';
end
if strcmp(interpreter,'latex')
    dollar='$';    
else
    dollar='';    
end

if h==0
    sf=num2str(h);
else
    sf=  strrep(num2str(h,format),'00','');
    sf=  strrep(sf,'e-0','e-');
    sf=  strrep(sf,'e0','e');
    pos_e=strfind(sf,'e');
    if ~isempty(pos_e)
        if ~strcmp(sf(end),'+')
            sf=[strrep(sf,'e','\cdot 10^{'),'}'];
        else
            sf=sf(1:max(1,end-2));
        end
        sf=[' ',dollar,sf,' ',dollar];
    else
        sf=strrep(sf,'Inf',[' ',dollar,'\infty',dollar,' ']);
    end
end

end