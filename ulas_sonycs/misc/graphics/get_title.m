function str = get_title(ah)
% str=~([ah]) get title of axis handle ah (default is current axis handle of
% current figure).
if nargin < 1 || isempty(ah)
    ah=gca;
end
h=get(ah,'Title');
str=get(h,'String'); 

end

