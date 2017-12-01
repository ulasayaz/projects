function vecstr = vec2str( vec, format, delim )
%vec2str vector to string
% transforms a vector of numerical values into a string
% Autor: G. Troll
% SMT 2011

assert(nargin<3 || (ischar(delim) && length(delim)>2),'if given, delim must be string of at least length 3');
if nargin<2 
    format=[];
end

if nargin<3 
    delim='[], ';
end   

vecstr=delim(1);
 
for i=1:length(vec)
    if i>1
        vecstr=[vecstr,delim(3:end)];
    end
    if ~isempty(format)
        shortstr= @(h,format) strrep(strrep(num2str(h,format),'e-00','e-'),'e+00','e+');
        vecstr=[vecstr,shortstr(vec(i),format)];  
    else
        shortstr= @(h,format) strrep(strrep(num2str(h),'e-00','e-'),'e+00','e+');
        vecstr=[vecstr,shortstr(vec(i))];          
    end    
end

vecstr=[vecstr,delim(2)];

end

