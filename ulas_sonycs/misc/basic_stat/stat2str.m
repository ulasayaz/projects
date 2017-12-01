function strres = stat2str( data,varargin )
% stres=~(data,[opts]) returns statistical quantities of data (e.g. mean and std)
% data ... array or matrix of data
% opts.statfunc1/2 ... 1st/2nd statistical function chosen
% opts.format1/2 ... formats for 2 statistical values

if nargin>1
    if ~isstruct(varargin{1})
        opts=struct;
        nopts=size(varargin,2);
        if nopts>=1, opts.format1=varargin{1}; end
        if nopts>=2, opts.format2=varargin{2}; end
        if nopts>=3, opts.interpreter= varargin{3}; end
    else
        opts=varargin{1};
    end
end

if ~exist('opts','var')
    opts=struct;
end

if ~isfield(opts,'statfunc1') || isempty(opts.statfunc1)
    % 1st statistical function chosen
    opts.statfunc1=@mean;
end
if ~isfield(opts,'statfunc2') || isempty(opts.statfunc2)
    % 2nd statistical function chosen
    opts.statfunc2=@std;
end
if ~isfield(opts,'format1') || isempty(opts.format1)
    opts.format1='%3.1e';
end
if ~isfield(opts,'format2') || isempty(opts.format2)
    opts.format2=opts.format1;
end
if ~isfield(opts,'interpreter') || isempty(opts.interpreter)    
    opts.interpreter='none';
end

stat1=opts.statfunc1(data(:));
stat2=opts.statfunc2(data(:));
strres=num2tex(stat1,opts.format1,opts.interpreter);
if abs(stat2)>0
    strres2=num2tex(stat2,opts.format2,opts.interpreter);
    if isequal(opts.statfunc2,@std)
        strres=[strres,'\pm',strres2];
    else
        strres=[strres,' (\delta=',strres2,')'];
    end
end


