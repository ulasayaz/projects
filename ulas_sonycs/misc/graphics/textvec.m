function textvec(x,y,varargin)
% Beschrifte alle Punkte xv, yv mit den Werten von x or y oder mit opt.data

% *** call structure 1: 
% up to 4 additional arguments in the order of opt fields below.
% 
% *** call structure 2:
% 1 additional argument:
% opt (struct)
%     .valign  ... (string) default 'bottom'
%     .halign  ... (string) default 'center'
%     .outputY ... (bool, default false) output is x (false) or y (true)
%     .data    ... (array or matrix) optional data to be used for text output
%     .format  ... (string) default []
%     .fontsize ... default 10
%
%
% Autor: G. Troll
% SMT 2011

if nargin<3
    opt=struct;
else
    nopt=size(varargin,2);
    if nopt==1 && isstruct(varargin{1})
        opt=varargin{1};  % new call
    else  % downward compatibility to old calls
        if nopt>=1, opt.valign=varargin{1}; end
        if nopt>=2, opt.halign=varargin{2}; end
        if nopt>=3, opt.outputY=varargin{3}; end
        if nopt>=4, opt.format=varargin{4}; end  
    end
end

if ~isfield(opt,'valign') || isempty(opt.valign)
    opt.valign='bottom';
end
if ~isfield(opt,'halign') || isempty(opt.halign)
    opt.halign='center';
end
if ~isfield(opt,'outputY') || isempty(opt.outputY)
    opt.outputY=false;
end
if ~isfield(opt,'format') % allow empty format
    opt.format='%3.2e';
end
if ~isfield(opt,'fontsize') || isempty(opt.fontsize)
    opt.fontsize=10;
end
if ~isfield(opt,'color') || isempty(opt.color)
    opt.color='k';
end

if isfield(opt,'data') && ~isempty(opt.data)
    data=opt.data;
elseif ~opt.outputY
    data=x;
else
    data=y;
end

if ~isempty(opt.format)
    df=@(i)num2str(data(i),opt.format);
else
    df=@(i)num2str(data(i));
end
for i=1:numel(x)
    text(x(i),y(i),df(i),'VerticalAlignment',opt.valign,'HorizontalAlignment',opt.halign,...
        'fontsize',opt.fontsize,'color',opt.color);
end

end



