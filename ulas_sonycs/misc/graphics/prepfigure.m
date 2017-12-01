function fighandle= prepfigure(fig, algo, opt)
% figure output in Fenster der durch algo bestimmten ID
% Die Zahl der pro ID offenen Fenster kann beschraenkt werden:
% Autor: G. Troll
% SMT 2008-2011
%
% INPUT
% fig .... (integer, default 1) max. Zahl offener Fenster pro ID
%           fig=1 ... immer in das Fenster derselben ID
%           fig>1 ... setzt maximale Zahl der Grafiken derselben ID
%           fig<0 ... immer neues Fenster
%           fig=0 ... keine Grafikausgabe
% algo     (struct, default mfilename der aufrufenden Funktion)  Info uber 
%          Algorithmus, Funktion oder Klasse.
%   .mfilename   (string, opt) mfilename oder sonstiger Name des Algorithmus
%   .version     (string, opt) Versionsinfo
%   .versiondate (string, opt) Datum der Version
%   .variant     (string, opt) Variante
% opt (struct)
%   .sizeX/Y     (double, default 6/7) Anteil der Fenstergröße am Bildschirm
%   .pixsizeX/Y  (int) statt rel. groesse, groesse in pixeln
%

if nargin <1
    fig=1;
end

if fig==0
    fighandle=0;
    return;
end

if nargin <2 || isempty(algo) || ~isfield(algo,'mfilename')
    try
        ST=dbstack;
        algo.mfilename=ST(2).file;       
    catch exception
        algo.mfilename='';        
    end   
end
if nargin<3 || isempty(opt)
    opt=struct;
end

if ~isfield(algo,'version')
    algo.version='';
else
    algo.version=[' V. ',num2str(algo.version)];
end
if ~isfield(algo, 'versiondate')
    algo.versiondate=[];
else
    algo.versiondate=['  vom ',algo.versiondate];
end
if ~isfield(algo, 'variant')
    algo.variant=[];
else
   algo.variant=[', ',num2str(algo.variant)];
end

% Options
if ~isfield(opt,'pixsizeX')
    opt.pixsizeX=[];
end
if ~isfield(opt,'pixsizeY')
    opt.pixsizeY=[];
end
if ~isfield(opt,'sizeX') || isempty(opt.sizeX)
    opt.sizeX=6/7;
end
if ~isfield(opt,'sizeY') || isempty(opt.sizeY)
    opt.sizeY=6/7;
end
 
% id fuer die offene Fenster gezaehlt werden:
figurename_gen=[algo.mfilename,algo.version,algo.versiondate,algo.variant];
figurename_reg=regexptranslate('escape',figurename_gen);

ScreenSize = get(0,'ScreenSize');
if isempty(opt.pixsizeX)
    opt.pixsizeX=ScreenSize(3)*opt.sizeX;
end
if isempty(opt.pixsizeY)
    opt.pixsizeY=ScreenSize(4)/ScreenSize(3)*opt.pixsizeX;
end

fig_timestamp=datestr(clock);
hh= datevec(fig_timestamp);
% verwende als tag Sek nach Mitternacht
fig_timetag=num2str(hh(4)*3600+hh(5)*60+hh(6),'%05i');
figurename=[figurename_gen,' @ ', fig_timestamp];


h=findobj('-regexp','Name',['(',figurename_reg,').*']);
h_timetags=get(h,'Tag');
fighandle=0;

if isempty(h) || (~isempty(fig) && (fig<0 || length(h)<fig))
    % neues Fenster
    h = figure('Name',figurename,'Color', 'White',...
        'Position',round([190,50,opt.pixsizeX,opt.pixsizeY]));
    fighandle=h;
elseif isempty(fig) || fig~=0
    % ueberschreibe aeltestes Fenster
    if iscell(h_timetags)
        [oldest_tag,idx]= min(cellfun(@str2num,h_timetags));
    else
        idx=1;
    end
    
    fighandle=h(idx);
end

if fighandle>0
    set(fighandle,'Name',figurename);
    set(fighandle,'Tag',fig_timetag);
    set(0,'currentFigure',fighandle); 
    figure(fighandle);
    clf;
end

end

