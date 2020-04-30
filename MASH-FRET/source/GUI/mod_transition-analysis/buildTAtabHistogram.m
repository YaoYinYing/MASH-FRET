function h = buildTAtabHistogram(h,p)
% h = buildTAtabHistogram(h,p)
%
% Builds tab "Histogram" in "Kinetic model" panel.
%
% h: structure to update with handles to new UI components and that must contain fields:
%   h.figure_MASH: handle to main figure
%   h.uitab_TA_histogram: handle to tab "Histogram"

% defaults
hpop0 = 22;
str0 = {'simulation','experimental','overlay'};
ttstr0 = wrapHtmlTooltipString('<b>Select data</b> to display in the histogram plot');

% parents
h_fig = h.figure_MASH;
h_tab = h.uitab_TA_histogram;

% dimensions
postab = get(h_tab,'position');
wpop0 = postab(3)-2*p.mg;
waxes0 = postab(3)-2*p.mg;
haxes0 = postab(4)-3*p.mg-hpop0;

% GUI
x = p.mg;
y = postab(4)-p.mg-hpop0;

h.popupmenu_TA_mdlHistDat = uicontrol('style','popupmenu','parent',h_tab,...
    'units',p.posun,'fontunits',p.fntun,'fontsize',p.fntsz1,'position',...
    [x,y,wpop0,hpop0],'string',str0,'tooltipstring',ttstr0,'callback',...
    {@popupmenu_TA_mdlDat_Callback,h_fig});

x = p.mg;
y = y-p.mg-haxes0;

h.axes_TA_mdlHist = axes('parent',h_tab,'units',p.posun,'fontunits',p.fntun,...
    'fontsize',p.fntsz1,'position',[x,y,waxes0,haxes0],'nextplot',...
    'replacechildren','box','on','xtick',[],'ytick',[]);

