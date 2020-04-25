function h = buildPanelTAkineticModel(h,p)
% h = buildPanelTAkineticModel(h,p)
%
% Builds panel "Kinetic model" in "Transition analysis" module.
%
% h: structure to update with handles to new UI components and that must contain fields:
%   h.figure_MASH: handle to main figure
%   h.uipanel_TA_kineticModel: handle to panel "Kinetic model"

% defaults
htxt0 = 14;
hbut0 = 20;
str0 = 'Refresh model';
str1 = 'Comparison with simulation:';
ttl0 = 'TDP';
ttl1 = 'Histogram';
ttl2 = 'Dwell times';
ttstr0 = wrapHtmlTooltipString('<b>Repeat simulation</b>');

% parent
h_fig = h.figure_MASH;
h_pan = h.uipanel_TA_kineticModel;

% dimensions
pospan = get(h_pan,'position');
wbut0 = getUItextWidth(str0,p.fntun,p.fntsz1,'normal',p.tbl)+p.wbrd;
wtxt0 = getUItextWidth(str1,p.fntun,p.fntsz1,'normal',p.tbl);
waxes0 = (pospan(3)-3*p.mg)/2;
haxes0 = pospan(4)-p.mgpan-hbut0-2*p.mg;
wtab = waxes0;
htab = pospan(4)-p.mgpan-htxt0-2*p.mg;

x = p.mg;
y = pospan(4)-p.mgpan-hbut0;

h.pushbutton_TA_refreshModel = uicontrol('style','pushbutton','parent',...
    h_pan,'units',p.posun,'fontunits',p.fntun,'fontsize',p.fntsz1,...
    'position',[x,y,wbut0,hbut0],'string',str0,'tooltipstring',ttstr0,...
    'callback',{@pushbutton_TA_refreshModel_Callback,h_fig});

y = y-p.mg-haxes0;

h.axes_TDPplot3 = axes('parent',h_pan,'units',p.posun,'fontunits',p.fntun,...
    'fontsize',p.fntsz1,'position',[x,y,waxes0,haxes0],'visible','off');

x = p.mg+waxes0+p.mg;
y = pospan(4)-p.mgpan-htxt0;

h.text_TA_simCompare = uicontrol('style','text','parent',h_pan,'units',...
    p.posun,'fontunits',p.fntun,'fontsize',p.fntsz1,'position',...
    [x,y,wtxt0,htxt0],'string',str1,'horizontalalignment','left');

y = y-p.mg-htab;

h.uitabgroup_TA_simModel = uitabgroup('parent',h_pan,'units',p.posun,...
    'position',[x,y,wtab,htab]);
h_tabgrp = h.uitabgroup_TA_simModel;

h.uitab_TA_tdp = uitab('parent',h_tabgrp,'units',p.posun,'title',ttl0);

h.uitab_TA_hist = uitab('parent',h_tabgrp,'units',p.posun,'title',ttl1);

h.uitab_TA_dwelltimes = uitab('parent',h_tabgrp,'units',p.posun,'title',...
    ttl2);

