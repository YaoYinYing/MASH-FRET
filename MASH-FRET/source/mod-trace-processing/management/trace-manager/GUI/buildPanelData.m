function q = buildPanelData(q,p,h_fig)

% default
fact = 5;
ud = 1;
str0 = 'x-data';
str1 = 'x-value';
str2 = {'original time traces (/frame)','state trajectories (/frame)',...
    'means (/mol.)','minima (/mol.)','maxima (/mol.)','medians (/mol.)',...
    'nb. of states (/mol.)','nb. of transitions (/mol.)',...
    'state values (/state)','state lifetime (/state)'};
str3 = 'x min';
str4 = 'x max';
str5 = 'nbins';
str6 = 'y-data';
str7 = 'y-value';
str8 = 'y min';
str9 = 'y max';
str10 = 'nbins';
ttstr0 = 'Select the data to represent on the x-axis of the histogram';
ttstr1 = 'Select the value to histogram on the x-axis of the histogram';
ttstr2 = 'Lower bound of x-axis';
ttstr3 = 'Upper bound of x-axis';
ttstr4 = 'Number of binning intervals in x-axis';
ttstr5 = 'Select the data to represent on the y-axis of the 2D-histogram';
ttstr6 = 'Select the value to histogram on the y-axis of the 2D-histogram';
ttstr7 = 'Lower bound of y-axis';
ttstr8 = 'Upper bound of y-axis';
ttstr9 = 'Number of binning intervals in y-axis';

% parent
h_pan = q.uipanel_data;

% dimensions
pospan = get(h_pan,'position');
wpop2 = (pospan(3)-2*p.mg-p.mg/fact)/2;
wedit2 = (pospan(3)-2*p.mg-2*p.mg/fact)/3;

% list strings
str_pop = getStrPlot_overall(h_fig);

x = p.mg;
y = pospan(4)-p.mgbig-p.htxt;

q.text_xdata = uicontrol('style','text','parent',h_pan,'units',p.posun,...
    'position',[x,y,wpop2,p.htxt],'fontunits',p.fntun,'fontsize',p.fntsz,...
    'string',str0);

x = x+wpop2+p.mg/fact;

q.text_xval = uicontrol('style','text','parent',h_pan,'units',p.posun,...
    'position',[x,y,wpop2,p.htxt],'fontunits',p.fntun,'fontsize',p.fntsz,...
    'string',str1);

x = p.mg;
y = y-p.hpop;

q.popupmenu_selectXdata = uicontrol('style','popupmenu','parent',h_pan,...
    'units',p.posun,'position',[x,y,wpop2,p.hpop],'fontunits',p.fntun,...
    'fontsize',p.fntsz,'string',str_pop{2},'tooltipstring',ttstr0,...
    'callback',{@popupmenu_selectData_Callback,h_fig});

x = x+wpop2+p.mg/fact;
 
q.popupmenu_selectXval = uicontrol('style','popupmenu','parent',h_pan,...
    'units',p.posun,'position',[x,y,wpop2,p.hpop],'fontunits',p.fntun,...
    'fontsize',p.fntsz,'string',str2,'tooltipstring',ttstr1,'callback',...
    {@popupmenu_selectData_Callback,h_fig},'userdata',ud);

x = p.mg;
y = y-p.mg/2-p.htxt;

q.text_xlow = uicontrol('style','text','parent',h_pan,'units',p.posun,...
    'string',str3,'position',[x,y,wedit2,p.htxt],'fontunits',p.fntun,...
    'fontsize',p.fntsz);

x = x+wedit2+p.mg/fact;

q.text_xup = uicontrol('style','text','parent',h_pan,'units',p.posun,...
    'string',str4,'position',[x,y,wedit2,p.htxt],'fontunits',p.fntun,...
    'fontsize',p.fntsz);

x = x+wedit2+p.mg/fact;

q.text_xniv = uicontrol('style','text','parent',h_pan,'units',p.posun,...
    'string',str5,'position',[x,y,wedit2,p.htxt],'fontunits',p.fntun,...
    'fontsize',p.fntsz);

x = p.mg;
y = y-p.hedit;

q.edit_xlow = uicontrol('style','edit','parent',h_pan,'units',p.posun,...
    'position',[x,y,wedit2,p.hedit],'fontunits',p.fntun,'fontsize',p.fntsz,...
    'string','','tooltipstring',ttstr2,'callback',...
    {@edit_xlow_Callback,h_fig});

x = x+wedit2+p.mg/fact;

q.edit_xup = uicontrol('style','edit','parent',h_pan,'units',p.posun,...
    'position',[x,y,wedit2,p.hedit],'fontunits',p.fntun,'fontsize',...
    p.fntsz,'string','','tooltipstring',ttstr3,'callback',...
    {@edit_xup_Callback,h_fig});

x = x+wedit2+p.mg/fact;

q.edit_xniv = uicontrol('style','edit','parent',h_pan,'units',p.posun,...
    'string','','tooltipstring',ttstr4,'position',[x,y,wedit2,p.hedit],...
    'fontunits',p.fntun,'fontsize',p.fntsz,'callback',...
    {@edit_xniv_Callback,h_fig});

y = y-p.mg-p.htxt;
x = p.mg;

q.text_ydata = uicontrol('style','text','parent',h_pan,'units',p.posun,...
    'position',[x,y,wpop2,p.htxt],'fontunits',p.fntun,'fontsize',p.fntsz,...
    'string',str6);

x = x+wpop2+p.mg/fact;

q.text_yval = uicontrol('style','text','parent',h_pan,'units',p.posun,...
    'position',[x,y,wpop2,p.htxt],'fontunits',p.fntun,'fontsize',p.fntsz,...
    'string',str7);

x = p.mg;
y = y-p.hpop;

q.popupmenu_selectYdata = uicontrol('style','popupmenu','parent',h_pan,...
    'units',p.posun,'position',[x,y,wpop2,p.hpop],'fontunits',p.fntun,...
    'fontsize',p.fntsz,'string',str_pop{2},'tooltipstring',ttstr5,...
    'callback',{@popupmenu_selectData_Callback,h_fig});

x = x+wpop2+p.mg/fact;
 
q.popupmenu_selectYval = uicontrol('style','popupmenu','parent',...
    h_pan,'units',p.posun,'string',str2,'tooltipstring',ttstr6,'position',...
    [x,y,wpop2,p.hpop],'fontunits',p.fntun,'fontsize',p.fntsz,'callback',...
    {@popupmenu_selectData_Callback,h_fig},'userdata',ud);

x = p.mg;
y = y-p.mg/2-p.htxt;

q.text_ylow = uicontrol('style','text','parent',h_pan,'units',p.posun,...
    'string',str8,'position',[x,y,wedit2,p.htxt],'fontunits',p.fntun,...
    'fontsize',p.fntsz,'enable','off');

x = x+wedit2+p.mg/fact;

q.text_yup = uicontrol('style','text','parent',h_pan,'units',p.posun,...
    'string',str9,'position',[x,y,wedit2,p.htxt],'fontunits',p.fntun,...
    'fontsize',p.fntsz,'enable','off');

x = x+wedit2+p.mg/fact;

q.text_yniv = uicontrol('style','text','parent',h_pan,'units',p.posun,...
    'string',str10,'position',[x,y,wedit2,p.htxt],'fontunits',p.fntun,...
    'fontsize',p.fntsz,'enable','off');

x = p.mg;
y = y-p.hedit;

q.edit_ylow = uicontrol('style','edit','parent',h_pan,'units',p.posun,...
    'string','','tooltipstring',ttstr7,'position',[x,y wedit2,p.hedit],...
    'fontunits',p.fntun,'fontsize',p.fntsz,'callback',...
    {@edit_ylow_Callback,h_fig},'enable','off');

x = x+wedit2+p.mg/fact;

q.edit_yup = uicontrol('style','edit','parent',h_pan,'units',p.posun,...
    'string','','tooltipstring',ttstr8,'position',[x,y,wedit2,p.hedit],...
    'fontunits',p.fntun,'fontsize',p.fntsz,'callback',...
    {@edit_yup_Callback,h_fig},'enable','off');

x = x+wedit2+p.mg/fact;

q.edit_yniv = uicontrol('style','edit','parent',h_pan,'units',p.posun,...
    'string','','tooltipstring',ttstr9,'position',[x,y,wedit2,p.hedit],...
    'fontunits',p.fntun,'fontsize',p.fntsz,'enable','off','callback',...
    {@edit_yniv_Callback,h_fig});


