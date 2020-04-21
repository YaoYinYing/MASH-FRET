function h = buildContextmenu_proj(h)

% default
lbl0 = 'Edit...';
lbl1 = 'Close';
lbl2 = 'Save As...';
lbl3 = 'Export ASCII files...';
lbl4 = 'Merge projects';

h_fig = h.figure_MASH;

% build and add context menus
h.cm_proj = uicontextmenu('parent',h_fig);
h_cm = h.cm_proj;
h.projMenu_edit = uimenu('parent',h_cm,'label',lbl0,'callback',...
    {@pushbutton_editParam_Callback,h_fig});
h.projMenu_close = uimenu('parent',h_cm,'label',lbl1,'callback',...
    {@pushbutton_remTraces_Callback,h_fig});
h.projMenu_save = uimenu('parent',h_cm,'label',lbl2,'callback',...
    {@pushbutton_expProj_Callback,h_fig});
h.projMenu_export= uimenu('parent',h_cm,'label',lbl3,'callback',...
    {@menu_projMenu_export_Callback,h_fig});
h.projMenu_merge = uimenu('parent',h_cm,'label',lbl4,'callback',...
    {@menu_projMenu_merge_Callback,h_fig});

set(h.listbox_traceSet,'uicontextmenu',h_cm);

