function popupmenu_TDPcolour_Callback(obj, evd, h)
p = h.param.TDP;
if ~isempty(p.proj)
    val = get(obj, 'Value');
    proj = p.curr_proj;
    tpe = p.curr_type(proj);
    trans = p.proj{proj}.prm{tpe}.clst_start{1}(4);
    p.proj{proj}.prm{tpe}.clst_start{3}(trans,:) = p.colList(val,:);
    h.param.TDP = p;
    guidata(h.figure_MASH, h);
    updateFields(h.figure_MASH, 'TDP');
end