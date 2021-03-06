function checkbox_thm_BS_Callback(obj, evd, h)
p = h.param.thm;
if ~isempty(p.proj)
    proj = p.curr_proj;
    tpe = p.curr_tpe(proj);
    p.proj{proj}.prm{tpe}.thm_start{1}(2) = get(obj, 'Value');
    p.proj{proj}.prm{tpe}.thm_res(1,1:3) = {[] [] []};
    p.proj{proj}.prm{tpe}.thm_res(2,1:3) = {[] [] []};
    h.param.thm = p;
    guidata(h.figure_MASH, h);
    updateFields(h.figure_MASH, 'thm');
end