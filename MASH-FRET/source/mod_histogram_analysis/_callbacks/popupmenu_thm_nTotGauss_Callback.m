function popupmenu_thm_nTotGauss_Callback(obj, evd, h_fig)
h = guidata(h_fig);
p = h.param.thm;
if ~isempty(p.proj)
    updateFields(h_fig, 'thm');
end