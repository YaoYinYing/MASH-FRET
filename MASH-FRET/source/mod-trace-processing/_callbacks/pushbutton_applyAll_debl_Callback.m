function pushbutton_applyAll_debl_Callback(obj, evd, h)
p = h.param.ttPr;
if ~isempty(p.proj)
    choice = questdlg( {['Overwriting debleaching parameters of all ' ...
        'molecules erases previous traces processing'], ...
        'Overwrite debleaching parameters of all molecule?'}, ...
        'Overwrite parameters', 'Apply', 'Cancel', 'Apply');

    if strcmp(choice, 'Apply')
        proj = p.curr_proj;
        mol = p.curr_mol(proj);
        nMol = size(p.proj{proj}.coord_incl,2);
        for m = 1:nMol
            p.proj{proj}.curr{m}{2} = p.proj{proj}.curr{mol}{2};
        end
        p.proj{proj}.def.mol{2} = p.proj{proj}.curr{mol}{2};
        h.param.ttPr = p;
        guidata(h.figure_MASH, h);
    end
end