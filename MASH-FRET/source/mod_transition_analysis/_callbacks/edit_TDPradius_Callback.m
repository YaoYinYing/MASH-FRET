function edit_TDPradius_Callback(obj, evd, h)
p = h.param.TDP;
if ~isempty(p.proj)
    val = str2num(get(obj, 'String'));
    set(obj, 'String', num2str(val));
    proj = p.curr_proj;
    tpe = p.curr_type(proj);
    meth = p.proj{proj}.prm{tpe}.clst_start{1}(1);

    if ~(numel(val)==1 && ~isnan(val) && val >= 0)
        set(obj, 'BackgroundColor', [1 0.75 0.75]);
        switch meth
            case 1 % kmean
                str = 'Tolerance radii';
            case 2 % gaussian model based
                str = 'Min. Gaussian standard deviation.';
        end
        setContPan([str ' must be >= 0.'], 'error', h.figure_MASH);

    else
        switch meth
            case 1 % kmean
                state = get(h.popupmenu_TDPstate, 'Value');
                p.proj{proj}.prm{tpe}.clst_start{2}(state,2) = val;

            case 2 % gaussian model based
                p.proj{proj}.prm{tpe}.clst_start{2}(1,2) = val;
        end
        h.param.TDP = p;
        guidata(h.figure_MASH, h);
        updateFields(h.figure_MASH, 'TDP');
    end
end
