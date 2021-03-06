function edit_photoblParam_03_Callback(obj, evd, h)
p = h.param.ttPr;
if ~isempty(p.proj)
    proj = p.curr_proj;
    mol = p.curr_mol(proj);
    method = p.proj{proj}.curr{mol}{2}{1}(2);
    if method == 2 % Threshold
        val = str2num(get(obj, 'String'));
        inSec = p.proj{proj}.fix{2}(7);
        nExc = p.proj{proj}.nb_excitations;
        len = nExc*size(p.proj{proj}.intensities,1);
        start = p.proj{proj}.prm{mol}{2}{1}(4);
        rate = p.proj{proj}.frame_rate;
        if inSec
            val = rate*round(val/rate);
            minVal = rate*start;
            maxVal = rate*len;
        else
            val = round(val);
            minVal = start;
            maxVal = len;
        end
        set(obj, 'String', num2str(val));

        if ~(~isempty(val) && numel(val) == 1 && ~isnan(val) && ...
                val >= minVal && val <= maxVal)
            set(obj, 'BackgroundColor', [1 0.75 0.75]);
            updateActPan(['Minimum cutoff must be >= ' num2str(minVal) ...
                ' and <= ' num2str(maxVal)], h.figure_MASH, 'error');

        else
            set(obj, 'BackgroundColor', [1 1 1]);
            if inSec
                val = val/rate;
            end
            chan = p.proj{proj}.curr{mol}{2}{1}(3);
            p.proj{proj}.curr{mol}{2}{2}(chan,3) = val;
            h.param.ttPr = p;
            guidata(h.figure_MASH, h);
            ud_bleach(h.figure_MASH);
        end
    end
end