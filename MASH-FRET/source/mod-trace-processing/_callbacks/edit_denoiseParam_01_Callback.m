function edit_denoiseParam_01_Callback(obj, evd, h)
p = h.param.ttPr;
if ~isempty(p.proj)
    val = round(str2num(get(obj, 'String')));
    set(obj, 'String', num2str(val));
    proj = p.curr_proj;
    mol = p.curr_mol(proj);
    method = p.proj{proj}.curr{mol}{1}{1}(1);

    if ~(~isempty(val) && numel(val) == 1 && ~isnan(val))
        set(obj, 'BackgroundColor', [1 0.75 0.75]);
        updateActPan('Denoising parameters must be >= 0', ...
            h.figure_MASH, 'error');

    else
        switch method

            case 1 % Rolling point: running average window size
                nmax = p.proj{proj}.nb_excitations * ...
                    size(p.proj{proj}.intensities,1);
                if ~(val > 0 && val < nmax)
                    set(obj, 'BackgroundColor', [1 0.75 0.75]);
                    updateActPan(['Running average window size must be' ...
                        ' > 0 and < ' num2str(nmax)], h.figure_MASH, ...
                        'error');
                    return;
                end

            case 2 % Haran filter: exponential factor for predictor weights

            case 3 % Wavelet analysis: shrink strength, firm/hard/soft
                if ~(sum(double(val == [1 2 3])))
                    set(obj, 'BackgroundColor', [1 0.75 0.75]);
                    updateActPan(['Shrink strength must be 1, 2 or 3 ' ...
                        '(firm, hard or soft)'], h.figure_MASH, 'error');
                    return;
                end

            otherwise
                return;
        end
        set(obj, 'BackgroundColor', [1 1 1]);
        p.proj{proj}.curr{mol}{1}{2}(method,1) = val;
        h.param.ttPr = p;
        guidata(h.figure_MASH, h);
        ud_denoising(h.figure_MASH);
    end
end