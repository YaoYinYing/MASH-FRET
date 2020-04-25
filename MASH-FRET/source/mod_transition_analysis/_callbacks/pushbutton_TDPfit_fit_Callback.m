function pushbutton_TDPfit_fit_Callback(obj, evd, h_fig)
% pushbutton_TDPfit_fit_Callback([],[],h_fig)
%
% h_fig: handle to main figure

% Last update by MH, 27.1.2020: move fitting script to separate function updateDtHistFit.m and plot dwell time histogram after fitting

% get interface parameters
h = guidata(h_fig);
p = h.param.TDP;
if isempty(p.proj)
    return
end

proj = p.curr_proj;
tpe = p.curr_type(proj);
tag = p.curr_tag(proj);
curr = p.proj{proj}.curr{tag,tpe};

% get processing parameters
v = curr.lft_start{2}(2);

% update histogram and fit
p = updateDtHistFit(p,tag,tpe,v,h_fig);

% save results
h.param.TDP = p;
guidata(h_fig,h);

% update plots and GUI
updateFields(h_fig, 'TDP');
