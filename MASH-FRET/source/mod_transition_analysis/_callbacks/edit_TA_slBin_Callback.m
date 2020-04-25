function edit_TA_slBin_Callback(obj,evd,h_fig)

h = guidata(h_fig);
p = h.param.TDP;
if isempty(p.proj)
    return
end

val = round(str2num(get(obj, 'String')));
set(obj, 'String', num2str(val));
if ~(numel(val)==1 && ~isnan(val) && val>=0)
    set(obj, 'BackgroundColor', [1 0.75 0.75]);
    setContPan('State binning must be null or positive','error',h_fig);
    return
end

proj = p.curr_proj;
tpe = p.curr_type(proj);
tag = p.curr_tag(proj);
prm = p.proj{proj}.prm{tag,tpe};
curr = p.proj{proj}.curr{tag,tpe};

curr.lft_start{2}(3) = val;
prm.lft_start{2}(3) = curr.lft_start{2}(3);

% recalculate histograms
V = size(prm.clst_res{4},2);
prm.clst_res{4} = cell(1,V);
for v = 1:V
    prm.clst_res{4}{v} = updateDtHist(prm,[],v);
end
curr.clst_res{4} = prm.clst_res{4};

p.proj{proj}.prm{tag,tpe} = prm;
p.proj{proj}.curr{tag,tpe} = curr;
h.param.TDP = p;
guidata(h_fig, h);

