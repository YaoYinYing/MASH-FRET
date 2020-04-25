function checkbox_TA_slExcl_Callback(obj,evd,h_fig)

h = guidata(h_fig);
p = h.param.TDP;
if isempty(p.proj)
    return
end

proj = p.curr_proj;
tpe = p.curr_type(proj);
tag = p.curr_tag(proj);
prm = p.proj{proj}.prm{tag,tpe};
curr = p.proj{proj}.curr{tag,tpe};

curr.lft_start{2}(4) = get(obj, 'Value');
prm.lft_start{2}(4) = curr.lft_start{2}(4);

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

ud_kinFit(h_fig);
updateTAplots(h_fig,'kin');
