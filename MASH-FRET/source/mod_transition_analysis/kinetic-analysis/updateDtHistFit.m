function p = updateDtHistFit(p,tag,tpe,v,h_fig)

% collect interface parameters
proj = p.curr_proj;

% collect prcoessing parameters
prm = p.proj{proj}.prm{tag,tpe};
curr = p.proj{proj}.curr{tag,tpe};

% make current settings the last applied settings
prm.lft_start = curr.lft_start;
prm.lft_res = curr.lft_res;

J = prm.lft_start{2}(1);
bin = prm.lft_start{2}(3);
excl = prm.lft_start{2}(4);
rearr = prm.lft_start{2}(5);
lft_k = prm.lft_start{1}(v,:);
stchExp = lft_k{1}(2);
boba = lft_k{1}(5);
p_boba = [];
if boba
    p_boba = lft_k{1}([6 7 8]);
end
mat = prm.clst_start{1}(4);
clstDiag = prm.clst_start{1}(9);
dat = prm.clst_res{1}.clusters{J};
stateVals = prm.clst_res{1}.mu{J};

% reset results
prm.lft_res(v,:) = p.proj{proj}.def{tag,tpe}.lft_res;

% get fit settings
if stchExp
    % amp, dec, beta
    p_fit.lower = lft_k{2}(1,[1 4 7]);
    p_fit.start = lft_k{2}(1,[2 5 8]);
    p_fit.upper = lft_k{2}(1,[3 6 9]);
else
    % amp1, dec1, amp2, dec2 ...
    p_fit.lower = reshape(lft_k{2}(:,[1 4])', [1 ...
        numel(lft_k{2}(:,[1 4]))]);
    p_fit.start = reshape(lft_k{2}(:,[2 5])', [1 ...
        numel(lft_k{2}(:,[2 5]))]);
    p_fit.upper = reshape(lft_k{2}(:,[3 6])', [1 ...
        numel(lft_k{2}(:,[3 6]))]);
end

% bin state values
nTrs = getClusterNb(J,mat,clstDiag);
[j1,j2] = getStatesFromTransIndexes(1:nTrs,J,mat,clstDiag);
[stateVals,js] = binStateValues(stateVals,bin,[j1,j2]);
V = numel(stateVals);
for val = 1:V
    for j = 1:numel(js{val})
        dat(dat(:,end-1)==js{val}(j),end-1) = val;
        dat(dat(:,end-1)==js{val}(j),end-3) = stateVals(val);
        dat(dat(:,end)==js{val}(j),end) = val;
        dat(dat(:,end)==js{val}(j),end-2) = stateVals(val);
    end
end

% re-arrange state sequences by cancelling transitions belonging to diagonal clusters
if rearr
    [mols,o,o] = unique(dat(:,4));
    dat_new = [];
    for m = mols'
        dat_m = dat(dat(:,4)==m,:);
        if isempty(dat_m)
            continue
        end
        dat_m = adjustDt(dat_m);
        if size(dat_m,1)==1
            continue
        end
        dat_new = cat(1,dat_new,dat_m);
    end
    dat = dat_new;
end

% histogram fitting
setContPan('Fitting in progress ...', 'process', h_fig);
res = fitDt(dat, v, excl, prm.clst_res{4}{v}, p_fit, p_boba, h_fig);
if isempty(res)
    return
end

if boba
    % update number of replicates
    prm.lft_start{1}{v,1}(6) = res.n_rep;
    prm.lft_res{v,1} = res.boba_mean;
    prm.lft_res{v,3} = res.boba_inf;
    prm.lft_res{v,4} = res.boba_sup;
    prm.lft_res{v,5} = {res.histspl,res.boba_fitres};
end

prm.lft_res{v,2} = res.fit_ref;

% update modifications of processing parameters to current settings
curr.lft_start = prm.lft_start;
curr.lft_res = prm.lft_res;

% save modifications
p.proj{proj}.prm{tag,tpe} = prm;
p.proj{proj}.curr{tag,tpe} = curr;
