function p = updateDtHistFit(p,tag,tpe,v,h_fig,lb)

% defaults
nSpl = 100;

% collect interface parameters
proj = p.curr_proj;

% collect processing parameters
prm = p.proj{proj}.prm{tag,tpe};
curr = p.proj{proj}.curr{tag,tpe};
def = p.proj{proj}.def{tag,tpe};

% make current settings the last applied settings
prm.lft_start = curr.lft_start;
prm.lft_res = curr.lft_res;

J = prm.lft_start{2}(1);
bin = prm.lft_start{2}(3);
excl = prm.lft_start{2}(4);
rearr = prm.lft_start{2}(5);
auto = prm.lft_start{1}{v,1}(1);
mat = prm.clst_start{1}(4);
clstDiag = prm.clst_start{1}(9);
dat = prm.clst_res{1}.clusters{J};
stateVals = prm.clst_res{1}.mu{J};

% bin state values
nTrs = getClusterNb(J,mat,clstDiag);
[j1,j2] = getStatesFromTransIndexes(1:nTrs,J,mat,clstDiag);
[stateVals,js] = binStateValues(stateVals,bin,[j1,j2]);
V = numel(stateVals);
dat_new = dat;
for val = 1:V
    for j = 1:numel(js{val})
        dat_new(dat(:,end-1)==js{val}(j),end-1) = val;
        dat_new(dat(:,end)==js{val}(j),end) = val;
    end
end
dat = dat_new;

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

% reset results
prm.lft_res(v,:) = def.lft_res;
prm.mdl_res = def.mdl_res;

% model selection
if auto
    [res,ok] = getOptExpFitModel(dat,v,excl,prm.clst_res{4}{v},nSpl,h_fig,...
        lb);
    if ~ok || isempty(res)
        return
    end
    
    fitprm = res.fit_ref;
    nExp = size(fitprm,1);
    boba = true;
    
    prm.lft_start{1}{v,1}(2) = false; % stretched exp
    prm.lft_start{1}{v,1}(3) = nExp; % exp nb.
    prm.lft_start{1}{v,1}(5) = boba; % boba
    prm.lft_start{1}{v,1}(7) = nSpl; % sample number
    prm.lft_start{1}{v,1}(8) = false; % weighing
    
    % lower, start, upper
    prm.lft_start{1}{v,2} = [zeros(nExp,6),...
        repmat(prm.lft_start{1}{v,2}(1,7:9),[nExp,1])];
    prm.lft_start{1}{v,2}(:,[2 5]) = fitprm;
    prm.lft_start{1}{v,2}(:,[3 6]) = Inf(nExp,2);
    
else
    % collect parameters
    lft_k = prm.lft_start{1}(v,:);
    stchExp = lft_k{1}(2);
    boba = lft_k{1}(5);
    p_boba = [];
    if boba
        p_boba = lft_k{1}([6 7 8]);
    end

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

    % histogram fitting
    setContPan('Fitting in progress ...', 'process', h_fig);
    res = fitDt(dat,v,excl,prm.clst_res{4}{v},p_fit,p_boba,h_fig,lb);
    if isempty(res)
        return
    end
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
curr.mdl_res = prm.mdl_res;

% save modifications
p.proj{proj}.prm{tag,tpe} = prm;
p.proj{proj}.curr{tag,tpe} = curr;
