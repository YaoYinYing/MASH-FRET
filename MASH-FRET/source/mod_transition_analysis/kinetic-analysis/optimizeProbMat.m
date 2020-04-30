function [mat_iter,mat_err,simdat] = optimizeProbMat(r,states,expPrm,opt,mat0)
% [mat,mat_err,simdat] = optimizeProbMat(r,states,expPrm,opt,mat0)
%
% Find the probability matrix that describes the input dwell time set
%
% r: [1-by-J] restricted rate matrix (second-1)
% states: [1-by-J] state values
% expPrm: experimental parameters used in simulation with fields:
%   expPrm.dt: [nDt-by-3] dwell times (seconds) and state indexes before and after transition
%   expPrm.w: [1-by-J] relative total number of transitions from each state
%   expPrm.Ls: [1-by-N] experimental trajectory lengths
%   expPrm.expT: binning time
%   expPrm.fitPrm: fitting parameters
%   expPrm.excl: (1) to exclude first and last dwell times of each sequence, (0) otherwise
%   expPrm.clstPop: [V-by-V] TDP cluster relative populations
% opt: 'n','w' or 'tp' to optimize the matrix of number of transitions (fixed restricted rates), the repartition probability matrix (fixed restricted rates) or the transition probability matrix respectively
% mat0: [J-by-J] matrix starting guess
% mat: [J-by-J] best inferrence matrix
% mat_err: [J-by-J-by-2] negative and positive absolute mean deviation
% simdat: structure containing simulated data with fields:
%  simdat.dt: [nDt-by-4] dwell times (molecule index, state durations (frames), state values before and after transition)
%  simdat.tp_exp: [J-by-J] transition probability matrix calculated from simulated data
%  simdat.r_exp: [J-by-J] restricted rate matrix calculated from simulated data
%  simdat.w_exp: [J-by-J] weighing factors calculated from simulated data

% default
T = 5; % number of matrix initialization
nSpl = 3; % number of simulated data set used in error calculation
varsteps = 1./logspace(0,2,5); % variation step (relative)
plotIt = true;

% create figure for plot
if plotIt
   h_fig1 = figure('windowstyle','docked');
else
    h_fig1 = [];
end

% identify degenerated states
vals = unique(states);
V = numel(vals);
degen = cell(1,V);
for v = 1:V
    degen{v} = find(states==vals(v));
end

% build starting guess for transition probabilities
J = size(r,2);
if isempty(mat0)
    switch opt
        case 'w' % diagonal terms = 0
            mat0 = double(~eye(J));
            mat0 = mat0./repmat(sum(mat0,2),[1,J]);
        case 'tp' % diagonal terms > 0
            mat0 = zeros(J);
            mat0(~~eye(J)) = exp(-expPrm.expT*r);
            for j1 = 1:J
                j2s = true(1,J);
                j2s(j1) = false;
                mat0(j1,j2s) = (1-mat0(j1,j1))/(J-1);
            end
    end
end

if strcmp(opt,'n')
    % identify which number of transitions are fixed
    expPrm.isFixed = false(J);
    expPrm.isFixed(~~eye(J)) = true; % diagonal terms
    for j1 = 1:J
        if sum(states==states(j1))==1 % no degenerated states
            expPrm.isFixed(:,j1) = true;
            continue
        end
        v1 = find(vals==states(j1));
        for j2 = 1:J
            v2 = find(vals==states(j2));
            if expPrm.clstPop(v1,v2)==0 % no density of TDP
                expPrm.isFixed(v1,v2) = true;
            end
        end
    end

    % calculates and associates fixed sums
    expPrm.sumConstr = cell(J);
    degenTrans = false(J);
    for v = 1:V
        degenTrans(degen{v}(1):degen{v}(end),degen{v}(1):degen{v}(end)) = true;
    end
    for j1 = 1:J
        for j2 = 1:J
            expPrm.sumConstr{j1,j2} = zeros(J);
            if j1==j2
                continue
            end

            if states(j1)==states(j2) % transition between states of same value
                j_lvl = find(states==states(j1));
                j_lvl1 = min(j_lvl);
                j_lvl2 = max(j_lvl);
                nLvl = numel(j_lvl);
                sumDegen = Inf(nLvl); % no restriction on sum (=Inf)
                sumDegen(~~eye(nLvl)) = 0;
                expPrm.sumConstr{j1,j2}(j_lvl1:j_lvl2,j_lvl1:j_lvl2) = ...
                    sumDegen;
            else % transition between states of different values
                j_lvl = find(states==states(j2));
                j_lvl1 = min(j_lvl);
                j_lvl2 = max(j_lvl);
                expPrm.sumConstr{j1,j2}(j1,j_lvl1:j_lvl2) = ...
                    sum(mat0(j1,j_lvl1:j_lvl2));
                expPrm.sumConstr{j1,j2}(~degenTrans(:,j2),j2) = ...
                    sum(mat0(~degenTrans(:,j2),j2));
            end

            expPrm.sumConstr{j1,j2}(j1,j2) = 0;
        end
    end

    % convert variation step into number of transitions
    Ntot = sum(sum(mat0));
    varsteps = [logspace(2,0,5),varsteps(2:end)]*Ntot;
    varsteps = round([varsteps,1]);
end

% start probability variation, data simulation and data comparison cycles
t = tic;

ncycle = size(varsteps,2);
mat_all = cell(1,T);
gof_all = -Inf(1,T);
for restart = 1:T
    
    disp(['restart: ',num2str(restart),' ...']);
    
    if ~strcmp(opt,'n') && restart>1
        mat0 = rand(J);
        if strcmp(opt,'w')
            mat0(~~eye(J)) = 0;
        end
        mat0 = mat0./repmat(sum(mat0,2),[1,J]);
        
    elseif strcmp(opt,'n')
        % generate random number of transitions between same state values
        for v = 1:V
            mat0(degen{v}(1):degen{v}(end),degen{v}(1):degen{v}(end)) = ...
                rand(1)*Ntot;
        end
        mat0(~~eye(J)) = 0;
    end
    mat_iter = mat0;
    gof_prev = -Inf;
    bestgof = -Inf;
    bestgof_prev = -Inf;

    if t==1
        [gof_iter,simdat] = calcGOF(degen,mat_iter,r,expPrm,opt,h_fig1);
    else
        [gof_iter,simdat] = calcGOF(degen,mat_iter,r,expPrm,opt,[]);
    end
    switch opt
        case 'w'
            mat_iter = simdat.w_exp;
        case 'tp'
            mat_iter = simdat.tp_exp;
        case 'n'
            mat_iter = simdat.n_exp;
    end

    % repeat all variation steps
    while isinf(bestgof) || bestgof>bestgof_prev

        mat = cell(1,ncycle);
        gof = -Inf(1,ncycle);

        bestgof_prev = bestgof;
        mat_iter_prev = mat_iter;

        % increase variation step
        for cycle = 1:ncycle

            % repeat up and down variation of matrix elements
            while isinf(gof_prev) || gof_iter>gof_prev

                gof_prev = gof_iter;
                mat_prev = mat_iter;

                % vary each matrix element one after the other
                [mat_iter,gof_iter] = varyProb(mat_iter,gof_iter,1:J,degen,...
                    varsteps(cycle),r,expPrm,opt,[]);

                if gof_iter>gof_prev
                    if gof_iter>=max(gof) && gof_iter>=max(gof_all)
                        disp(['>> improvement: ',num2str(gof_iter)]);
                        disp(mat_iter);

                        % display best iteration
                        if plotIt
                            calcGOF(degen,mat_iter,r,expPrm,opt,h_fig1);
                        end
                    end
                else
                    if isinf(gof_iter) || isequal(mat_iter,mat_prev)
                        break
                    end
                    gof_iter = gof_prev;
                    mat_iter = mat_prev;
                end
            end

            gof(cycle) = gof_iter;
            mat{cycle} = mat_iter;
            gof_iter = -Inf;
            gof_prev = -Inf;
        end

        [bestgof,bestcycle] = max(gof);
        mat_iter = mat{bestcycle};

        if bestgof<=bestgof_prev
            bestgof = bestgof_prev;
            mat_iter = mat_iter_prev;
        end
    end

    mat_all{restart} = mat_iter;
    gof_all(restart) = bestgof;
end

[bestgof,bestrestart] = max(gof_all);
mat_iter = mat_all{bestrestart};

% calculate error range (3*sigma)
mat_spl = zeros(J,J,nSpl);
for s = 1:nSpl
    [gof,simdat] = calcGOF(degen,mat_iter,r,expPrm,opt,[]);
    switch opt
        case 'w'
            mat_spl(:,:,s) = simdat.w_exp;
        case 'tp'
            mat_spl(:,:,s) = simdat.tp_exp;
        case 'n'
            mat_spl(:,:,s) = simdat.n_exp;
    end
end
mat_iter = mean(mat_spl,3);
mat_err = 3*std(mat_spl,[],3);
disp(['best fit: ',num2str(bestgof)]);
disp(['processing time: ',num2str(toc(t))]);
disp(mat_iter);
disp(mat_err);


function [mat,bestgof] = varyProb(mat,gof,j1s,degen,step,r,expPrm,opt,h_fig)

% defaults
valmax = 1;

J = size(mat,2);
if strcmp(opt,'n')
    degenTrans = false(J);
    V = numel(degen);
    for v = 1:V
        degenTrans(degen{v}(1):degen{v}(end),degen{v}(1):degen{v}(end)) = true;
    end
    Ntot = sum(mat(~degenTrans)); % total number of transitions excluding transitions between degenerated states
end

alliter = {mat};
allgof = gof;
for j1 = j1s
    for j2 = 1:J
        if ~strcmp(opt,'tp') && j2==j1
            continue
        end
        if strcmp(opt,'n')
            valmax = min(min(...
                expPrm.sumConstr{j1,j2}(expPrm.sumConstr{j1,j2}>0)));
            if isinf(valmax) % transition between degenerated states
                valmax = 100*Ntot;
            end
        end

        if strcmp(opt,'tp') && j1==j2
            valmin = 0.1;
        elseif strcmp(opt,'tp')
            valmax = 0.9;
        else
            valmin = 0;
        end
        
        vardown = false;
        varup = false;

        % chose sens of variation (increment or decrement)
        if mat(j1,j2)>valmin
            switch opt
                case 'n' % decreases number of transitions
                    mat_down = incrConstRowColSumMat(mat,j1,j2,-step,...
                        expPrm.isFixed,expPrm.sumConstr); 
                otherwise % decreases transition probability
                    mat_down = setProb(mat,j1,j2,-step,opt); 
            end
            
            if ~isempty(mat_down)
                [gof_down,res] = ...
                    calcGOF(degen,mat_down,r,expPrm,opt,h_fig);
                switch opt
                    case 'w'
                        mat_down = res.w_exp;
                    case 'tp'
                        mat_down = res.tp_exp;
                    case 'n'
                        mat_down = res.n_exp;
                end
                if gof_down>gof
                    vardown = true;
                end
            end
        end
        if mat(j1,j2)<valmax
            switch opt
                case 'n' % increases number of transitions
                    mat_up = incrConstRowColSumMat(mat,j1,j2,step,...
                        expPrm.isFixed,expPrm.sumConstr); 
                otherwise % increases transition probability
                    mat_up = setProb(mat,j1,j2,step,opt); 
            end
            if ~isempty(mat_up)
                [gof_up,res] = calcGOF(degen,mat_up,r,expPrm,opt,h_fig);
                switch opt
                    case 'w'
                        mat_up = res.w_exp;
                    case 'tp'
                        mat_up = res.tp_exp;
                    case 'n'
                        mat_up = res.n_exp;
                end
                if gof_up>gof
                    varup = true;
                end
            end
        end
        
        % up iteration
        if varup
            gof_up_prev = gof;
            while gof_up>gof_up_prev

                gof_up_prev = gof_up;
                mat_up_prev = mat_up;
                if mat_up(j1,j2)==valmax
                    break
                end
                switch opt
                    case 'n' % increases number of transitions
                        mat_up = incrConstRowColSumMat(mat_up,j1,j2,step,...
                            expPrm.isFixed,expPrm.sumConstr);
                    otherwise % increases transition probability
                        mat_up = setProb(mat_up,j1,j2,step,opt);
                end
                if isempty(mat_up)
                    mat_up = mat_up_prev;
                    gof_up = gof_up_prev;
                    break
                end
                [gof_up,res] = calcGOF(degen,mat_up,r,expPrm,opt,h_fig);
                switch opt
                    case 'w'
                        mat_up = res.w_exp;
                    case 'tp'
                        mat_up = res.tp_exp;
                    case 'n'
                        mat_up = res.n_exp;
                end
                if gof_up<=gof_up_prev
                    mat_up = mat_up_prev;
                    gof_up = gof_up_prev;
                end
            end
            alliter = cat(2,alliter,mat_up);
            allgof = cat(2,allgof,gof_up);
        end

        % down iteration
        if vardown
            gof_down_prev = gof;
            while gof_down>gof_down_prev

                gof_down_prev = gof_down;
                mat_down_prev = mat_down;
                if mat_down(j1,j2)==valmin
                    break
                end
                switch opt
                    case 'n' % decreases number of transitions
                        mat_down = incrConstRowColSumMat(mat_down,j1,j2,...
                            -step,expPrm.isFixed,expPrm.sumConstr);
                    otherwise % decreases transition probability
                        mat_down = setProb(mat_down,j1,j2,-step,opt);
                end
                if isempty(mat_down)
                    mat_down = mat_down_prev;
                    gof_down = gof_down_prev;
                    break
                end
                [gof_down,res] = ...
                    calcGOF(degen,mat_down,r,expPrm,opt,h_fig);
                switch opt
                    case 'w'
                        mat_down = res.w_exp;
                    case 'tp'
                        mat_down = res.tp_exp;
                    case 'n'
                        mat_down = res.n_exp;
                end
                if gof_down<=gof_down_prev
                    mat_down = mat_down_prev;
                    gof_down = gof_down_prev;
                end
            end
            alliter = cat(2,alliter,mat_down);
            allgof = cat(2,allgof,gof_down);
        end
    end
end

[bestgof,bestiter] = max(allgof);
mat = alliter{bestiter};


function [gof,res] = calcGOF(degen,mat,r,expPrm,opt,varargin)

% default
fntsz = 10;
plotIt = false;

if ~isempty(varargin)
    h_fig = varargin{1};
    if ~isempty(h_fig)
        plotIt = true;
    end
end
Ls = expPrm.Ls;
expT = expPrm.expT;
dt_exp = expPrm.dt;
fitPrm = expPrm.fitPrm;
excl = expPrm.excl;
J = numel(degen);

% simulate dwell times with degenerated states
switch opt
    case 'n'
        w = mat./repmat(sum(mat,2),[1,size(mat,2)]);
        tp = w.*repmat(expT*r',[1,size(mat,2)]);
        tp(~~eye(size(tp))) = 1-sum(tp,2);
        res = simStateSequences(tp,expT*r,Ls);
    case 'w'
        tp = mat.*repmat(expT*r',[1,size(mat,2)]);
        tp(~~eye(size(tp))) = 1-sum(tp,2);
        res = simStateSequences(tp,expT*r,Ls);
    case 'tp'
        tp = mat;
        res = simStateSequences(tp,expT*r,Ls);
end
dt_sim0 = res.dt;

% fusion dwells from degenerated states
for j = 1:J
    for j_degen = degen{j}'
        dt_sim0(dt_sim0(:,3)==j_degen,3) = j;
        dt_sim0(dt_sim0(:,4)==j_degen,4) = j;
    end
end
mols = unique(dt_sim0(:,2));
dt_sim = [];
for m = mols'
    dt_sim_m = adjustDt(dt_sim0(dt_sim0(:,2)==m,:));
    if excl
        dt_sim_m([1,end],:) = [];
    end
    dt_sim = cat(1,dt_sim,dt_sim_m);
end

% convert frame count to time
dt_sim(:,1) = dt_sim(:,1)*expT;

% save dwell time set
res.dt = dt_sim;

% remove molecule index
dt_sim = dt_sim(:,[1,3:end]);

gof = 0;
sumExp = 0;
sumSim = 0;
cum_counts_exp = cell(1,J);
cum_counts_sim = cell(1,J);
bins = cell(1,J);
ndtExp = cell(1,J);
ndtSim = cell(1,J);
maxDt = zeros(1,J);
minDt = zeros(1,J);
for j1 = 1:J
    dt_sim_j1 = dt_sim(dt_sim(:,2)==j1,1);
    dt_exp_j1 = dt_exp(dt_exp(:,2)==j1,1);

    maxDt(j1) = max([dt_sim_j1;dt_exp_j1]);
    minDt(j1) = min([dt_sim_j1(dt_sim_j1>0);dt_exp_j1(dt_exp_j1>0)]);
    bins{j1} = 0:expT:maxDt(j1);

    cum_counts_exp{j1} = buildDtHist(dt_exp_j1,bins{j1});
    cum_counts_sim{j1} = buildDtHist(dt_sim_j1,bins{j1});
    
    ndtExp{j1} = max(cum_counts_exp{j1});
    ndtSim{j1} = max(cum_counts_sim{j1});
    
    sumExp = sumExp+ndtExp{j1};
    sumSim = sumSim+ndtSim{j1};
end

for j1 = 1:J
    cum_exp = cum_counts_exp{j1}/sumExp;
    maxCumExp = max(cum_counts_exp{j1});
    cum_sim = cum_counts_sim{j1}/sumSim;
    maxCumSim = max(cum_counts_sim{j1});
    cmpl_exp = 1-cum_counts_exp{j1}/maxCumExp;
    cmpl_sim = 1-cum_counts_sim{j1}/maxCumSim;
    
    n = numel(degen{j1});
    
    % fit simulated dwell time histogram (time consuming)
%     fitbounds.start = fitPrm{j1};
%     fitbounds.lower = zeros(1,2*n);
%     fitbounds.upper = Inf(1,2*n);
%     [tau,amp,~] = fitMultiExp(cmpl_sim',bins{j1},fitbounds);
%     
%     fit_cum_sim = zeros(size(bins{j1}));
%     for i = 1:n
%         fit_cum_sim = fit_cum_sim + amp(i)*exp(-bins{j1}/tau(i));
%     end
%     fit_cum_sim = maxCumSim*(1-fit_cum_sim)/sumSim;
%     fit_cmpl_sim = 1-sumSim*fit_cum_sim/maxCumSim;

    % build fit function for experimental dwell time histogram
    fit_cum_exp = zeros(size(bins{j1}));
    for i = 1:n
        fit_cum_exp = fit_cum_exp + ...
            fitPrm{j1}(2*(i-1)+1)*exp(-bins{j1}/fitPrm{j1}(2*i));
    end
    fit_cum_exp = maxCumExp*(1-fit_cum_exp)/sumExp;
    fit_cmpl_exp = 1-sumExp*fit_cum_exp/maxCumExp;
    
    % plot histograms and fit functions
    if plotIt
        lim_y = [min([cum_exp(cum_exp>0)',...
            cum_sim(cum_sim>0)']),...
            max([cum_exp(cum_exp>0)',...
            cum_sim(cum_sim>0)'])];
        dat_id = cmpl_exp>0;
        max_x = find(dat_id==0);
        max_x = max_x(1);
        lim_x = [bins{j1}(1),bins{j1}(max_x)];
        
        h_a = subplot(2,J,j1,'replace','parent',h_fig);
        plot(h_a,bins{j1},cum_exp,'+b');
        hold(h_a,'on');
        plot(h_a,bins{j1},cum_sim,'+r');
        plot(h_a,bins{j1},fit_cum_exp,'-k');
%         plot(h_a,bins{j1},fit_cum_sim,'-k'); % fit of simulation
        hold(h_a,'off');
        ylim(h_a,lim_y);
        xlim(h_a,lim_x);
        
        if j1==1
            switch opt
                case 'n'
                    k = res.n_exp./...
                        repmat(sum(res.n_exp,2)./r',1,size(res.n_exp,2));
                case 'w'
                    k = repmat(r',1,size(res.n_exp,2)).*res.w_exp;
                case 'tp'
                    k = res.tp_exp;
                    k(~~eye(size(k))) = 0;
                    k = k/expT;
            end
            str = repmat('%0.3f  \t',[1,size(k,2)]);
            str(end) = 'n';
            str = sprintf(str,k');
            t = text(h_a,0,0,str);
            ext = get(t,'extent');
            set(t,'position',[lim_x(2)-ext(3),lim_y(1)+ext(4)],'fontunits',...
                'points','fontsize',fntsz)
        end
        
        lim_y = [min([cmpl_exp(cmpl_exp>0)',...
            cmpl_sim(cmpl_sim>0)']),...
            max([cmpl_exp(cmpl_exp>0)',...
            cmpl_sim(cmpl_sim>0)'])];

        h_a = subplot(2,J,J+j1,'replace','parent',h_fig);
        plot(h_a,bins{j1},cmpl_exp,'+b');
        hold(h_a,'on');
        plot(h_a,bins{j1},cmpl_sim,'+r');
        plot(h_a,bins{j1},fit_cmpl_exp,'-k');
%         plot(h_a,bins{j1},fit_cmpl_sim,'-k'); % fit of simulation
        hold(h_a,'off');
        set(h_a,'yscale','log');
        ylim(h_a,lim_y);
        xlim(h_a,lim_x);

        drawnow;
    end
    
    % calculate log probability
    incl = cmpl_sim>0 & fit_cmpl_exp'>0;
    Y_cmpl = log(cmpl_sim(incl));
    Yfit_cmpl = log(fit_cmpl_exp(incl)');
    
    incl = cum_sim>0 & fit_cum_exp'>0;
    Y_cum = log(cum_sim(incl));
    Yfit_cum = log(fit_cum_exp(incl)');
    
    % calculate GOF = sum(P)/RSSE
    gof_lin_cum = sum(cum_sim)/sqrt(sum((fit_cum_exp'-cum_sim).^2));
    gof_lin_cmpl = sum(cmpl_sim)/sqrt(sum((fit_cmpl_exp'-cmpl_sim).^2));
    gof_log_cmpl = -sum(Y_cmpl)/sqrt(sum(((Yfit_cmpl-Y_cmpl).^2)));
    
    % increment GOF
    gof = gof + log(gof_lin_cum);
end


function dtHist = buildDtHist(dt,bins)

dtHist = hist(dt,bins)';
dtHist = cumsum(dtHist);
% dtHist = 1-dtHist/dtHist(end);


function mat = setProb(mat,j1,j2,step,opt)

% default
pmax = 1;
if j1==j2
    pmin = 0.1;
else
    pmin = 0;
end

J = size(mat,2);

switch opt
    case 'w'
        j2s = find((1:J)~=j2 & (1:J)~=j1);
    case 'tp'
        if j1~=j2
            pmax = 1-mat(j1,j1);
            j2s = find((1:J)~=j2 & (1:J)~=j1);
        else
            j2s = find((1:J)~=j2);
        end
end

% increment current transition probability
if step>0
    if (mat(j1,j2)+step)>=pmax
        mat(j1,j2) = pmax;
    else
        mat(j1,j2) = mat(j1,j2)+step;
    end
else
    if (mat(j1,j2)+step)<=pmin
        mat(j1,j2) = pmin;
    else
        mat(j1,j2) = mat(j1,j2)+step;
    end
end

if sum(mat(j1,j2s))==0
    mat(j1,j2s) = (pmax-mat(j1,j2))/numel(j2s);
else
    mat(j1,j2s) = (pmax-mat(j1,j2))*mat(j1,j2s)/sum(mat(j1,j2s));
end

% for j1 = 1:J
%     for j2 = 1:J
%         if sum(~~mat(j1,:))==1 && sum(~~mat(:,j2))==1
%             mat = [];
%             return
%         end
%     end
% end


