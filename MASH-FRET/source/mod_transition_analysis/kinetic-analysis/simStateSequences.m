function res = simStateSequences(tp,r,Ls,varargin)
% res = simStateSequences(tp,r,Ls)
% res = simStateSequences(tp,r,Ls,ip)
%
% Generate state sequences from probability-based parameters
%
% tp: [J-by-J] transition probability matrix
% r: [1-by-J] restricted rate constants (in frame-1)
% Ls: [1-by-N] trajectory lengths (in frame)
% ip: [1-by-J] initial state probabilities
% res: resturned results as a structure containing fields:
%  res.dt: [nDt-by-4] dwell times (molecule index, state durations (frames), state values before and after transition)
%  res.tp_exp: [J-by-J] transition probability matrix calculated from data
%  res.r_exp: [J-by-J] restricted rate matrix calculated from data
%  res.w_exp: [J-by-J] weighing factors calculated from data
%  res.n_exp: [J-by-J] numbers of transitions calculated from data

% initialize output
res = [];

J = size(tp,2);
N = numel(Ls);

% get pre-defined initial state probabilities
ip = [];
for arg = 1:numel(varargin)
    if numel(varargin{arg})==J
        ip = varargin{arg};
    end
end

% calculate initial state probabilities
if isempty(r)
    r = -log(tp(~~eye(size(tp))));
end
r(r==0) = 1E-6;
if isempty(ip)
    ip = zeros(size(r));
    ip(~isinf(r)) = (1./r(~isinf(r)))./sum(1./r(~isinf(r)));
    ip(isnan(ip)) = Inf;
end

% calculate transition probabilities
w = zeros(size(tp));
w(~eye(size(tp))) = tp(~eye(size(tp)));
w = w./repmat(sum(w,2),[1,J]);
w(isnan(w)) = 0;
w(:,ip==0) = 0;

dt = [];
res.traces = {};
for n = 1:N
    L = Ls(n);
    state1 = randsample(1:J,1,true,ip);
%     seq = NaN(L,1);
%     seq(1,1) = state1;
%     for l = 2:L
%         state2 = randsample(1:J,1,true,tp(state1,:));
%         seq(l,1) = state2;
%         state1 = state2;
%     end
%     dt_n = getDtFromDiscr(seq,1);
%     dt = cat(1,dt,[dt_n(:,1),repmat(n,[size(dt_n,1),1]),dt_n(:,2:end)]);
    
    dt_n = [];
    l = 0;
    while l<L
        if ~isinf(r(state1))
            dl = random('exp',1/(r(state1)));
            if (l+dl)>L
                dl = L-l;
            end
            l = l+dl;
        end
        if sum(w(state1,:))==0
            state2 = state1;
        else
            state2 = randsample(1:J,1,true,w(state1,:));
        end
        dt_n = cat(1,dt_n,[dl,state1,state2]);
        state1 = state2;
    end
    
    dt = cat(1,dt,[dt_n(:,1),repmat(n,[size(dt_n,1),1]),dt_n(:,2:end)]);
end
dt(dt(:,1)<1,:) = [];

w_exp = zeros(J);
tp_exp = zeros(J);
n_exp = zeros(J);
for j1 = 1:J
    dt_j1 = dt(dt(:,3)==j1,:);
    for j2 = 1:J
        if j1==j2
            continue
        end
        dt_j1j2 = dt_j1(dt_j1(:,4)==j2,1);
        w_exp(j1,j2) = numel(dt_j1j2)/size(dt_j1,1);
        tp_exp(j1,j2) = numel(dt_j1j2)/sum(dt_j1(:,1));
        n_exp(j1,j2) = numel(dt_j1j2);
    end
end
tp_exp(~~eye(size(tp_exp))) = 1-sum(tp_exp,2);

res.dt = dt;
res.tp_exp = tp_exp;
res.w_exp = w_exp;
res.n_exp = n_exp;

