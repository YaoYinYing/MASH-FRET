function [dt,w_vect,mol_incl,ok,errmsg] = getDtFromState(dat,state,excl)
% [dt,w_vect,mol_incl,ok,errmsg] = getDtFromState(dat,state,excl)
%
% Gather dwell times of a same state and exclude/include first and last
% dwell times of each sequences.
%
% dt: [nDt-by-8] dwell times, state value before and after transition, molecule index, x- and y- TDP cooridnates, state indexes before and after transition
% sate: state index in list
% excl: (1) to exclude first & last dwell times (0) otherwise
% dt: {1-by-N} dwell times sorted by molecule
% w_vect: [1-by-N] molecule weight
% mol_incl: [1-by-N] true if molecule shows the state, false otherwise
% ok: computation success (1) or failure (0)
% errmsg: error message in case of failure

% Created, 25.4.2020 by MH.

ok = 0;
errmsg = '';
dt = {};
mol_incl = [];
w_vect = [];

% get dwell times concerning state only
dat_state = dat(dat(:,7)==state,1:end-2);

mols = unique(dat_state(:,4));
N = numel(mols);

for n = 1:N
    dt_m = dat(dat(:,4)==mols(n),:);
    
    % exclude molecules without dwell times concerning state->... transitions
    if isempty(dt_m)
        disp(cat(2,'molecule ',num2str(mols(n)),' excluded: no dwell time'));
        continue
    end
    
    if excl && isempty(dt_m(2:end-1,:))
        disp(cat(2,'molecule ',num2str(mols(n)),' excluded: no dwell time',...
            ' left after exclusion.'));
        continue
        
    elseif excl
        dt_m = dt_m(2:end-1,:);
    end

    dt_m_state = dt_m(dt_m(:,7)==state,:);
    
    if isempty(dt_m_state)
        disp(cat(2,'molecule ',num2str(mols(n)),' excluded: no dwell time',...
            ' left after exclusion.'));
        continue
    end
    
    dt{size(dt,2)+1} = dt_m_state;
    mol_incl = cat(2,mol_incl,mols(n));

    w_vect = cat(1,w_vect,sum(dt{end}(:,1)));

    if n == N
        w_vect = w_vect/sum(w_vect);
    end
end

if size(dt,2)==0
    errmsg = cat(2,'No dwell time is left after excluding first last ',...
        'occurrence in each trajectory.');
    return
end

ok = 1;
