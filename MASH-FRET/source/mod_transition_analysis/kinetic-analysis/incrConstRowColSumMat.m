function mat = incrConstRowColSumMat(mat0,j1,j2,incr,isFixed,sumConstr)
% mat = incrConstRowColSumMat(mat,j1,j2,incr,isFixed,sumConstr)
%
% Increment matrix element located at (j1,j2) and recalculate other elements to respect the constant cum of elements in each row and column.
%
% mat: [J-by-J] input matrix
% j1: row index in matrix of element to increment
% j2: columns index in matrix of element to increment
% incr: increment
% isFixed: [J-by-J] matrix indicating if an element is fixed (1) or can be varied (0)
% sumConstr: {J-by-J-by} cell matrix containg:
%   sumCnstr{j1,j2}: [J-by-J] matrix containing in cell (a1,a2) the sum of all elements to which the element (j1,j2) is constrainted to or (0) if not constrainted
% mat: re-calculated matrix

% initialize output
mat = mat0;

% defaults
Nmin = 0;

if isFixed(j1,j2)
    return
end

% get state indexes other than element
J = size(mat0,2);
js = (1:J);

% increment element
% | 0 0 x 0 0 0 |
% | 0 0 0 0 0 0 |
% | 0 0 0 0 0 0 |
% | 0 0 0 0 0 0 |
% | 0 0 0 0 0 0 |
% | 0 0 0 0 0 0 |
Nmax = min(min(sumConstr{j1,j2}(sumConstr{j1,j2}>0)));
if (mat(j1,j2)==Nmax && incr>0) || (mat(j1,j2)==Nmin && incr<0)
    return
end
if mat(j1,j2)<Nmax && (mat(j1,j2)+incr)>Nmax
    incr = Nmax-mat(j1,j2);
elseif mat(j1,j2)>Nmin && (mat(j1,j2)+incr)<Nmin
    incr = -(mat(j1,j2)-Nmin);
end
mat(j1,j2) = mat(j1,j2)+incr;

% varies an element contrained by sums
if ~isinf(Nmax)
    % re-calculate elements constraint by the particular row-sum
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | x 1 1 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    cnstr = ~~sumConstr{j1,j2}(j1,:);
    rowSum = unique(sumConstr{j1,j2}(j1,cnstr));
    if sum(mat(j1,cnstr))==0
        mat(j1,cnstr) = (rowSum-mat(j1,j2))/sum(cnstr);
    else
        mat(j1,cnstr) = (rowSum-mat(j1,j2)).*mat(j1,cnstr)/sum(mat(j1,cnstr));
    end
    
    % re-calculate elements constraint by the particular column-sum
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | x 0 0 0 0 0 |
    % | 1 0 0 0 0 0 |
    % | 1 0 0 0 0 0 |
    cnstr = ~~sumConstr{j1,j2}(:,j2);
    colSum = unique(sumConstr{j1,j2}(cnstr,j2));
    if sum(mat(cnstr,j2))==0
        mat(cnstr,j2) = (colSum-mat(j1,j2))/sum(cnstr);
    else
        mat(cnstr,j2) = ...
            (colSum-mat(j1,j2)).*mat(cnstr,j2)./sum(mat(cnstr,j2));
    end
    
    % re-calculate elements constrainted by the particular row-sum of each element previously varied
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | x 0 0 0 0 0 |
    % | 0 1 1 0 0 0 |
    % | 0 1 1 0 0 0 |
    j2s = js(~~sumConstr{j1,j2}(:,j2));
    n = 0;
    for j = j2s
        n = n+1;
        cnstr = ~~sumConstr{j,j2}(j,:);
        rowSum = unique(sumConstr{j,j2}(j,cnstr));
        if sum(mat(j,cnstr))==0
            mat(j,cnstr) = (rowSum-mat(j,j2))/sum(cnstr);
        else
            mat(j,cnstr) =...
                (rowSum-mat(j,j2)).*mat(j,cnstr)./sum(mat(j,cnstr));
        end
    end
    
    % re-calculate elements constrainted by the particular column-sum of each element previously varied
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | x 1 1 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    j2s = js(~~sumConstr{j1,j2}(j1,:));
    for j = j2s
        cnstr = ~~sumConstr{j1,j}(:,j);
        colSum = unique(sumConstr{j1,j}(cnstr,j));
        mat(j1,j) = colSum-sum(mat(cnstr,j));
    end

% varies an element not constrained by sum (row- and column-sums are re-defined after increment)
else
    % incr. elements constrained by row- and column-sum
    % | 0 0 x 0 0 0 | > N
    % | 1 0 0 0 0 0 |
    % | 1 1 0 1 0 0 |
    % | 1 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    %   v
    %   N
    isFixed = ~sumConstr{j1,j2};
    isFixed(j1,:) = true;
    isFixed(:,j2) = true;
    Nsum = zeros(1,size(mat,2));
    Nsum(j1) = sum(mat(j1,:));
    Nsum(j2) = sum(mat(:,j2));
    
    row = j2;
    col = j1;
    mat = reorgConstRowColSumMat(mat,isFixed,row,col,Nsum);
    
    % re-calculate all other elements 
    % | 0 0 x 0 0 0 |
    % | 0 0 0 1 0 0 | > N1
    % | 0 0 0 0 0 0 |
    % | 0 1 0 0 0 0 | > N2
    % | 0 0 0 0 0 0 |
    % | 0 0 0 0 0 0 |
    %     v   v
    %     N1  N2
    isFixed(row,:) = true;
    isFixed(:,col) = true;
    js = find(sum(sumConstr{j1,j2},1));
    js = js(js~=j1 & js~=j2);
    Nsum0 = sum(mat0,1);
    Nsum(js) = sum(sum(mat(js,:)))*Nsum0(js)/sum(Nsum0(js));
    
    for row = js
        for col = js
            if row==col
                continue
            end
            mat = reorgConstRowColSumMat(mat,isFixed,row,col,Nsum);
        end
    end
end

mat(mat<0) = 0;

if ~isempty(mat)
    for j1 = 1:J
        if sum(mat(j1,:)>0)==1 && sum(mat(mat(j1,:)>0,:)>0)==1 && ...
                find(mat(mat(j1,:)>0,:)>0)==j1
            mat = [];
            return
        end
    end
end

