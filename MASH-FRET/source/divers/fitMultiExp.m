function [tau,amp,gof] = fitMultiExp(y,x,varargin)
% [tau,amp,gof] = fitMultiExp(y,x)
% [tau,amp,gof] = fitMultiExp(y,x,p)
%
% x: x-data
% y: y-data
% p: structure containing fields:
%  p.upper: [1-by-2*n] upper bounds for fitting parameters
%  p.lower: [1-by-2*n] lower counds for fittinf parameters
% tau: [n-by-1] decay contants of the n expoenntial functions
% amp: [n-by-1] relative amplitudes of the n exponential functions
% gof: goodness of fit

% defaults
n = 0;
plotIt = false;
fmt_hist = '-r';
fmt_fit = '-b';
ampmin = 0;
taumin = 0;
ampmax = 2;
taumax = Inf;

if ~isempty(varargin)
    fitbounds = varargin{1};
    n = size(fitbounds.upper,2)/2;
end

Y = log(y);
excl = isinf(Y);
Y(excl) = [];
x(excl) = [];

if plotIt
    h_fig = figure();
    h_axes = axes('parent',h_fig);
else
    h_axes = [];
end

% find inflexion points
cp = findChangePoint(Y,x,h_axes,fmt_hist,fmt_fit);

datrange = 1:numel(x);
if n==0
    cp = cp(1,:);
    n = size(cp,2)+1;
    minVal = repmat([ampmin;taumin],[1,n]);
    maxVal = repmat([ampmax;taumax],[1,n]);
    
else
    if ~isempty(cp)
        if size(cp,2)<(n-1)
            cp = cat(2,cp,repmat(cp(:,end),[1,(n-1)-size(cp,2)]));
        elseif size(cp,2)>(n-1)
            [~,id] = sort(cp(2,:),'descend');
            if n==1
                datrange = 1:cp(1,1);
                cp = [];
            else
                 cp = cp(:,id(1:n-1));
                [~,id] = sort(cp(1,:));
                cp = cp(:,id);
            end
        end
        if ~isempty(cp)
            cp = cp(1,:);
        end
    else
        cp = (1:(n-1))*numel(x)/n;
    end
    minVal = reshape(fitbounds.lower,[2,n]);
    maxVal = reshape(fitbounds.upper,[2,n]);
end

[tau0,amp0] = estimateStartGuess(x(datrange),Y(datrange),cp,n);

[prm,gof] = optimizeMultiExpModel(x(datrange),Y(datrange),[amp0;tau0],...
    minVal,maxVal,h_axes,fmt_hist,fmt_fit);
amp = prm(1,:)/sum(prm(1,:));
tau = prm(2,:);

if plotIt
    hold(h_axes,'off');
    plot(h_axes,x,Y,fmt_hist);
    hold(h_axes,'on');

    % plot fit
    yfit = zeros(size(Y));
    for ni = 1:n
        yfit = yfit + amp(ni)*exp(-x/tau(ni))/sum(amp);
    end
    plot(h_axes,x,log(yfit),fmt_fit);
    close(h_fig);
end


function cp = findChangePoint(Y,x,h_axes,fmt_max,fmt_fit)

% defaults
plotIt = true;

% default
cp = [];
tol = 0.5;

% determine parameters of base-line exponential function
A0 = (Y(end)-Y(1))/(x(end)-x(1));
B0 = Y(end)-A0*x(end);
Y0 = B0 + A0*x;

% find maximum difference between base line and data
dk = Y0-Y;
[dkmax,id] = max(dk);
[dkmin1,id1] = min(dk(1:id));
[dkmin2,id2] = min(dk(id:end));
id2 = id2+id-1;

% calculate data slope at each side of change point
A1 = (dk(id)-dkmin1)/(x(id)-x(id1));
A1mean = mean((dk(id)-dk(1:id-1))./(x(id)-x(1:id-1)));
A1std = std((dk(id)-dk(1:id-1))./(x(id)-x(1:id-1)));
B1 = dk(id)-A1*x(id);
dk1 = B1+A1*x(1:id);

A2 = (dkmin2-dk(id))/(x(id2)-x(id));
A2mean = mean((dk(id+1:end)-dk(id))./(x(id+1:end)-x(id)));
A2std = std((dk(id+1:end)-dk(id))./(x(id+1:end)-x(id)));
B2 = dk(id)-A2*x(id);
dk2 = B2+A2*x(id:end);

% plot 
if plotIt && ~isempty(h_axes)
    plot(h_axes,x,dk,'+');
    y_lim = get(h_axes,'ylim');
    hold(h_axes,'on');
    plot(h_axes,[x(id),x(id)],y_lim,fmt_max);
    set(h_axes,'ylim',[min(dk),max(dk)]);
    plot(h_axes,x(1:id),dk1,fmt_fit);
    plot(h_axes,x(id:end),dk2,fmt_fit);
    hold(h_axes,'off');
end

if A1>0 && A2<0
    cp = [id;dkmax];
    if ~((A1mean-3*A1std)<0 && (A1mean+3*A1std)>0) && ...
            sum(dk(1:id)>dk1)/numel(dk(1:id))>tol 
        cp = cat(2,cp,...
            findChangePoint(Y(1:id),x(1:id),h_axes,fmt_max,fmt_fit));
    end
    if ~((A2mean-3*A2std)<0 && (A2mean+3*A2std)>0) && ...
            sum(dk(id:end)>dk2)/numel(dk(id:end))>tol
        cp = cat(2,cp,...
            findChangePoint(Y(id:end),x(id:end),h_axes,fmt_max,fmt_fit));
    end
end

if ~isempty(cp)
    [~,id] = sort(cp(1,:));
    cp = cp(:,id);
end

if plotIt && ~isempty(h_axes)
    hold(h_axes,'off');
end


function [tau0,amp0] = estimateStartGuess(x,Y,cp,n)

cp = [1,cp,numel(x)];
tau0 = NaN(1,n);
amp0 = NaN(1,n);
for ni = 1:n
    xstart = cp(ni);
    xstop = cp(ni+1);
    Yn = Y(xstart:xstop);
    xn = x(xstart:xstop);
    
    % linear regression
    mdl = fitlm(xn,Yn);
    amp0(ni) = exp(mdl.Coefficients.Estimate(1,1));
    tau0(ni) = -1/mdl.Coefficients.Estimate(2,1);
end


function [prm,gof] = optimizeMultiExpModel(x,Y,prm0,minVal,maxVal,h_axes,...
    fmt_hist,fmt_fit)

% default
plotIt = false;
minVar = 1E-6;
gof_prev = -Inf;
gof = 0;
prm = prm0;
prm_prev = prm0;
n = size(prm,2);
varsteps = [0.01,0.1];
ncycle = 3;

for cycle = 1:ncycle
    while gof>(gof_prev+minVar)

        % plot interation
        if plotIt && ~isempty(h_axes)
            plot(h_axes,x,Y,fmt_hist);
            hold(h_axes,'on');
            yfit = zeros(size(Y));
            for ni = 1:n
                yfit = yfit + prm(1,ni)*exp(-x/prm(2,ni))/sum(prm(1,:));
            end
            plot(h_axes,x,log(yfit),fmt_fit);
            hold(h_axes,'off');
            drawnow;
        end

        gof_prev = gof;
        prm_prev = prm;

        % move each function after the other
        for ni = 1:n
            prm = varyParam(x,Y,prm,ni,varsteps(cycle,:),maxVal,minVal,...
                minVar);
        end

        gof = calcGOF(prm,x,Y);
    end
    prm = prm_prev;
    gof = gof_prev;
end


function prm = varyParam(x,Y,prm,n,step,prm_max,prm_min,minVar)

gof = calcGOF(prm,x,Y);

nPrm = size(prm,1);

for id = 1:nPrm
    % up iteration
    prm_up = prm;
    prm_up_prev = prm;
    gof_up_prev = -Inf;
    gof_up = gof;
    while gof_up>(gof_up_prev+minVar) && ...
            (prm_up(id,n)+step(id))<prm_max(id,n)
        gof_up_prev = gof_up;
        prm_up_prev = prm_up;

        prm_up(id,n) = prm_up(id,n)+step(id);

        gof_up = calcGOF(prm_up,x,Y);
    end

    % down iteration
    prm_down = prm;
    prm_down_prev = prm;
    gof_down_prev = -Inf;
    gof_down = gof;
    while gof_down>(gof_down_prev+minVar) && ...
            (prm_down(id,n)-step(id))>prm_min(id,n)
        gof_down_prev = gof_down;
        prm_down_prev = prm_down;

        prm_down(id,n) = prm_down(id,n)-step(id);

        gof_down = calcGOF(prm_down,x,Y);
    end
    
    alliter = {prm,prm_up_prev,prm_down_prev};
    [~,j] = max([gof,gof_up,gof_down]);
    prm = alliter{j};
end


function gof = calcGOF(prm,x,Y)

a = prm(1,:);
b = prm(2,:);
N = numel(Y);
yfit = zeros(size(Y));
n = numel(a);
for ni = 1:n
    yfit = yfit + a(ni)*exp(-x/b(ni))/sum(a);
end
w = Y/sum(Y);

gof = N/sqrt(sum(((log(yfit)-Y).^2).*w));


