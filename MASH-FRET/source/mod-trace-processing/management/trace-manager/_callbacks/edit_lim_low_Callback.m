function edit_lim_low_Callback(obj,evd,h_fig,row)

% Last update: by RB, 4.1.2018
% >> adapted for FRET-S-histograms
% >> hist2 rather slow replaced by hist2D
%
%

h = guidata(h_fig);
p = h.param.ttPr;
proj = p.curr_proj;
nChan = p.proj{proj}.nb_channel;
nExc = p.proj{proj}.nb_excitations;
FRET = p.proj{proj}.FRET;
nFRET = size(FRET,1);
S = p.proj{proj}.S;
nS = size(S,1);

dat1 = get(h.tm.axes_ovrAll_1, 'UserData');
dat2 = get(h.tm.axes_ovrAll_2, 'UserData');
plot2 = get(h.tm.popupmenu_axes2, 'Value');
lim_low = str2num(get(obj,'String'));

if lim_low >= dat1.lim{plot2}(row,2)
    setContPan('Lower bound must be lower than higher bound.','error',...
        h_fig);
    return;
end

dat1.lim{plot2}(row,1) = lim_low;

if plot2 <= nChan*nExc+nFRET+nS
    [dat2.hist{plot2},dat2.iv{plot2}] = getHistTM(dat1.trace{plot2},...
        dat1.lim{plot2},dat1.niv(plot2,1));
    
else
    ind = plot2-nChan*nExc-nFRET-nS;
    ind_e = ceil(ind/nS) + nChan*nExc;
    ind_s = ind - (ceil(ind/nS)-1)*nS + nChan*nExc + nFRET;
    ES = [dat1.trace{ind_e},dat1.trace{ind_s}];
    [dat2.hist{plot2},dat2.iv{plot2}] = getHistTM(ES,...
        dat1.lim{plot2},dat1.niv(plot2,[1,2]));
end

set(h.tm.axes_ovrAll_1, 'UserData', dat1);
set(h.tm.axes_ovrAll_2, 'UserData', dat2);
plotData_overall(h_fig);
    