function pushbutton_TDPremProj_Callback(obj, evd, h)
p = h.param.TDP;
if ~isempty(p.proj)
    
    % collect selected projects
    slct = get(h.listbox_TDPprojList, 'Value');
    
    % build confirmation message box
    str_proj = ['"' p.proj{slct(1)}.exp_parameters{1,2} '"'];
    for pj = 2:numel(slct)
        str_proj = [str_proj ', "' p.proj{slct(pj)}.exp_parameters{1,2} ...
            '"'];
    end
    del = questdlg(['Remove project ' str_proj ' from the list?'], ...
        'Remove project', 'Yes', 'No', 'No');
    
    if strcmp(del, 'Yes')
        
        % build action
        list_str = get(h.listbox_TDPprojList,'String');
        str_act = '';
        for i = slct
            str_act = cat(2,str_act,'"',list_str{i},'" (',...
                p.proj{i}.proj_file,')\n');
        end
        str_act = str_act(1:end-2);
        
        % delete projects and reorganize project and current data 
        % structures
        projLst = {};
        curr_type = [];
        for i = 1:size(p.proj,2)
            if prod(double(i ~= slct))
                projLst{size(projLst,2)+1} = p.proj{i};
                curr_type(size(curr_type,2)+1) = p.curr_type(i);
            end
        end
        p.proj = projLst;
        p.curr_type = curr_type;
        
        % set new current project
        if size(projLst,2) <= 1
            p.curr_proj = 1;
        elseif slct(end) < size(p.proj,2)
            p.curr_proj = slct(end)-numel(slct) + 1;
        else
            p.curr_proj = slct(end)-numel(slct);
        end
        
        % update project list
        p = ud_projLst(p, h.listbox_TDPprojList);
        h.param.TDP = p;
        guidata(h.figure_MASH, h);
        
        % clear axes
        cla(h.axes_TDPplot1);
        
        % update calculations and GUI
        updateFields(h.figure_MASH, 'TDP');
    end
end