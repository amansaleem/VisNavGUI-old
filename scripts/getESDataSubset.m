function es_out = getESDataSubset(es_in,Name,minValue,maxValue)
if nargin < 2
    Name = '';
    Fshowdialog = true;
else
    Fshowdialog = false;
end

clear global field_vals es_g
global field_vals
global es_g

es_g = es_in;

field_list = fieldnames(es_g);
processComplete = 0;
icount = 1;
for i = 1:length(field_list)
    if isnumeric(eval(['es_g.' field_list{i}])) && ~strcmp(field_list{i}, 'sampleRate')
        field_vals(icount).minVal = nanmin(eval(['es_g.' field_list{i} '(:)']));
        field_vals(icount).maxVal = nanmax(eval(['es_g.' field_list{i} '(:)']));
        field_vals(icount).name = field_list{i};
        if strcmp(Name,field_list{i})
            if ~isempty(minValue)
                field_vals(icount).minVal = minValue;
            end
            if ~isempty(maxValue)
                field_vals(icount).maxVal = minValue;
            end            
        end            
        icount = icount + 1;
    end
end
if Fshowdialog
    figNum = createUI(es_g, field_list);
    uiwait(figNum.figInfo.fig);

    es_out = es_g;
    close(figNum.figInfo.fig)
else
    cutProcess0;
    es_out = es_g;
end

clear global field_vals es_g
    function cutProcess0
        for ifield = 1:length(field_vals)
            temp = eval(['es_g.' field_vals(ifield).name]);
            temp(temp<field_vals(ifield).minVal) = NaN;
            temp(temp>field_vals(ifield).maxVal) = NaN;
            eval(['es_g.' field_vals(ifield).name '= temp;']);
        end
        processComplete = 1;
    end
        
    function cutProcess(ESsubsetGUI)
        for ifield = 1:length(field_vals)
            temp = eval(['es_g.' field_vals(ifield).name]);
            temp(temp<field_vals(ifield).minVal) = NaN;
            temp(temp>field_vals(ifield).maxVal) = NaN;
            eval(['es_g.' field_vals(ifield).name '= temp;']);
        end
        processComplete = 1;
        uiresume(ESsubsetGUI.figInfo.fig)
%         es_out = es_g;
%         close(ESsubsetGUI.figInfo.fig)
        
%         assignin('base','es_new',es2)
    end
    function clipProcess(ESsubsetGUI)
        for ifield = 1:length(field_vals)
            temp = eval(['es_g.' field_vals(ifield).name]);
            temp(temp<field_vals(ifield).minVal) = field_vals(ifield).minVal;
            temp(temp>field_vals(ifield).maxVal) = field_vals(ifield).maxVal;
            eval(['es_g.' field_vals(ifield).name '= temp;']);
        end
        processComplete = 1;
        uiresume(ESsubsetGUI.figInfo.fig)
%         es_out = es_g;
%         close(ESsubsetGUI.figInfo.fig)
%         assignin('base','es_g',es_g)
    end
    function setMax(idx)
        field_vals(idx).maxVal = str2num(get(field_vals(idx).maxFig,'String'));
    end
    function setMin(idx)
        field_vals(idx).minVal = str2num(get(field_vals(idx).minFig,'String'));
    end
    function ESsubsetGUI = createUI(es_g, field_list)
        ESsubsetGUI.figInfo.fig = figure('Name', 'Data subset',...
            'Position', [200 50 500 1000], ...
            'MenuBar', 'none', ...
            'ToolBar', 'none', ...
            'NumberTitle', 'off');
        
        ESsubsetGUI.figInfo.layout = uiextras.HBox('Parent',ESsubsetGUI.figInfo.fig,'Spacing',10,'Padding',5);
        ESsubsetGUI.figInfo.textCol = uiextras.VBox('Parent',ESsubsetGUI.figInfo.layout,'Spacing',5);
        ESsubsetGUI.figInfo.minCol = uiextras.VBox('Parent',ESsubsetGUI.figInfo.layout,'Spacing',5);
        ESsubsetGUI.figInfo.maxCol = uiextras.VBox('Parent',ESsubsetGUI.figInfo.layout,'Spacing',5);
        ESsubsetGUI.figInfo.prosCol = uiextras.VBox('Parent',ESsubsetGUI.figInfo.layout,'Spacing',5);
        uicontrol('Parent',ESsubsetGUI.figInfo.prosCol, 'String', 'CUT','fontsize',14,'HorizontalAlignment','center', 'Callback', @(src,evt) cutProcess(ESsubsetGUI));
        uicontrol('Parent',ESsubsetGUI.figInfo.prosCol, 'String', 'CLIP','fontsize',14,'HorizontalAlignment','center', 'Callback', @(src,evt) clipProcess(ESsubsetGUI));
        idx = 1;
        
        uicontrol('Style','text','Parent',ESsubsetGUI.figInfo.textCol,'String','FIELD NAME', 'fontsize',10, 'HorizontalAlignment','center');
        uicontrol('Style','text','Parent',ESsubsetGUI.figInfo.minCol,'String','MIN VALUE', 'fontsize',10, 'HorizontalAlignment','center');
        uicontrol('Style','text','Parent',ESsubsetGUI.figInfo.maxCol,'String','MAX VALUE', 'fontsize',10, 'HorizontalAlignment','center');
        
        for ifield = 1:length(field_list)
            if isnumeric(eval(['es_g.' field_list{ifield}])) & ~strcmp(field_list{ifield}, 'sampleRate')
%                 field_vals(idx).minVal = nanmin(eval(['es_g.' field_list{ifield} '(:)']));
%                 field_vals(idx).maxVal = nanmax(eval(['es_g.' field_list{ifield} '(:)']));
%                 field_vals(idx).name = field_list{ifield};
                
                % Min value
                field_vals(idx).minFig = uicontrol('Parent', ESsubsetGUI.figInfo.minCol, 'Style','edit', 'Callback',@(src,evt) setMin(idx));
                set(field_vals(idx).minFig, 'String',num2str(field_vals(idx).minVal))
                % Max value
                field_vals(idx).maxFig = uicontrol('Parent', ESsubsetGUI.figInfo.maxCol, 'Style','edit', 'Callback',@(src,evt) setMax(idx));
                set(field_vals(idx).maxFig, 'String',num2str(field_vals(idx).maxVal));
                uicontrol('Style','text','Parent',ESsubsetGUI.figInfo.textCol,'String',[field_vals(idx).name ':'], 'fontsize',10, 'HorizontalAlignment','right');
                
                idx = idx + 1;
            end
        end
        
    end
end