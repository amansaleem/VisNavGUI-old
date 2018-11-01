function obj = SelectExpt(obj, Fall)
if nargin < 2
    Fall = false;
end
infoAll = getDataInfo(obj.animal);
date_list = [];
for n = 1:length(infoAll)
    date_list = [date_list, {infoAll(n).date}];
end
whichdate = obj.series;
for i=1:numel(infoAll)
    if strcmpi(infoAll(i).date,whichdate)==1
        expt_info = infoAll(i);
        break;
    end
end
for iexp = 1:length(expt_info.sessions)
    expt_list{iexp} = num2str(expt_info.sessions(iexp));
end
if ~Fall
    [Selection,ok] = listdlg('PromptString','Select experiment(s):',...
        'SelectionMode','multiple',...
        'ListString',expt_list);
else
    ok = true;
    Selection = 1:numel(expt_info.sessions);
end
if ok
    obj.exptList = expt_info.sessions(Selection);
end
end