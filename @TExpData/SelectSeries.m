function obj = SelectSeries(obj,selmode)
if nargin < 2
    selmode = 'single';
end
infoAll = getDataInfo(obj.animal);
date_list = [];
for n = 1:length(infoAll)
    date_list = [date_list, {infoAll(n).date}];
end

if ismac
    date_list{1} = 'bloddy mac';
end
[Selection,ok] = listdlg('PromptString','Select a series:',...
    'SelectionMode',selmode,...
    'ListString',date_list);
if ok
    obj.series = [];
    obj.iseries = [];
    date = date_list(Selection);
    for i = 1:numel(date)
        obj.series = [obj.series date(i)];
        obj.iseries = [obj.iseries str2num(date{i})];
    end
end
end