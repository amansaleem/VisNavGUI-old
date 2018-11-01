function obj = SelectAnimal(obj, DIRS, strcell)
if nargin < 3
    strcell = [];
end
dir_list = dir(DIRS.ball);
animal_list = [];
for n = 3:length(dir_list)
    mdate(n-2) = dir_list(n).datenum;
end
[~, dorder] = sort(mdate);
for n = length(dorder):-1:1
    if ~isempty(strcell)
        for s = 1:numel(strcell)
            pos = strfind(dir_list(dorder(n)+2).name,strcell{s});
            if ~isempty(pos)
                animal_list = [animal_list, {dir_list(dorder(n)+2).name}];
            end
        end
    else
        animal_list = [animal_list, {dir_list(dorder(n)+2).name}];
    end
end
[Selection,ok] = listdlg('PromptString','Select an animal:',...
    'SelectionMode','single',...
    'ListString',animal_list);
if ok
    obj.animal = animal_list(Selection);
    obj.animal = obj.animal{1};
end
end