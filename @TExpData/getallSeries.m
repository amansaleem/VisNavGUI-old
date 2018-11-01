function iseries_list = getallSeries(obj)
    infoAll = getDataInfo(obj.animal);
    date_list = [];
    for n = 1:length(infoAll)
        date_list = [date_list, {infoAll(n).date}];
    end
    
    if ismac
        date_list{1} = 'bloddy mac';
    end
    
    iseries_list = zeros(1,numel(date_list));
    for n = 1:numel(date_list)
        iseries_list(n) = str2num(date_list{n});
    end    
end