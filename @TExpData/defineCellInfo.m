function obj = defineCellInfo(obj, spikeIDs, chanIDs, probeIDs)
obj.CellInfo.NumCells = length(spikeIDs);
obj.CellInfo.CellList = 1:obj.CellInfo.NumCells;
obj.CellInfo.CellListFull = 1:obj.CellInfo.NumCells;
obj.CellInfo.NumCellsFull = length(spikeIDs);
obj.CellInfo.Shank = cell2mat(chanIDs);
obj.CellInfo.CellListString = [];
obj.CellInfo.MUAcluster = false(1,numel(spikeIDs));
obj.CellInfo.Goodcluster = false(1,numel(spikeIDs));
obj.CellInfo.Unsortedcluster = false(1,numel(spikeIDs));
obj.CellInfo.Probe = zeros(1,numel(spikeIDs));
for icell = 1:obj.CellInfo.NumCells
    str = spikeIDs{icell};
    if iscell(str)
        str = cell2mat(str);
    end
    stridx = strfind(str,'MUA');
    if ~isempty(stridx)
        obj.CellInfo.MUAcluster(icell) = true;
    end
    if isempty(stridx)
        stridx = strfind(str,'mua');
        if ~isempty(stridx)
            obj.CellInfo.Goodcluster(icell) = true;
        end
    end
    if isempty(stridx)
        stridx = strfind(str,'Good');
        if ~isempty(stridx)
            obj.CellInfo.Goodcluster(icell) = true;
        end
    end 
    if isempty(stridx)
        stridx = strfind(str,'good');
        if ~isempty(stridx)
            obj.CellInfo.Goodcluster(icell) = true;
        end
    end 
    if isempty(stridx)
        stridx = strfind(str,'Unsorted');
        if ~isempty(stridx)
            obj.CellInfo.Unsortedcluster(icell) = true;
        end
    end
    if isempty(stridx)
        stridx = strfind(str,'unsorted');
        if ~isempty(stridx)
            obj.CellInfo.Unsortedcluster(icell) = true;
        end
    end
    if strcmp(probeIDs{icell},'CA1') || strcmp(probeIDs{icell},'1')
        obj.CellInfo.Probe(icell) = 1;
    end
    if strcmp(probeIDs{icell},'V1') || strcmp(probeIDs{icell},'2')
        obj.CellInfo.Probe(icell) = 2;
    end
%     obj.CellInfo.CellListString{icell} = num2str(icell);
    if strcmp(probeIDs{icell},'CA1') || strcmp(probeIDs{icell},'1')
        obj.CellInfo.CellListString{icell} = [num2str(icell) ': tet#' num2str(obj.CellInfo.Shank(icell)) '_' str(stridx:end)] ; 
    else
        obj.CellInfo.CellListString{icell} = [num2str(icell) ': V1_' str(stridx:end)] ; 
    end
end
end