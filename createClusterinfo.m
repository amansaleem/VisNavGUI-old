function createClusterinfo(animal, iseries, Probe, NCh, group)
% get the average waveform of each cluster and compute some measure on
% clusters before saving it to a mat file
SetDefaultDirs;
global DIRS


if nargin < 3
    Probe = 'CA1';
    group = 'tet0';
    NCh = 4;
end
if nargin < 5
    group = [];
end
for s = 1:numel(iseries)
    exptInfo.animal  = animal;
    exptInfo.iseries = iseries(s);
    
    %% Load the data file
    if ~isempty(group)
        expName = [exptInfo.animal '_s' num2str(exptInfo.iseries) '_' Probe '_' group]
        datfilePath = [DIRS.multichanspikes filesep exptInfo.animal filesep num2str(exptInfo.iseries) filesep Probe filesep group filesep expName '.dat'];
        kwikfilePath = [DIRS.multichanspikes filesep exptInfo.animal filesep num2str(exptInfo.iseries) filesep Probe filesep group filesep expName '.kwik'];
        savedfilePath = [DIRS.multichanspikes filesep exptInfo.animal filesep num2str(exptInfo.iseries) filesep Probe filesep group filesep expName '_clusterinfo.mat'];
    else
        expName = [exptInfo.animal '_s' num2str(exptInfo.iseries) '_' Probe]
        datfilePath = [DIRS.multichanspikes filesep exptInfo.animal filesep num2str(exptInfo.iseries) filesep Probe filesep expName '.dat'];
        kwikfilePath = [DIRS.multichanspikes filesep exptInfo.animal filesep num2str(exptInfo.iseries) filesep Probe filesep expName '.kwik'];
        savedfilePath = [DIRS.multichanspikes filesep exptInfo.animal filesep num2str(exptInfo.iseries) filesep Probe filesep expName '_clusterinfo.mat'];
    end
    fid = fopen(datfilePath, 'r');
    fseek(fid, 0, 'bof');
    Vdata = fread(fid, '*int16');
    Vdata = reshape(Vdata, NCh, numel(Vdata)/NCh);
    fclose(fid);
    
    spkTimes = hdf5read(kwikfilePath, ['/channel_groups/0/spikes/time_samples']);
    spkClus = hdf5read(kwikfilePath, ['/channel_groups/0/spikes/clusters/main']);
    cellIDs = unique(spkClus);
    ncells  = length(cellIDs);
    

    for icell = 1:ncells
        ispktimes = spkTimes(spkClus == cellIDs(icell));
        clus(icell).waveform = zeros(NCh,61);
        for dt = -30:30
            clus(icell).waveform(:,dt+31) = mean(Vdata(:,max(1,min(ispktimes + dt,end))),2);
        end
        
        [~, clus(icell).bestchan] = max(max(abs(clus(icell).waveform), [], 2));
        clus(icell).igroup = group;
        clus(icell).icell = cellIDs(icell);
        
        temp = h5info(kwikfilePath, ['/channel_groups/0/clusters/main/' num2str(cellIDs(icell)) '/']);

        for idx = 1:length(temp.Attributes)
            if strcmp(temp.Attributes(idx).Name, 'cluster_group')
                ilabel = temp.Attributes(idx).Value;
                break
            end
            idx = idx + 1;
        end
        if idx <= length(temp.Attributes)
            ilabel = temp.Attributes(idx).Value;
            temp = h5info(kwikfilePath, ['/channel_groups/0/cluster_groups/main']);
            clus(icell).labelType = eval(['temp.Groups(' num2str(ilabel+1) ').Attributes(end).Value']);
        else
            clus(icell).labelType = {'undefined'};
        end
    end
    save(savedfilePath,'clus');
end


% expt = getExperimentList;
% for ianimal = 1:numel(expt)
%     for iseries = 1:numel(expt(ianimal).series)
%         if expt(ianimal).CA1{iseries}
%             for itet = 0:7
%                 createClusterinfo(expt(ianimal).animal, expt(ianimal).series{iseries}, 'CA1', 4, ['tet' num2str(itet)]);
%             end
%         end
%         if expt(ianimal).V1{iseries}
%             createClusterinfo(expt(ianimal).animal, expt(ianimal).series{iseries}, 'V1', 32);
%         end
%     end
% end
