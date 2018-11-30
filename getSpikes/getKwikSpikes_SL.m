function chans = getKwikSpikes_SL(animal, iseries, iexp, igroup, addInfo)

SetDefaultDirs2018

DIRname  = [DIRS.multichanspikes filesep num2str(iseries) filesep];

% check if the chans data structure for each recording and every shank
% (igroup) have been saved already
if  exist([DIRname 'SpikeTimes_Matlab' filesep animal '_' num2str(iseries) '_' num2str(iexp) '_' num2str(igroup) '.mat'],'file')
    load([DIRname 'SpikeTimes_Matlab' filesep animal '_' num2str(iseries) '_' num2str(iexp) '_' num2str(igroup) '.mat']);
else 
    display('Spiketimes have not been converted yet or some error reading the saved file. Loading them now...');
    
    if nargin>4 & ~isempty(addInfo)
        basename = [animal '_' num2str(iseries) '_' addInfo];
    else
        basename = [animal '_' num2str(iseries)];
    end
    
%     load([DIRname basename '.dat_meta.mat']);
%     kwikFile = [DIRname basename '.kwik'];
    
    FileNameMetaDat = dir([DIRname '*.dat_meta.mat']);
    load([DIRname FileNameMetaDat.name]);
    FileNameKwik = dir([DIRname '*.kwik']);
    kwikFile = [DIRname FileNameKwik.name];
    
    expIdx = [];
%     while isempty(expIdx)
%         userinput = input('Please insert number of i-th number of ePhys recording coincidental with visual stimulus of interest: \n','s');
%         expIdx = str2num(userinput);
%     end
    %expIdx = find(SELECTED_EXPERIMENTS==iexp);
%     expEnds = cumsum(lims);
%     if expIdx>1
%         startTime = expEnds(expIdx-1);
%     else
%         startTime = 1;
%     end
%     endTime = expEnds(expIdx);
    
    spkTimes = hdf5read(kwikFile, ['/channel_groups/' num2str(igroup) '/spikes/time_samples']);
    spkClus = hdf5read(kwikFile, ['/channel_groups/' num2str(igroup) '/spikes/clusters/main']);
    cellIDs = unique(spkClus);
    ncells  = length(cellIDs);
    
    
%     spkClus(spkTimes>endTime) = [];
%     spkTimes(spkTimes>endTime) = [];
%     
%     spkTimes = double(spkTimes);
%     spkTimes = spkTimes - startTime;
%     spkClus(spkTimes<0) = [];
%     spkTimes(spkTimes<0) = [];
    
    temp =  h5info(kwikFile,  '/recordings/0/');
    sampleRate = double((temp.Attributes(3).Value));
    spkTimes = double(spkTimes)./sampleRate;
    
    noise_list = [];
    
    for icell = 1:ncells
        chans(icell).spiketimes = spkTimes(spkClus==cellIDs(icell));
        chans(icell).ichan = igroup;
        chans(icell).iexp = iexp;
        chans(icell).icell = cellIDs(icell);
        chans(icell).sampleRate = sampleRate;
        
        temp = h5info(kwikFile, ['/channel_groups/' num2str(igroup) '/clusters/main/' num2str(cellIDs(icell)) '/']);
        
        for idx = 1:length(temp.Attributes)
            if strcmp(temp.Attributes(idx).Name, 'cluster_group')
                ilabel = temp.Attributes(idx).Value;
                break
            end
            idx = idx + 1;
        end
        if idx <= length(temp.Attributes)
            ilabel = temp.Attributes(idx).Value;
            temp = h5info(kwikFile, ['/channel_groups/' num2str(igroup) '/cluster_groups/main']);
            labelType = eval(['temp.Groups(' num2str(ilabel+1) ').Attributes(end).Value']);
            if strcmp(labelType,'Noise')
                noise_list = [noise_list icell];
            end
            chans(icell).id = [basename '_c' num2str(igroup) '_' labelType '_cluster' num2str(cellIDs(icell))];
        else
            chans(icell).id = [basename '_c' num2str(igroup) '_none_cluster' num2str(cellIDs(icell))];
        end
    end
    chans(noise_list) = [];
    if ~exist([DIRname 'SpikeTimes_Matlab'],'dir')
        mkdir([DIRname 'SpikeTimes_Matlab']);
    end
    save([DIRname 'SpikeTimes_Matlab' filesep animal '_' num2str(iseries) '_' num2str(iexp) '_' num2str(igroup) '.mat'],'chans');
end

