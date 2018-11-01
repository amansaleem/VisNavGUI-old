function chans = getKiloSortSpikes(animal, iseries, iexp, igroup, addInfo)
% [filenum;spkevtTime;chind;elind;crossmax;unitID;unitID;chi2;sqrsum;zeros(1,numel(spkevtTime));zeros(1,numel(spkevtTime));zeros(1,numel(spkevtTime))]
SetDefaultDirs
% DIRS.spikes = '\\zserver\Data\Spikes';
% DIRS.multichanspikes = '\\zserver\Data\multichanspikes';

basename = [animal '_s' num2str(iseries) '_1'];
basenameKiloSort = [animal '_s' num2str(iseries) '_1'];
if nargin>4
    if ~isempty(addInfo)
        basename = [animal '_s' num2str(iseries) '_' addInfo];
        basenameKiloSort = [animal '_s' num2str(iseries) '_' addInfo];
    end
end

DIRname  = [DIRS.multichanspikes filesep animal filesep num2str(iseries) filesep];
if isdir([DIRname addInfo filesep 'tet' num2str(igroup)])
    load([DIRname addInfo filesep basename]);
    KiloSortFile_times = [DIRname addInfo filesep 'Kilosort' filesep 'tet' num2str(igroup) filesep 'spike_times.npy'];
    KiloSortFile_templates = [DIRname addInfo filesep 'Kilosort' filesep 'tet' num2str(igroup) filesep 'spike_templates.npy'];
    KiloSortFile_clusters = [DIRname addInfo filesep 'Kilosort' filesep 'tet' num2str(igroup) filesep 'spike_clusters.npy'];
else
    load([DIRname addInfo filesep basename]);
    KiloSortFile_times = [DIRname addInfo filesep 'Kilosort' filesep 'spike_times.npy'];
    KiloSortFile_templates = [DIRname addInfo filesep 'Kilosort' filesep 'spike_templates.npy'];
    KiloSortFile_clusters = [DIRname addInfo filesep 'Kilosort' filesep 'spike_clusters.npy'];
end

expIdx = find(SELECTED_EXPERIMENTS==iexp);
expEnds = cumsum(lims);
if expIdx>1
    startTime = expEnds(expIdx-1);
else
    startTime = 1;
end
endTime = expEnds(expIdx);

times = readNPY(KiloSortFile_times);
templates = readNPY(KiloSortFile_templates);
clusters = readNPY(KiloSortFile_clusters);
try
T = readtable([DIRname addInfo filesep 'Kilosort' filesep 'tet' num2str(igroup) filesep 'cluster_groups.csv']);
Cluster_labels = table2cell(T);
Cluster_labels = Cluster_labels(:,2);
catch
    Cluster_labels = {'unsorted'};
end

cellIDs = sort(unique(clusters),'ascend');
ncells  = numel(cellIDs);
sampleRate = 30000;
times = double(times)./sampleRate;

for icell = 1:ncells
    chans(icell).spiketimes = times(clusters == cellIDs(icell));
    chans(icell).ichan = igroup;
    chans(icell).iexp = iexp;
    chans(icell).icell = cellIDs(icell);
    chans(icell).sampleRate = sampleRate;
    chans(icell).id = [basename '_c' num2str(igroup) '_' Cluster_labels{min(icell,size(Cluster_labels,1))} '_cluster' num2str(cellIDs(icell))];
end

end