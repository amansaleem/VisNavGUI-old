function chans = getConsensusSpikes(animal, iseries, iexp, igroup, addInfo)
% [filenum;spkevtTime;chind;elind;crossmax;unitID;unitID;chi2;sqrsum;zeros(1,numel(spkevtTime));zeros(1,numel(spkevtTime));zeros(1,numel(spkevtTime))]
SetDefaultDirs
% DIRS.spikes = '\\zserver\Data\Spikes';
% DIRS.multichanspikes = '\\zserver\Data\multichanspikes';

basename = [animal '_s' num2str(iseries) '_1'];
basenameConsensus = [animal '_s' num2str(iseries)];
if nargin>4
    if ~isempty(addInfo)
        basename = [animal '_s' num2str(iseries) '_' addInfo];
        basenameConsensus = [animal '_s' num2str(iseries) '_' addInfo];
    end
end

DIRname  = [DIRS.multichanspikes filesep animal filesep num2str(iseries) filesep];
if isdir([DIRname addInfo])    
    load([DIRname addInfo filesep basename]);
    ConsensusFile = [DIRname addInfo filesep 'tet' num2str(igroup) filesep basenameConsensus '_tet' num2str(igroup) '_Spkevent.mat'];
    igroupfile = 0;
else
    load([DIRname basename]);
    ConsensusFile = [DIRname basenameConsensus '.mat'];
    igroupfile = igroup;
end

expIdx = find(SELECTED_EXPERIMENTS==iexp);
expEnds = cumsum(lims);
if expIdx>1
    startTime = expEnds(expIdx-1);
else
    startTime = 1;
end
endTime = expEnds(expIdx);

S = load(ConsensusFile);
fileidx = find(diff(S.spkevent(2,:)) < 0);
fileidx = [fileidx size(S.spkevent,2)];
Spkevent = S.spkevent;
clusterID = S.spkclustIDmulti;
clusterID(sum(S.spkclustIDmulti,2) == 0,:) = [];
if numel(fileidx) > 1
    for f = 1:numel(fileidx)-1
        Spkevent(2,fileidx(f)+1:fileidx(f+1)) = S.spkevent(2,fileidx(f)+1:fileidx(f+1)) + sum(S.spkevent(2,fileidx(1:f)));
    end
end
Spkevent(2,:) = Spkevent(2,:)*10^-6*30000;
cellIDs = 1:S.nbunitmultich;
ncells  = S.nbunitmultich;


Spkevent(:,Spkevent(2,:) > endTime) = [];
Spkevent(2,:) = Spkevent(2,:) - startTime;
Spkevent(:,Spkevent(2,:) < 0) = [];

sampleRate = 30000;
Spkevent(2,:) = double(Spkevent(2,:))./sampleRate;

for icell = 1:ncells
    cellID = clusterID(icell,clusterID(icell,:) > 0);
    cellidx = ismember(Spkevent(6,:),cellID);
    chans(icell).spiketimes = Spkevent(2, cellidx);
    chans(icell).ichan = Spkevent(4, find(cellidx,1, 'first'));
    chans(icell).iexp = Spkevent(1, find(cellidx,1, 'first'));
    chans(icell).icell = icell;
    chans(icell).sampleRate = sampleRate;
    chans(icell).id = num2str(icell);
end

end