function [pow_es] = getPowSpd(es, t, spd_bins, variable, numIter)

if nargin>3
    spd = variable;
else
    spd = es.smthBallSpd;
end

if nargin<5
    numIter = 100;
end
if nargin<3 | isempty(spd_bins)
    spd_bins = [0 1 3 8 14 22 32];
end
if nargin<2
    spd_bins = [0 1 3 8 14 22 32];
    t = true(size(es.ballspeed));
end

while max(spd)<spd_bins(end)
    spd_bins(end) = [];
end
numSpdBins = length(spd_bins);
skip_bins = [];
for iSpdBin = 1:numSpdBins
    if iSpdBin==numSpdBins
        if sum(t & (spd>spd_bins(iSpdBin)))==0
            skip_bins = [skip_bins iSpdBin];
        end
    else
        if sum(t & (spd>spd_bins(iSpdBin)) & (spd<=spd_bins(iSpdBin+1)))==0
            skip_bins = [skip_bins iSpdBin];
        end
    end 
end
spd_bins(skip_bins) = [];

numSpdBins = length(spd_bins);
pow_es.freq = es.freq;

pow_es.booterror.coherence  = zeros(numIter, numSpdBins, length(es.freq));
pow_es.booterror.powA       = zeros(numIter, numSpdBins, length(es.freq));
pow_es.booterror.powB       = zeros(numIter, numSpdBins, length(es.freq));
pow_es.booterror.powAB      = zeros(numIter, numSpdBins, length(es.freq));

if numSpdBins>1
    if isfield(es,'pairs')
       for ibin = 1:numSpdBins
        if ibin==numSpdBins
            t_bin = t & (spd>spd_bins(ibin));
        else
            t_bin = t & (spd>spd_bins(ibin) & spd<=spd_bins(ibin+1));
        end
        pow_es.coherence(ibin,:) = nanmean(es.pairs(1).C(t_bin,:),1);
        pow_es.powA(ibin,:)      = nanmean(es.pairs(1).SA(t_bin,:),1);
        pow_es.powB(ibin,:)      = nanmean(es.pairs(1).SB(t_bin,:),1);
        pow_es.powAB(ibin,:)     = nanmean(es.pairs(1).SAB(t_bin,:),1);
        pow_es.spdBins(ibin)     = nanmean(spd(t_bin));
        pow_es.timeInBins(ibin)  = nansum(t_bin);
        
        numBins = nansum(t_bin);
        binPos = find(t_bin);
        
        resample_idx = ceil(rand(numIter,numBins)*numBins);
        resample_distr = binPos(resample_idx); 
        display(['SpdBin: ' num2str(ibin)]);
        for biter = 1:numIter
%             display(['Iter: ' num2str(biter)]);
            sample = resample_distr(biter,:);
            pow_es.booterror.coherence(biter,ibin,:) = nanmean(es.pairs(1).C(sample,:),1);
            pow_es.booterror.powA(biter,ibin,:)      = nanmean(es.pairs(1).SA(sample,:),1);
            pow_es.booterror.powB(biter,ibin,:)      = nanmean(es.pairs(1).SB(sample,:),1);
            pow_es.booterror.powAB(biter,ibin,:)     = nanmean(es.pairs(1).SAB(sample,:),1);
        end
    end 
    else
        for ibin = 1:numSpdBins
        if ibin==numSpdBins
            t_bin = t & (spd>spd_bins(ibin));
        else
            t_bin = t & (spd>spd_bins(ibin) & spd<=spd_bins(ibin+1));
        end
        pow_es.coherence(ibin,:) = nanmean(es.coherence(t_bin,:),1);
        pow_es.powA(ibin,:)      = nanmean(es.powA(t_bin,:),1);
        pow_es.powB(ibin,:)      = nanmean(es.powB(t_bin,:),1);
        pow_es.powAB(ibin,:)     = nanmean(es.powAB(t_bin,:),1);
        pow_es.spdBins(ibin)     = nanmean(spd(t_bin));
        pow_es.timeInBins(ibin)  = nansum(t_bin);
        
        pow_es.corr(ibin).AA(:,:) = corr(es.powA(t_bin,:),es.powA(t_bin,:));
        pow_es.corr(ibin).BB(:,:) = corr(es.powB(t_bin,:),es.powB(t_bin,:));
        pow_es.corr(ibin).AB(:,:) = corr(es.powA(t_bin,:),es.powB(t_bin,:));
        
        numBins = nansum(t_bin);
        binPos = find(t_bin);
        
        resample_idx = ceil(rand(numIter,numBins)*numBins);
        resample_distr = binPos(resample_idx); 
        display(['SpdBin: ' num2str(ibin)]);
        for biter = 1:numIter
%             display(['Iter: ' num2str(biter)]);
            sample = resample_distr(biter,:);
            pow_es.booterror.coherence(biter,ibin,:) = nanmean(es.coherence(sample,:),1);
            pow_es.booterror.powA(biter,ibin,:)      = nanmean(es.powA(sample,:),1);
            pow_es.booterror.powB(biter,ibin,:)      = nanmean(es.powB(sample,:),1);
            pow_es.booterror.powAB(biter,ibin,:)     = nanmean(es.powAB(sample,:),1);
        end
        end
    end
end