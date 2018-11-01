function [trig_pow, all_trig_pow] = calculateTrigPower(input, trig, pre_bins, post_bins)
% gray_starts = find(k<-50 & ~isnan(es.traj));

if nargin<4
    pre_bins = 10;
    post_bins = 15;
end

numBins = (pre_bins+post_bins);
start_idx = max([1 max(find(trig<pre_bins))+1]);
end_idx   = max(find((size(input,1) - trig)>post_bins));

% subtracting the mean power
% input = input - repmat(nanmean(input,1),size(input,1),1);

for ibin = 1:numBins
    trig_pow(ibin,:) = nanmean(input(trig(start_idx:end_idx)+ibin-pre_bins,:));
end

if nargout>1
    newTrig = trig(start_idx:end_idx);
    for ibin = 1:numBins
        for itrig = 1:length(newTrig)
            all_trig_pow(itrig, ibin, :) = input(newTrig(itrig)+ibin-pre_bins,:);
        end
    end
end