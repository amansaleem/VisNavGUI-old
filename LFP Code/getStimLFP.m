function [out] = getStimLFP(es, p, stimTimes, dT)

if nargin<4
    dT = 1; % 1 sec before and after the stimulus
end
sampleRate = 30000;
% if nargin<4
%     var = 6; % 6 for contrast in vs stimuli
% end
numStim = length(p.pfiledurs);
numReps = p.nrepeats;
for iStim = 1:numStim
    for iRep = 1:numReps
        startT = stimTimes.on(iStim, iRep)./sampleRate - dT;
        endT   = stimTimes.off(iStim, iRep)./sampleRate + dT;
        
        startBin = max(find(es.t<startT));
        endBin   = min(find(es.t>endT));
        
        try
            powA(iStim, iRep,:,:) = es.powA(startBin:endBin,:);
        catch
            powA(iStim, iRep,:,:) = nan;
        end
        %             out(iStim, iRep).powB = es.powB(startBin:endBin,:);
        
    end
end
out = squeeze(nanmean(powA,2));
