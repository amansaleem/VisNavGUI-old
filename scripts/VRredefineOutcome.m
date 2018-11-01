function es = VRredefineOutcome(es)
tolerancediff = 0.1;

es.outcome = zeros(size(es.outcome)) + 2;
roomlength = round(max(es.traj));
es.firstgoodlick = 0 * es.lick;
es.firstrewlick = 0 * es.lick;
es.goodlick = 0 * es.lick;
es.preRewlick = 0 * es.lick;
es.postRewlick = 0 * es.lick;
es.passivelick = 0 * es.lick;
es.firstbadlick = 0 * es.lick;
es.badlick = 0 * es.lick;
rewPos = min(es.rewardPos);
unwraptraj = unwrap(es.traj/roomlength*2*pi)*roomlength/(2*pi);
difflicks = -1 * ones(size(es.traj,1) + 1, 1);
difflicks(end) = inf;
difflicks((es.lick == 1)) = [inf ; abs(diff(unwraptraj((es.lick == 1))))];
for rewval = 2:-1:0
    rewidx = find(ismember(es.reward,rewval));%we take into account both passive and active reward as passive are usually in the reward location
    rewidx = [rewidx;inf(10,1)];
    for i = 1:numel(rewidx)-10
        try
        if es.reward(rewidx(i)) == 2
            firstrewlickidx = max(1,rewidx(i) - es.sampleRate) + find(es.lick(max(1,rewidx(i) - es.sampleRate):rewidx(i)), 1, 'last') - 1;
        else
            firstrewlickidx = rewidx(i) + find(es.lick(rewidx(i):min(end, rewidx(i) + 2*es.sampleRate)), 1, 'first') - 1;
        end
        if ~isempty(firstrewlickidx)
            lastlickidx = firstrewlickidx + find(difflicks(firstrewlickidx+1:min(end,rewidx(i+5))) > tolerancediff * roomlength, 1, 'first') - 1;
            if isempty(lastlickidx)
                lastlickidx = min(numel(es.lick),rewidx(i+1));
            end
            if i > 1
                firstlickidx = rewidx(i-1) + find(difflicks(rewidx(i-1)+1:firstrewlickidx) > tolerancediff * roomlength, 1, 'last');
            else
                firstlickidx = find(difflicks(1:firstrewlickidx) > tolerancediff * roomlength, 1, 'last');
            end
            if rewval < 2
                firstlickidx = firstrewlickidx;
            end
            if es.reward(rewidx(i)) == 2
                es.firstrewlick(firstrewlickidx) = 1;
                es.firstgoodlick(firstlickidx) = 1;
                es.preRewlick(firstlickidx:firstrewlickidx) = es.lick(firstlickidx:firstrewlickidx);
                es.postRewlick(firstrewlickidx + 1:lastlickidx) = es.lick(firstrewlickidx + 1:lastlickidx);
                es.goodlick(firstlickidx:lastlickidx) = es.lick(firstlickidx:lastlickidx);
            else
                es.passivelick(firstlickidx:lastlickidx) = es.lick(firstlickidx:lastlickidx);
            end        
        end
        catch
            keyboard
        end
    end
end
es.passivelick(es.goodlick > 0) = 0;
es.badlick(es.lick & ~es.goodlick & ~es.passivelick) = 1;
badlickidx = find(es.badlick == 1);
if ~isempty(badlickidx)
    diffbadlicks = [inf ; abs(diff(unwraptraj(badlickidx)))];
    es.firstbadlick(badlickidx(diffbadlicks > tolerancediff * roomlength)) = 1;
end

badlickidx = find(es.badlick == 1);
badlickidx = [badlickidx ; numel(es.badlick)];
rewlick = es.goodlick + es.passivelick;
for i = 1:numel(badlickidx) - 1
%     prevgoodlickidx = badlickidx(i) - find(rewlick(badlickidx(i):-1:max(1,badlickidx(i)-900)), 1, 'first') + 1;
%     if isempty(prevgoodlickidx)
        prevlickidx = badlickidx(i) - find(es.badlick(badlickidx(i)-1:-1:max(1,badlickidx(i)-60*es.sampleRate)) + rewlick(badlickidx(i)-1:-1:max(1,badlickidx(i)-60*es.sampleRate)), 1, 'first');        
%     end
    nextgoodlickidx = badlickidx(i) + find(rewlick(badlickidx(i)+1:badlickidx(i+1)), 1, 'first');
    nextlapidx = badlickidx(i) + find(unwraptraj(badlickidx(i)+1:badlickidx(i+1)) - unwraptraj(badlickidx(i)) > roomlength);
    if isempty(nextlapidx)
        nextlapidx = inf;
    end
%     prevlapidx = badlickidx(i) + find(unwraptraj(badlickidx(i)) - unwraptraj(badlickidx(i)-1:-1:badlickidx(i-1)) > roomlength);
%     if isempty(prevlapidx)
        prevlapidx = 1;
%     end
    if ~isempty(nextgoodlickidx)
        if nextgoodlickidx - 1 < nextlapidx
            es.outcome(max(prevlapidx,badlickidx(i)):min(nextgoodlickidx - 1,nextlapidx)) = 1;
        else
            es.outcome(max(prevlapidx,badlickidx(i)):min(nextgoodlickidx - 1,nextlapidx)) = 0;
        end
    else
        es.outcome(max(prevlapidx,badlickidx(i)):min(badlickidx(i+1),nextlapidx)) = 0;
    end
    if isempty(prevlickidx)
        trialstartidx = find(es.trialID(badlickidx(i)-1:-1:max(1,badlickidx(i)-60*es.sampleRate)) == es.trialID(badlickidx(i)), 1, 'last');
        es.outcome(badlickidx(i)-1:-1:max(1,badlickidx(i)-trialstartidx)) = 0;
    else    
        if rewlick(prevlickidx)
            if unwraptraj(badlickidx(i)) - unwraptraj(prevlickidx) < roomlength
                es.outcome(prevlickidx+1:badlickidx(i)) = 3;
            elseif unwraptraj(badlickidx(i)) - unwraptraj(prevlickidx) < 2*roomlength
                es.outcome(prevlickidx+1:badlickidx(i)) = 4;
            else
                es.outcome(prevlickidx+1:badlickidx(i)) = 5; 
            end
        elseif es.badlick(prevlickidx)
            if unwraptraj(badlickidx(i)) - unwraptraj(prevlickidx) < 2*roomlength
                es.outcome(prevlickidx+1:badlickidx(i)) = 0;
            else
                es.outcome(prevlickidx+1:badlickidx(i)) = 5; 
            end
        end
    end
end

passivelickidx = find(es.passivelick == 1);
passivelickidx = [passivelickidx ; numel(es.passivelick)];
rewlick = es.goodlick + es.badlick;
for i = 1:numel(passivelickidx) - 1
    nextgoodlickidx = passivelickidx(i) + find(rewlick(passivelickidx(i)+1:passivelickidx(i+1)), 1, 'first');
    nextlapidx = passivelickidx(i) + find(unwraptraj(passivelickidx(i)+1:passivelickidx(i+1)) - unwraptraj(passivelickidx(i)) > roomlength);
    if isempty(nextlapidx)
        nextlapidx = inf;
    end
    if ~isempty(nextgoodlickidx)
        es.outcome(passivelickidx(i):min(nextgoodlickidx - 1,nextlapidx)) = 2;%6
    else
        es.outcome(passivelickidx(i):min(passivelickidx(i+1),nextlapidx)) = 2;%6
    end
end

Alltrials = 1:max(es.trialID);
noLickstrials = Alltrials(~ismember(Alltrials, unique(es.trialID(es.lick == 1))) & ismember(Alltrials, unique(es.trialID(es.outcome == 2))));
es.outcome(ismember(es.trialID, noLickstrials)) = 5;


blanks = false(size(es.traj));
blankstart = find(diff([0;es.blanks]) > 0);
maxblank = round(max(es.blanks));
for i = 1:numel(blankstart)
    blankEnd = blankstart(i)+1 + find(es.trajspeed(blankstart(i)+1:min(end,blankstart(i)+1 + 50*maxblank)) > 1, 1, 'first');
    if isempty(blankEnd)
        blankEnd = numel(es.traj);
    end
    blankStart = max(1,blankstart(i) - floor(es.sampleRate/12) + 1 - find(abs(es.traj(max(1,blankstart(i) - floor(es.sampleRate/12):-1:blankstart(i)-es.sampleRate)) - es.traj(blankEnd)) < floor(es.sampleRate/12), 1, 'first'));
    blanks(blankStart:blankEnd) = true;
end
% blanks(es.traj == 0 & es.blanks > 1) = true;
es.blanks = blanks;

es.afterblanks = false(size(es.traj));
trialID = unique(es.trialID);
for trial = 1:numel(trialID)
    trialidx = find(es.trialID == trialID(trial));
    if sum(es.blanks(trialidx))>0
        lastblank = find(es.blanks(trialidx)>0, 1, 'last');
        es.afterblanks(trialidx(lastblank:end)) = true;
    end
end

es.completetrial = false(size(es.traj));
mintraj = min(es.traj);
maxtraj = max(es.traj);
trialID = unique(es.trialID);
for trial = 1:numel(trialID)
    minxtrial = min(es.traj(es.trialID == trialID(trial) & ~es.blanks));
    maxxtrial = max(es.traj(es.trialID == trialID(trial) & ~es.blanks));
    if ~isempty(minxtrial) && ~isempty(maxxtrial)
        if minxtrial<mintraj + 0.1*(maxtraj - mintraj) && maxxtrial>maxtraj - 0.1*(maxtraj - mintraj)
            es.completetrial(es.trialID == trialID(trial)) = true;
        end
    end
end

trialnum = unique(es.trialID);
difftraj = [es.traj(1);diff(unwrap(es.traj/round(max(es.traj))*2*pi)*round(max(es.traj))/(2*pi))];
difftraj(abs(difftraj) > max(es.traj)/2) = 0;
difftraj(1) = es.traj(1);
Samples = 1/es.sampleRate*ones(size(es.traj));
ballgain = 1;%1.1;

difftraj(es.blanks) = 0;
Samples(es.blanks) = 0;

es.currTime = zeros(size(es.traj));
es.runTime = zeros(size(es.traj));
for tt = 1:numel(trialnum)
    es.distTrav(es.trialID == trialnum(tt)) = es.traj(es.trialID == trialnum(tt))./es.gain(es.trialID == trialnum(tt))/ballgain;%cumsum(difftraj(es.trialID == trialnum(tt)).*(1./es.gain(es.trialID == trialnum(tt)))/ballgain);
    es.currTime(es.trialID == trialnum(tt)) = cumsum(Samples(es.trialID == trialnum(tt)));
    es.runTime(es.trialID == trialnum(tt)) = cumsum(Samples(es.trialID == trialnum(tt)).*(difftraj(es.trialID == trialnum(tt)) > 0));
end
es.distTrav(es.distTrav < 0) = 0;

% for tt = 1:numel(trialnum)
%     es.distfromStart(es.trialID == trialnum(tt)) = es.traj
% end

end