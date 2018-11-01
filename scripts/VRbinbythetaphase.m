function es_out = VRbinbythetaphase(es_in,nthetaphsbins,chn)
sampleFrequency = 1/es_in.sampleRate;
if isfield(es_in, 'LFPphase')
    if isempty(nthetaphsbins) || nthetaphsbins == 0
        es_out = es_in;
        es_out.sampleSize = sampleFrequency * ones(size(es_out.sampleTimes));
        phsval = mod(es_in.LFPphase(:,min(end,chn))-180,360);
        es_out.LFPbinphase = phsval;
    else        
        LFPbinphase = zeros(size(es_in.LFPphase(:,min(end,chn))));
        if nthetaphsbins > 1
            phsval = mod(es_in.LFPphase(:,min(end,chn))-180,360);
            for phs = 1:nthetaphsbins
                phs0 = (phs-1)*360/nthetaphsbins;        
                phsidx = ismember(phsval,mod(phs0:(phs0 + 360/nthetaphsbins),360));
                LFPbinphase(phsidx) = phs0;
            end
            idxphschange = [1 ; find(diff(LFPbinphase)~=0)+1 ; numel(LFPbinphase)+1];
        else
            phsval = mod(es_in.LFPphase(:,min(end,chn))-180,360);
            idxphschange = [1 ; find(diff(phsval)< -180)+1 ; numel(LFPbinphase)+1];
            LFPbinphase(:) = 0;
        end
        
%         maxtraj = round(max(es_in.traj));
%         bintraj = zeros(size(es_in.traj));
%         if nthetaphsbins > 1
%             phsval = round(es_in.traj);
%             for x = 1:nthetaphsbins
%                 phs0 = (x-1)*maxtraj/nthetaphsbins;        
%                 phsidx = ismember(phsval,mod(phs0:(phs0 + maxtraj/nthetaphsbins),maxtraj));
%                 bintraj(phsidx) = phs0;
%             end
%             idxphschange = [1 ; find(diff(bintraj)~=0)+1 ; numel(bintraj)+1];
%             LFPbinphase = mod(es_in.LFPphase(:,chn)-180,360);
%         else
%             phsval = mod(es_in.LFPphase(:,chn)-180,360);
%             idxphschange = [1 ; find(diff(phsval)< -180)+1 ; numel(LFPbinphase)+1];
%             LFPbinphase(:) = 0;
%         end


    %     phsval = mod(es_in.LFPphase(:,chn)-180,360);
    %     LFPbinphase = phsval;
    %     idxphschange = [(1:2:numel(LFPbinphase))' ; numel(LFPbinphase)+1];

        nbins = numel(idxphschange)-1;
        maxtraj = round(max(es_in.traj));
        maxtrajPercent = round(max(es_in.trajPercent));                 
         
        es_out.traj = zeros(nbins,1);
        es_out.sampleTimes = zeros(nbins,1);
        es_out.lick = zeros(nbins,1);
        es_out.reward = zeros(nbins,1);
        es_out.trialID = zeros(nbins,1);
        es_out.traj = zeros(nbins,1);
        es_out.trajspeed = zeros(nbins,1);
        es_out.ballspeed = zeros(nbins,1);
        es_out.currTime = zeros(nbins,1);
        es_out.distTrav = zeros(nbins,1);
        es_out.totDistTrav = zeros(nbins,1);
        es_out.trajPercent = zeros(nbins,1);
        es_out.smthBallSpd = zeros(nbins,1);
        es_out.smthBallAcc = zeros(nbins,1);
        es_out.smthTrajSpd = zeros(nbins,1);
        es_out.contrast = zeros(nbins,1);
        es_out.start = zeros(nbins,1);
        es_out.gain = zeros(nbins,1);
        es_out.blanks = zeros(nbins,1);
        es_out.afterblanks = zeros(nbins,1);
        es_out.active = zeros(nbins,1);
        es_out.rewardPos = zeros(nbins,1);
        es_out.outcome = zeros(nbins,1);
        es_out.roomLength = zeros(nbins,1);
        es_out.firstgoodlick = zeros(nbins,1);
        es_out.firstrewlick = zeros(nbins,1);
        es_out.goodlick = zeros(nbins,1);
        es_out.preRewlick = zeros(nbins,1);
        es_out.postRewlick = zeros(nbins,1);
        es_out.passivelick = zeros(nbins,1);
        es_out.firstbadlick = zeros(nbins,1);
        es_out.badlick = zeros(nbins,1);
        es_out.runTime = zeros(nbins,1);
        es_out.iexp = zeros(nbins,1);
        es_out.series = zeros(nbins,1);
        es_out.LFPphase = zeros(nbins,size(es_in.LFPphase,2));
        es_out.LFPpower = zeros(nbins,size(es_in.LFPpower,2));
        es_out.LFPtheta = zeros(nbins,size(es_in.LFPtheta,2));
        es_out.spikeTrain = zeros(nbins,size(es_in.spikeTrain,2));
        es_out.mua = zeros(nbins,1);

        es_out.LFPbinphase = zeros(nbins,1);

        es_out.sampleSize = zeros(nbins,1);
        
        es_out.trialgainchange = zeros(nbins,1);
        es_out.eyeXpos = zeros(nbins,1);
        es_out.eyeYpos = zeros(nbins,1);
        es_out.pupilSize = zeros(nbins,1);
        es_out.halfID = zeros(nbins,1);

        for i = 1:numel(idxphschange)-1
            idxmerge = idxphschange(i):idxphschange(i+1)-1;

            es_out.traj(i) = mod(nanmean(unwrap(es_in.traj(idxmerge)/maxtraj*2*pi)*maxtraj/(2*pi)),maxtraj);
            es_out.trajPercent(i)  = mod(nanmean(unwrap(es_in.trajPercent(idxmerge)/maxtrajPercent*2*pi)*maxtrajPercent/(2*pi)),maxtrajPercent);

            es_out.sampleSize(i) = numel(idxmerge) * sampleFrequency;
            es_out.sampleTimes(i) = nanmean(es_in.sampleTimes(idxmerge));
            es_out.lick(i) = max(es_in.lick(idxmerge));
            es_out.reward(i) = max(es_in.reward(idxmerge));
            es_out.trialID(i) = es_in.trialID(idxmerge(1));
            es_out.trajspeed(i) = nanmean(es_in.trajspeed(idxmerge));
            es_out.ballspeed(i) = nanmean(es_in.ballspeed(idxmerge));
            es_out.currTime(i) = nanmean(es_in.currTime(idxmerge));
            es_out.distTrav(i) = nanmean(es_in.distTrav(idxmerge));
            es_out.totDistTrav(i) = nanmean(es_in.totDistTrav(idxmerge));
            es_out.smthBallSpd(i) = nanmean(es_in.smthBallSpd(idxmerge));
            es_out.smthBallAcc(i) = nanmean(es_in.smthBallAcc(idxmerge));
            es_out.smthTrajSpd(i) = nanmean(es_in.smthTrajSpd(idxmerge));
            es_out.eyeXpos = nanmean(es_in.eyeXpos(idxmerge));
            es_out.eyeYpos = nanmean(es_in.eyeYpos(idxmerge));
            es_out.pupilSize(i) = nanmean(es_in.pupilSize(idxmerge));

            es_out.contrast(i) = es_in.contrast(idxmerge(1));
            es_out.start(i) = es_in.start(idxmerge(1));
            es_out.gain(i) = es_in.gain(idxmerge(1));
            es_out.blanks(i) = es_in.blanks(idxmerge(1));
            es_out.afterblanks(i) = es_in.afterblanks(idxmerge(1));
            es_out.active(i) = es_in.active(idxmerge(1));
            es_out.rewardPos(i) = es_in.rewardPos(idxmerge(1));
            es_out.outcome(i) = max(es_in.outcome(idxmerge));
            es_out.roomLength(i) = es_in.roomLength(idxmerge(1));
            es_out.trialgainchange(i) = es_in.trialgainchange(idxmerge(1));

            es_out.firstgoodlick(i) = max(es_in.firstgoodlick(idxmerge));
            es_out.firstrewlick(i) = max(es_in.firstrewlick(idxmerge));
            es_out.goodlick(i) = max(es_in.goodlick(idxmerge));
            es_out.preRewlick(i) = max(es_in.preRewlick(idxmerge));
            es_out.postRewlick(i) = max(es_in.postRewlick(idxmerge));
            es_out.passivelick(i) = max(es_in.passivelick(idxmerge));
            es_out.firstbadlick(i) = max(es_in.firstbadlick(idxmerge));
            es_out.badlick(i) = max(es_in.badlick(idxmerge));

            es_out.runTime(i) = nanmean(es_in.runTime(idxmerge));
            es_out.iexp(i) = es_in.iexp(idxmerge(1));
            es_out.series(i) = es_in.series(idxmerge(1));

            es_out.LFPphase(i,:) = mod(nanmean(unwrap(es_in.LFPphase(idxmerge,:)/360*2*pi)*360/(2*pi),1),360);
            es_out.LFPpower(i,:) = nanmean(es_in.LFPpower(idxmerge,:),1);
            es_out.LFPtheta(i,:) = nanmean(es_in.LFPtheta(idxmerge,:),1);
            es_out.spikeTrain(i,:) = nansum(es_in.spikeTrain(idxmerge,:),1);
            es_out.mua(i) = nansum(es_in.mua(idxmerge));

            es_out.LFPbinphase(i) = nanmean(LFPbinphase(idxmerge)); 
            es_out.halfID = round(nanmean(LFPbinphase(idxmerge)));
        end

        es_out.isolDist = es_in.isolDist;
        es_out.spikeIDs = es_in.spikeIDs;
        es_out.chanIDs = es_in.chanIDs;
        es_out.ProbeIDs = es_in.ProbeIDs;
        es_out.CircularMaze = es_in.CircularMaze;
        es_out.sampleRate = es_in.sampleRate;
    end
else
    es_out = es_in;
    es_out.sampleSize = sampleFrequency * ones(size(es_out.sampleTimes));
    es_out.LFPbinphase = zeros(size(es_out.sampleSize ));
end
end