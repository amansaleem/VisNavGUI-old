function [VRdata, VRdata_o, es] = VRWheelLoad(animal, iseries, iexp, SmthTimeWindow)

% [VRdata, VRdata_o] = VRBallLoad_v2(animal, iseries, iexp)
% To load the ball data for the specified session (VRdata) and
% if a replay, of the replayed session
% 04/11: AS created
% 12/12: AS converted to read a wheel instead of the ball and process that
% data
% 03/13: AS added the even sampling option within.
% 09/15: JUL (change relative to the max number of trials)

global DIRS
global serverName server2Name server3Name
if isempty(DIRS); SetDefaultDirs; end

VRDIRname = [DIRS.ball filesep animal filesep num2str(iseries)];
dataDIRname = [DIRS.multichanspikes filesep animal filesep num2str(iseries)];
VRfname = [animal '_' num2str(iseries) '_session_' num2str(iexp) '_trial001'];

eyeDIRname = [DIRS.EyeCamera filesep animal filesep num2str(iseries) filesep num2str(iexp)];
eyetrackfname = [num2str(iseries) '_' num2str(iexp) '_' animal '_eye_proc.mat'];
eyelogfname = [num2str(iseries) '_' num2str(iexp) '_' animal '_eye.mat'];

if exist([VRDIRname filesep VRfname '.mat'],'file')
    VRdata = load([VRDIRname filesep VRfname]);


    if nargin < 4
        SmthTimeWindow = 150;
    end

    %% To truncate the trials to only to max trial lengths
    goodTimes = ~(VRdata.TRIAL.time==0); % non-zero times
    goodTimes(:,1) = true;
    [~,colGT,~] = find(goodTimes); % find columns with non-zero times
    maxTrialLength = max(colGT);

    % JUL - 15.09.2015: changed the way the total number of trials is defined
    % by looking for the last trials that is longer than half a second instead
    % of :
    % totTrials = min([min(find(sum(double(goodTimes'))<=60)) size(VRdata.TRIAL.traj,1)]); % the first trial that is less than one second
    totTrials = min([min(find(sum(double(goodTimes'))>30,1,'last')) size(VRdata.TRIAL.traj,1)]); % the last trial that is longer than half second
    if isfield(VRdata.TRIAL,'trialContr')
        totTrials = min([totTrials sum(~isnan(VRdata.TRIAL.trialContr))]);
    end
    if isempty(totTrials)
        totTrials = size(goodTimes,1);
    end
    goodTimes = goodTimes(1:totTrials,1:maxTrialLength);

    % Truncating all relevant variables
    VRdata.TRIAL.posdata  = VRdata.TRIAL.posdata(1:totTrials,1:maxTrialLength,:);
    VRdata.TRIAL.traj     = VRdata.TRIAL.traj(1:totTrials,1:maxTrialLength);
    VRdata.TRIAL.time     = VRdata.TRIAL.time(1:totTrials,1:maxTrialLength);
    VRdata.TRIAL.balldata = VRdata.TRIAL.balldata(1:totTrials,1:maxTrialLength,:);
    VRdata.TRIAL.trialIdx = VRdata.TRIAL.trialIdx(1:totTrials,1:maxTrialLength);

    VRdata.TRIAL.eyeXpos = zeros(size(VRdata.TRIAL.traj,1),size(VRdata.TRIAL.traj,2));
    VRdata.TRIAL.eyeYpos = zeros(size(VRdata.TRIAL.traj,1),size(VRdata.TRIAL.traj,2));
    VRdata.TRIAL.pupilSize = zeros(size(VRdata.TRIAL.traj,1),size(VRdata.TRIAL.traj,2));
    VRdata.TRIAL.eyetime = zeros(size(VRdata.TRIAL.traj,1),size(VRdata.TRIAL.traj,2));

    VRdata.REWARD.count(VRdata.REWARD.TRIAL>totTrials) = [];
    VRdata.REWARD.TYPE(VRdata.REWARD.TRIAL>totTrials) = [];
    VRdata.REWARD.TRIAL(VRdata.REWARD.TRIAL>totTrials) = [];

    if exist([eyeDIRname filesep eyetrackfname],'file')
        V = load([eyeDIRname filesep eyetrackfname]);
        if ~isempty(V.proc.data.pupil.com)
            if sum(~isnan(V.proc.data.pupil.com(:,2)))>0
                idx = 1:numel(V.proc.data.pupil.com(:,2));
                Xpupilpos = interp1(idx(~isnan(V.proc.data.pupil.com(:,2))), V.proc.data.pupil.com(~isnan(V.proc.data.pupil.com(:,2)),2), idx, 'linear');
                Ypupilpos = interp1(idx(~isnan(V.proc.data.pupil.com(:,1))), V.proc.data.pupil.com(~isnan(V.proc.data.pupil.com(:,1)),1), idx, 'linear');
                Ypupilpos = -(Ypupilpos - nanmean(Ypupilpos));%we want positive Y deflections when the eye goes up
                Xpupilpos = Xpupilpos - nanmean(Xpupilpos);
                pupilarea = V.proc.data.pupil.area;
                pupilarea(abs(pupilarea - nanmean(pupilarea)) > 5*std(pupilarea)) = nanmean(pupilarea); % coarse removal of blinks
                pupilarea(~isnan(pupilarea)) = medfilt1(pupilarea(~isnan(pupilarea)),5);
                meanPupilSize = mean(pupilarea);
                if ~exist([dataDIRname filesep eyelogfname])
                    S = load([eyeDIRname filesep eyelogfname]);
                    TriggerData = S.eyeLog.TriggerData;
                    udpEventTimes = S.eyeLog.udpEventTimes;
                    save([dataDIRname filesep eyelogfname],'TriggerData','udpEventTimes');
                end
                S = load([dataDIRname filesep eyelogfname]);
                frameAbstime = reshape([S.TriggerData(1:numel(S.TriggerData)).AbsTime],[numel(S.TriggerData(1).AbsTime) numel(S.TriggerData)]);
                frameAbstime = convertdate2time(frameAbstime');
                for tt = 1:totTrials
                    if 2 + (tt-1)*2 + 1 <= numel(S.udpEventTimes)
                        udpStart = convertdate2time(S.udpEventTimes{2 + (tt-1)*2 + 1});
                        udpEnd =  convertdate2time(S.udpEventTimes{2 + (tt-1)*2 + 2});

                        frameStart = max(1,find(frameAbstime >= udpStart, 1, 'first')-5);
                        frameStarttime = frameAbstime(frameStart);
                        frameEnd = min(numel(Xpupilpos),find(frameAbstime <= udpEnd, 1, 'last')+5);
                        frameEndtime = frameAbstime(frameEnd);

                        VRtime = VRdata.TRIAL.time(tt,VRdata.TRIAL.traj(tt,:) > 0 & VRdata.TRIAL.time(tt,:) > 0);
                        if ~isempty(VRtime)
                            VRtime = VRtime - VRtime(1);
                            eyetime = frameAbstime(frameStart:frameEnd);
                            eyetime = eyetime - eyetime(1);
                            eyetime = eyetime + (max(VRtime) - max(eyetime))/2;
                            VRdata.TRIAL.eyeXpos(tt,VRdata.TRIAL.traj(tt,:) > 0  & VRdata.TRIAL.time(tt,:) > 0) = interp1(eyetime, Xpupilpos(frameStart:frameEnd), VRtime', 'linear')';
                            VRdata.TRIAL.eyeYpos(tt,VRdata.TRIAL.traj(tt,:) > 0  & VRdata.TRIAL.time(tt,:) > 0) = interp1(eyetime, Ypupilpos(frameStart:frameEnd), VRtime', 'linear')';
                            VRdata.TRIAL.pupilSize(tt,VRdata.TRIAL.traj(tt,:) > 0  & VRdata.TRIAL.time(tt,:) > 0) =  interp1(eyetime, pupilarea(frameStart:frameEnd), VRtime', 'linear')';
                            VRdata.TRIAL.pupilSize(tt,VRdata.TRIAL.pupilSize(tt,:) == 0) = meanPupilSize;
                            VRdata.TRIAL.eyetime(tt,VRdata.TRIAL.traj(tt,:) > 0  & VRdata.TRIAL.time(tt,:) > 0) = VRdata.TRIAL.time(tt,VRdata.TRIAL.traj(tt,:) > 0 & VRdata.TRIAL.time(tt,:) > 0);
                            if sum(isnan(VRdata.TRIAL.eyeXpos(tt,VRdata.TRIAL.time(tt,:) > 0  & VRdata.TRIAL.time(tt,:) > 0))) > 0
                                warning('more than 5 frames missing in eye camera recording');
                            end
                        end
                    end
                end
            end
        end
    end

    VRdata.TRIAL.trialBlanks = VRdata.TRIAL.trialBlanks(1:totTrials);
    if isfield(VRdata.TRIAL,'lick')
        VRdata.TRIAL.lick  = VRdata.TRIAL.lick(1:totTrials,1:maxTrialLength);
    end
    if isfield(VRdata.TRIAL,'trialRL')
        VRdata.TRIAL.trialRL  = VRdata.TRIAL.trialRL(1:totTrials);
    end

    % Making the bad points NaNs
    VRdata.TRIAL.time(~goodTimes)       = NaN;
    VRdata.TRIAL.traj(~goodTimes)       = NaN;
    VRdata.TRIAL.trialIdx(~goodTimes)   = NaN;
    VRdata.TRIAL.time(:,1)              = NaN;
    VRdata.TRIAL.time       = VRdata.TRIAL.time - VRdata.TRIAL.time(1,2); %repmat(VRdata.TRIAL.time(:,2),1,size(VRdata.TRIAL.time(2,:),2));
    VRdata.TRIAL.currTime   = VRdata.TRIAL.time - repmat(VRdata.TRIAL.time(:,2),1,size(VRdata.TRIAL.time,2));

    VRdata.TRIAL.eyeXpos(~goodTimes)    = NaN;
    VRdata.TRIAL.eyeYpos(~goodTimes)    = NaN;
    VRdata.TRIAL.pupilSize(~goodTimes)  = NaN;
    VRdata.TRIAL.eyetime(:,1)           = NaN;
    VRdata.TRIAL.eyetime    = VRdata.TRIAL.eyetime - VRdata.TRIAL.eyetime(1,2);

    VRdata.TRIAL.reward  = NaN*ones(size(VRdata.TRIAL.traj));
    numRewards = length(VRdata.REWARD.count);

    for irew = 1:numRewards
        VRdata.TRIAL.reward(VRdata.REWARD.TRIAL(irew),VRdata.REWARD.count(irew)) = VRdata.REWARD.TYPE(irew);
    end

    for n = 1:size(goodTimes,1)
        VRdata.TRIAL.posdata(n,min(find(~goodTimes(n,:))):end,:)  = NaN;
        VRdata.TRIAL.balldata(n,min(find(~goodTimes(n,:))):end,:)  = NaN;
    end

    %%
    VRdata.traj = VRdata.TRIAL.traj;

    trajspeed = zeros(size(VRdata.traj));
    trajspeed(:,2:end) = (diff(VRdata.traj')');
    trajspeed(trajspeed > (max(VRdata.TRIAL.traj(:))/2)) = NaN;
    trajspeed(trajspeed < (-max(VRdata.TRIAL.traj(:))/2)) = NaN;
    trajspeed(:,2:end) = trajspeed(:,2:end)./(diff(VRdata.TRIAL.time')');
    % trajspeed = trajspeed

    trajspeed(trajspeed<-50) = NaN;
    trajspeed(trajspeed>200) = NaN;
    % idx = find(trajspeed);
    % trajspeed(idx(1)) = trajspeed(idx(1)) - 1000;

    balldata = VRdata.TRIAL.balldata(:,:,3);

    ballspeed = zeros(size(balldata));
    ballspeed(:,2:end) = [balldata(:,2:end)./(diff(VRdata.TRIAL.time')')];

    balldata(balldata<0) = 0;
    distTrav = zeros(size(balldata));
    totDistTrav = zeros(size(balldata));
    trajPercent = zeros(size(VRdata.traj));
    for itrial = 1:size(distTrav,1)
        distTrav(itrial,VRdata.traj(itrial,:)>0 & ~isnan(VRdata.traj(itrial,:))) = balldata(itrial,VRdata.traj(itrial,:)>0 & ~isnan(VRdata.traj(itrial,:)));
        totDistTrav(itrial,VRdata.traj(itrial,:)>=0 & ~isnan(VRdata.traj(itrial,:))) = balldata(itrial,VRdata.traj(itrial,:)>=0 & ~isnan(VRdata.traj(itrial,:)));

        trajPercent(itrial,:) = VRdata.traj(itrial,:)./(VRdata.TRIAL.trialRL(itrial));
    end

    distTrav = cumsum(distTrav,2);
    distTrav(isnan(VRdata.traj)) = nan;

    totDistTrav = cumsum(totDistTrav,2);
    totDistTrav(isnan(VRdata.traj)) = nan;

    % try
    %     VRdata.trajspeed = fastsmooth(trajspeed,10,3,1);
    %     VRdata.ballspeed = fastsmooth(ballspeed,10,3,1);
    % catch
    if isfield(VRdata.EXP, 'wheelToVR')
        VRdata.trajspeed = trajspeed./VRdata.EXP.wheelToVR;
    else
        VRdata.trajspeed = trajspeed./2;
        % Dividing by 2 because I always use gain of 2. so that should be 1.
    end
    VRdata.ballspeed = ballspeed;
    VRdata.distTrav  = distTrav;
    VRdata.totDistTrav  = totDistTrav;
    VRdata.trajPercent  = trajPercent;

    numTrials = size(VRdata.traj,1);
    % Truncating all the incomplete trial info
    VRdata.TRIAL.trialContr     = VRdata.TRIAL.trialContr(1:numTrials);
    VRdata.TRIAL.trialStart     = VRdata.TRIAL.trialStart(1:numTrials);
    VRdata.TRIAL.trialGain      = VRdata.TRIAL.trialGain(1:numTrials);
    VRdata.TRIAL.trialActive    = VRdata.TRIAL.trialActive(1:numTrials);
    VRdata.TRIAL.trialRewPos    = VRdata.TRIAL.trialRewPos(1:numTrials);
    VRdata.TRIAL.trialOutcome   = VRdata.TRIAL.trialOutcome(1:numTrials);
    VRdata.TRIAL.trialRL        = VRdata.TRIAL.trialRL(1:numTrials);

    VRdata.TRIAL.trialOutcome(VRdata.TRIAL.trialOutcome==0 & ...
        max(VRdata.TRIAL.currTime')>=VRdata.EXP.maxTrialDuration) = 5;

    if iexp > 200
        [VRdata_o] = load([VRDIRname filesep VRdata.EXP.replayed_filename]);
        %     [VRdata_o] = load([VRDIRname filesep '..' filesep '731' filesep VRdata.EXP.replayed_filename]);
        %     display('WARNING!!! reward times are wrong');
        %     VRdata.REWARD = VRdata_o.REWARD;
    else
        VRdata_o = [];
    end
    if isfield(VRdata.EXP, 'CircularMaze')
        VRdata.CircularMaze = VRdata.EXP.CircularMaze;
    else
        VRdata.CircularMaze = false;
    end
    VRdata.SmthTimeWindow = SmthTimeWindow;

    if nargout>2
        es = VREvenSample(VRdata);    
        es.iexp = zeros(size(es.sampleTimes));
        es.iexp(:) = iexp;
        es.series = ones(size(es.sampleTimes));
        es.rewardTolerance = VRdata.EXP.rew_tol;
    end
else
    VRdata = [];
    VRdata_o = [];
    es = [];
end