function es = VREvenSample_SL(VRdata, sampleTimes, screenTimes)


[es.screenTimes2 es.sortIDX] = sort(VRdata.TRIAL.time(:));
maxTime = min(find(isnan(es.screenTimes2)))-1;
es.screenTimes2 = es.screenTimes2(1:maxTime);

if nargin<3
    es.screenTimes = es.screenTimes2;
else
    evenSampleTime = sampleTimes;
    es.screenTimes = screenTimes;
    display('WARNING!!!! temporary fix. Check the photodiode signalling');
    es.screenTimes = es.screenTimes'-screenTimes(1);
    %condition added to deal with cases where the recording has more
    %screentimes than VR at the end. Please check that additional
    %screentimes are indeed at the end of the recording and not anywhere
    %else.
    if numel(es.screenTimes) < numel(es.screenTimes2)
        maxTime  = min(find((es.screenTimes2-max(es.screenTimes))>0));
        es.screenTimes = es.screenTimes(1:maxTime);
    else
        display('WARNING!!!! TOO MANY SCREENTIMES IN THE RECORDING: please check that they are at the end')
        maxTime  = min(find((es.screenTimes-max(es.screenTimes2))>0));
        es.screenTimes = es.screenTimes(1:maxTime);
    end
end

evenSampleTime = (1/60):(1/60):max(es.screenTimes);


es.sortIDX = es.sortIDX(1:maxTime);

ballspeed = VRdata.ballspeed(es.sortIDX);
ballspeed(1) = 0;

trajspeed = VRdata.trajspeed(es.sortIDX);
trajspeed(1) = 0;

distTrav  = VRdata.distTrav(es.sortIDX);
distTrav(1)  = 0;

totDistTrav = VRdata.totDistTrav(es.sortIDX);
totDistTrav(1)  = 0;

trajPercent = VRdata.trajPercent(es.sortIDX);
trajPercent(1) = 0;

traj      = VRdata.traj(es.sortIDX);
currTime  = VRdata.TRIAL.currTime(es.sortIDX);
trialIdx  = VRdata.TRIAL.trialIdx(es.sortIDX);
lick      = VRdata.TRIAL.lick(es.sortIDX);
reward    = VRdata.TRIAL.reward(es.sortIDX);

SyncPulse = VRdata.TRIAL.balldata(:,:,4);
SyncPulse = SyncPulse(es.sortIDX);

%JUL - 15.09.2015: to deal with cases where 2 sucessive indexes have
%a trialID = 0 value.
%trialIdx(find(trialIdx==0))= trialIdx(find(trialIdx==0)-1);
trialIdx(find(trialIdx==0))= max(trialIdx(find(trialIdx==0)-1),trialIdx(find(trialIdx==0)-2));

es.sampleRate = '60 Hz';
es.sampleSize = 60;
es.sampleTimes = evenSampleTime';

lickTimes = es.screenTimes(lick>0);
[~,idx] = (min(abs((repmat(es.sampleTimes,1,length(lickTimes)))...
    - (repmat(lickTimes',size(es.sampleTimes,1),1)))));
es.lick      = zeros(size(es.sampleTimes));
es.lick(idx) = 1;


rewTimes = es.screenTimes(~isnan(reward));
[~,idx] = (min(abs((repmat(es.sampleTimes,1,length(rewTimes)))...
    - (repmat(rewTimes',size(es.sampleTimes,1),1)))));
es.reward      = NaN*ones(size(es.sampleTimes));
try
    es.reward(idx) = VRdata.REWARD.TYPE;
catch
    display('Rewards are wrong')
end

es.trialID       = round(interp1(es.screenTimes, trialIdx, es.sampleTimes));
es.traj          = interp1(es.screenTimes, traj,      es.sampleTimes);
es.trajspeed     = interp1(es.screenTimes, trajspeed, es.sampleTimes);
es.ballspeed     = interp1(es.screenTimes, ballspeed, es.sampleTimes);
es.currTime      = interp1(es.screenTimes, currTime,  es.sampleTimes);

es.distTrav      = interp1(es.screenTimes, distTrav,  es.sampleTimes);
es.totDistTrav      = interp1(es.screenTimes, totDistTrav,  es.sampleTimes);
es.trajPercent      = interp1(es.screenTimes, trajPercent,  es.sampleTimes);

es.SyncPulse     = interp1(es.screenTimes, SyncPulse,  es.sampleTimes); % syncpulse for ephys recordings

es.smthBallSpd   = NaN*ones(size(es.traj));
% es.smthBallSpd(~isnan(es.ballspeed))   = smthInTime(es.ballspeed(~isnan(es.ballspeed)), 60, 500);

es.smthTrajSpd   = NaN*ones(size(es.traj));
% es.smthTrajSpd(~isnan(es.trajspeed))   = smthInTime(es.trajspeed(~isnan(es.trajspeed)), 60, 500);

es.contrast = NaN*ones(size(es.traj));
es.start    = NaN*ones(size(es.traj));
es.gain     = NaN*ones(size(es.traj));
es.blanks   = NaN*ones(size(es.traj));
es.active   = NaN*ones(size(es.traj));
es.rewardPos= NaN*ones(size(es.traj));
es.outcome  = NaN*ones(size(es.traj));
es.roomLength= NaN*ones(size(es.traj));

if isfield(VRdata.TRIAL, 'tex1pos')
    es.tex1pos= NaN*ones(size(es.traj));
    es.tex2pos= NaN*ones(size(es.traj));
    es.tex3pos= NaN*ones(size(es.traj));
    es.tex4pos= NaN*ones(size(es.traj));
    es.waveLength= NaN*ones(size(es.traj));
    es.currList= NaN*ones(size(es.traj));
end

trialIDs = unique(es.trialID);
numTrials = length(trialIDs);

for itrial = 1:numTrials
    es.contrast(es.trialID==trialIDs(itrial) & es.traj~=0)    = VRdata.TRIAL.trialContr(trialIDs(itrial));
    es.start(es.trialID==trialIDs(itrial) & es.traj~=0)       = VRdata.TRIAL.trialStart(trialIDs(itrial));
    es.gain(es.trialID==trialIDs(itrial) & es.traj~=0)        = VRdata.TRIAL.trialGain(trialIDs(itrial));
    es.blanks(es.trialID==trialIDs(itrial) & es.traj~=0)      = VRdata.TRIAL.trialBlanks(trialIDs(itrial));
    es.active(es.trialID==trialIDs(itrial) & es.traj~=0)      = VRdata.TRIAL.trialActive(trialIDs(itrial));
    es.rewardPos(es.trialID==trialIDs(itrial) & es.traj~=0)   = VRdata.TRIAL.trialRewPos(trialIDs(itrial));
    es.outcome(es.trialID==trialIDs(itrial) & es.traj~=0)     = VRdata.TRIAL.trialOutcome(trialIDs(itrial));
    es.roomLength(es.trialID==trialIDs(itrial) & es.traj~=0)  = VRdata.TRIAL.trialRL(trialIDs(itrial));
    if isfield(VRdata.TRIAL, 'tex1pos')
        es.tex1pos(es.trialID==trialIDs(itrial) & es.traj~=0)  = VRdata.TRIAL.tex1pos(trialIDs(itrial));
        es.tex2pos(es.trialID==trialIDs(itrial) & es.traj~=0)  = VRdata.TRIAL.tex2pos(trialIDs(itrial));
        es.tex3pos(es.trialID==trialIDs(itrial) & es.traj~=0)  = VRdata.TRIAL.tex3pos(trialIDs(itrial));
        es.tex4pos(es.trialID==trialIDs(itrial) & es.traj~=0)  = VRdata.TRIAL.tex4pos(trialIDs(itrial));
        es.waveLength(es.trialID==trialIDs(itrial) & es.traj~=0)  = VRdata.TRIAL.waveLength(trialIDs(itrial));
        es.currList(es.trialID==trialIDs(itrial) & es.traj~=0)  = VRdata.TRIAL.currList(trialIDs(itrial));
    end
end


if nargin==3
    display('WARNING!!!! temporary fix. Check the photodiode signalling');
    evenSampleTime = evenSampleTime + screenTimes(1);
end