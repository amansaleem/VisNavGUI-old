function es = VREvenSample(VRdata, sampleTimes, screenTimes,samplerate)


[es.screenTimes2 es.sortIDX] = sort(VRdata.TRIAL.time(:));
maxTime = min(find(isnan(es.screenTimes2)))-1;
es.screenTimes2 = es.screenTimes2(1:maxTime);
if nargin < 4
    samplerate = 60;
end
% es.sampleRate = '60 Hz';
es.sampleRate = samplerate;%60;

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
        if isempty(maxTime)
            maxTime = numel(es.screenTimes);
        end
        es.screenTimes = es.screenTimes(1:maxTime);
    elseif numel(es.screenTimes) > numel(es.screenTimes2)
        display('WARNING!!!! TOO MANY SCREENTIMES IN THE RECORDING: will check that they are at the end')
        maxTime  = min(find((es.screenTimes-max(es.screenTimes2))>0));
        es.screenTimes = es.screenTimes(1:maxTime);
    end    
end
[vmax,imax] = max(xcorr(es.screenTimes2(1:min(numel(es.screenTimes2),numel(es.screenTimes))),es.screenTimes(1:min(numel(es.screenTimes2),numel(es.screenTimes))),30,'coeff'));
if vmax < 0.999 || abs(imax - 31) > 2
    warning('additional screen times may be not at the end');
end

evenSampleTime = (1/es.sampleRate):(1/es.sampleRate):max(es.screenTimes);
es.sampleTimes = evenSampleTime';

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

for i = 1:size(VRdata.TRIAL.trialIdx,1)
    VRdata.TRIAL.trialIdx(i,VRdata.TRIAL.trialIdx(i,:) == 0) = max(VRdata.TRIAL.trialIdx(i,:));
end

traj      = VRdata.traj(es.sortIDX);
currTime  = VRdata.TRIAL.currTime(es.sortIDX);
trialIdx  = VRdata.TRIAL.trialIdx(es.sortIDX);
lick      = VRdata.TRIAL.lick(es.sortIDX);
reward    = VRdata.TRIAL.reward(es.sortIDX);

eyeXpos   = VRdata.TRIAL.eyeXpos(es.sortIDX);
eyeYpos   = VRdata.TRIAL.eyeYpos(es.sortIDX);
pupilSize   = VRdata.TRIAL.pupilSize(es.sortIDX);

%JUL - 15.09.2015: to deal with cases where 2 sucessive indexes have 
%a trialID = 0 value.
%trialIdx(find(trialIdx==0))= trialIdx(find(trialIdx==0)-1);

% trialIdx(find(trialIdx==0))= max(trialIdx(find(trialIdx==0)-1),trialIdx(find(trialIdx==0)-2));


lickTimes = es.screenTimes(lick>0);
%JUL 15.07.16: modified to avoid repmat related memory issues when number of
%samples x licktimes gets too large
% idx = zeros(1,numel(lickTimes));
% [~,idx] = (min(abs((repmat(es.sampleTimes,1,length(lickTimes)))...
% - (repmat(lickTimes',size(es.sampleTimes,1),1)))));
idx = zeros(1,numel(lickTimes));
Tsample = min(diff(es.sampleTimes));
for l = 1:numel(lickTimes)
    idx(l) = find((es.sampleTimes - lickTimes(l)) > -Tsample,1,'first');
end
es.lick      = zeros(size(es.sampleTimes));
es.lick(idx) = 1;

% [~,idx] = (min(abs((repmat(es.sampleTimes,1,length(lickTimes)))...
% - (repmat(lickTimes',size(es.sampleTimes,1),1)))));
% es.lick      = zeros(size(es.sampleTimes));
% es.lick(idx) = 1;


rewTimes = es.screenTimes(~isnan(reward));
[~,idx] = (min(abs((repmat(es.sampleTimes,1,length(rewTimes)))...
- (repmat(rewTimes',size(es.sampleTimes,1),1)))));
es.reward      = NaN*ones(size(es.sampleTimes));
try
    es.reward(idx) = VRdata.REWARD.TYPE;
catch
    disp('Rewards are wrong')
    try
        es.reward(idx) = VRdata.REWARD.TYPE(1:numel(rewTimes));
    catch
        rewTimes = es.screenTimes(~isnan(reward) & reward == 2);
        es.reward(idx) = VRdata.REWARD.TYPE(1:numel(rewTimes));
    end
end

es.trialID       = round(interp1(es.screenTimes, trialIdx, es.sampleTimes, 'nearest'));
%JUL: changed interpolation method on traj and trajpercent to avoid
%articfact due to missing frame at the beginning of every trial
% es.traj          = interp1(es.screenTimes, traj,es.sampleTimes, 'nearest');
es.traj          = mod((interp1(es.screenTimes, unwrap(traj/(max(traj))*2*pi)*(max(traj))/(2*pi),es.sampleTimes, 'linear')),(max(traj)));
% es.traj          = mod((interp1(es.screenTimes, unwrap(traj/round(max(traj))*2*pi)*round(max(traj))/(2*pi),es.sampleTimes, 'linear')),round(max(traj)));
es.trajspeed     = interp1(es.screenTimes(~isnan(trajspeed)), trajspeed(~isnan(trajspeed)), es.sampleTimes, 'nearest');
es.ballspeed     = interp1(es.screenTimes(~isnan(ballspeed)), ballspeed(~isnan(ballspeed)), es.sampleTimes, 'nearest');
es.currTime      = interp1(es.screenTimes, currTime,  es.sampleTimes, 'linear');

es.distTrav      = interp1(es.screenTimes, distTrav,  es.sampleTimes, 'nearest');
es.totDistTrav      = interp1(es.screenTimes, totDistTrav,  es.sampleTimes, 'nearest');
es.trajPercent      = interp1(es.screenTimes, trajPercent,  es.sampleTimes, 'nearest');
% es.trajPercent          = mod((interp1(es.screenTimes, unwrap(trajPercent/round(max(trajPercent))*2*pi)*round(max(trajPercent))/(2*pi),es.sampleTimes, 'linear')),round(max(trajPercent)));

es.smthBallSpd   = NaN*ones(size(es.traj));
es.smthBallSpd(~isnan(es.ballspeed))   = smthInTime(es.ballspeed(~isnan(es.ballspeed)), es.sampleRate, VRdata.SmthTimeWindow, [], [], 'median');

es.smthBallAcc   = NaN*ones(size(es.traj));
es.smthBallAcc(~isnan(es.smthBallSpd))   = [diff(es.smthBallSpd(~isnan(es.smthBallSpd))) ; 0] ;
es.smthBallAcc(~isnan(es.smthBallAcc))   = smthInTime(es.smthBallAcc(~isnan(es.smthBallAcc)), es.sampleRate, VRdata.SmthTimeWindow);

es.smthTrajSpd   = NaN*ones(size(es.traj));
es.smthTrajSpd(~isnan(es.trajspeed))   = smthInTime(es.trajspeed(~isnan(es.trajspeed)), es.sampleRate, VRdata.SmthTimeWindow, [], [], 'median');

es.eyeXpos = interp1(es.screenTimes, eyeXpos, es.sampleTimes, 'linear');
es.eyeYpos = interp1(es.screenTimes, eyeYpos, es.sampleTimes, 'linear');
es.pupilSize = interp1(es.screenTimes, pupilSize, es.sampleTimes, 'linear');

es.contrast = NaN*ones(size(es.traj));
es.start    = NaN*ones(size(es.traj));
es.gain     = NaN*ones(size(es.traj));
es.blanks   = NaN*ones(size(es.traj));
es.active   = NaN*ones(size(es.traj));
es.rewardPos= NaN*ones(size(es.traj));
es.outcome  = NaN*ones(size(es.traj));
es.roomLength= NaN*ones(size(es.traj));

trialIDs = unique(es.trialID);
numTrials = length(trialIDs);

for itrial = 1:numTrials
    es.contrast(es.trialID==itrial)    = VRdata.TRIAL.trialContr(trialIDs(itrial));
    es.start(es.trialID==itrial)       = VRdata.TRIAL.trialStart(trialIDs(itrial));
    es.gain(es.trialID==itrial)        = VRdata.TRIAL.trialGain(trialIDs(itrial));
    es.blanks(es.trialID==itrial)      = VRdata.TRIAL.trialBlanks(trialIDs(itrial));
    es.active(es.trialID==itrial)      = VRdata.TRIAL.trialActive(trialIDs(itrial));
    es.rewardPos(es.trialID==itrial)   = VRdata.TRIAL.trialRewPos(trialIDs(itrial));
    es.outcome(es.trialID==itrial)     = VRdata.TRIAL.trialOutcome(trialIDs(itrial));
    es.roomLength(es.trialID==itrial)  = VRdata.TRIAL.trialRL(trialIDs(itrial));
end

if nargin==3
    display('WARNING!!!! temporary fix. Check the photodiode signalling');
    evenSampleTime = evenSampleTime + screenTimes(1);
end

if VRdata.CircularMaze
    wrapfactor = 2;%1;% 
    disp(['wrapping up circular maze by a factor of ' num2str(wrapfactor)]);    
    es = VRwrap(es,wrapfactor,true);
    es = VRredefineOutcome(es);
else
    es = redefineOutcome_AS(es);
end

es.CircularMaze = VRdata.CircularMaze;