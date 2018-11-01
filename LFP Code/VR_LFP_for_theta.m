function [es procSpec] = VR_LFP_for_theta(animal, iseries, iexp, beh, ChnA, ChnB, sampleUP, shank_list)
% Calculate the power spectrum at the different frequencies and sample at
% the  es rate

global pepNEV

addpath(genpath('..\Behaviour analysis\'))
addpath(genpath('..\General\'))
addpath(genpath('..\Spike Analysis\'))
if nargin<8 | isempty(shank_list)
    getSpikes = 0;
else
    getSpikes = 1;
end
if nargin<4
    ChnA = 32;
    ChnB = 16;
    beh  = 0;
elseif nargin<5
    ChnB = 32;
    ChnA = 16;
elseif nargin<6
    ChnB = 16;
end
if nargin<7
    sampleUP = 1;
end

exptInfo.animal  = animal;
exptInfo.iseries = iseries;
exptInfo.iexp    = iexp;

%% Load the behaviour data
if getSpikes
    try
        es = getVRspikes(animal,iseries,iexp,1,100,1,0,shank_list,'CA1');
%         es = getVRspikes(animal,iseries,iexp,1,100,1,0,shank_list,'');
    catch
        display('Trying Ball load')
        addpath('\\zserver\Code\MouseRoom\BallTools\VRanalysis');
        es = getVRspikes(animal, iseries, iexp);
    end
else
    [~, ~ , es] = VRWheelLoad(animal, iseries, iexp);
end
runSpd = es.ballspeed(~isnan(es.ballspeed));
time   = es.sampleTimes(~isnan(es.ballspeed));
runSmthWin = 500;
smthRunSpd = smthInTime(runSpd,60,runSmthWin);

%% Load the data file
expName = [exptInfo.animal '_' num2str(exptInfo.iseries) '_' num2str(exptInfo.iexp) '.ns5'];
ns5file = ['\\ZSERVER\Data\Cerebus\' exptInfo.animal filesep num2str(exptInfo.iseries) filesep num2str(exptInfo.iexp) filesep expName]

[~,exptInfo.SamplingRateInKHZ,exptInfo.nchan] = nsopen2(ns5file);

[chan_list] = MichiganGetLayout(animal, iseries);

exptInfo.downsample = 1000;

exptInfo.ChnA = chan_list(ChnA);
exptInfo.ChnB = chan_list(ChnB);

%% Load the particular channels
timestamps.fulllength = length(pepNEV.ns.Data.data(exptInfo.ChnA,:));
timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;

% decimate(double(pepNEV.ns.Data.data(exptInfo.Chn(ichn).ChnLoc,:)), timestamps.samplingRateOrig/exptInfo.downsample);
ChnA = decimate(double(pepNEV.ns.Data.data(exptInfo.ChnA,:)), timestamps.samplingRateOrig/exptInfo.downsample);
timestamps.ChnAsampInt = 1./exptInfo.downsample;
timestamps.ChnAorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
timestamps.ChnA = 0:timestamps.ChnAsampInt:...
    (timestamps.ChnAsampInt*(timestamps.fulllength*(exptInfo.downsample/timestamps.samplingRateOrig)));

ChnB = decimate(double(pepNEV.ns.Data.data(exptInfo.ChnB,:)), timestamps.samplingRateOrig/exptInfo.downsample);
timestamps.ChnBsampInt = 1./exptInfo.downsample;
timestamps.ChnBorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
timestamps.ChnB = 0:timestamps.ChnBsampInt:...
    (timestamps.ChnBsampInt*(timestamps.fulllength*(exptInfo.downsample/timestamps.samplingRateOrig)));

timestamps.fulllength = length(pepNEV.ns.Data.data(exptInfo.ChnB,:));
timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;

procSpec.ChnA = ChnA';
procSpec.ChnB = ChnB';
procSpec.t    = timestamps.ChnB;

%% Hilbert info
procSpec.theta.frange  = [4 9];
procSpec.theta.hilbertA = LFPsignals(procSpec.ChnA', 1000, procSpec.theta.frange(1), procSpec.theta.frange(2));
procSpec.theta.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.theta.frange(1), procSpec.theta.frange(2));

%% Load the screen times on the recording
nevSamplingRateInKHZ = 30;
nevDir = ['\\ZSERVER\Data\Cerebus\' exptInfo.animal filesep num2str(exptInfo.iseries) filesep num2str(exptInfo.iexp) filesep];
period = nevSamplingRateInKHZ*1000 / nevSamplingRateInKHZ;
nevFileName = [nevDir animal '_' num2str(iseries) '_' num2str(exptInfo.iexp) '.nev'];
nevopen(nevFileName);

% screenTimes  = pepNEV.sync.timestamps;
% screenTimes(find(diff(screenTimes)<100) + 1) = [];
% procSpec.t = procSpec.t - screenTimes(1)./(1000*nevSamplingRateInKHZ);
nevclose;

global DIRS
fname = [animal '_' num2str(iseries) '_' num2str(iexp)];
dDIRname = [DIRS.multichanspikes filesep animal filesep num2str(iseries)];
load([dDIRname filesep fname '_screenTimes']);

procSpec.t = procSpec.t - screenTimes(1)./(1000*nevSamplingRateInKHZ);
timestamps.ChnA = timestamps.ChnA - screenTimes(1)./(1000*nevSamplingRateInKHZ);
timestamps.ChnB = timestamps.ChnB - screenTimes(1)./(1000*nevSamplingRateInKHZ);

% %  Correcting and matching the recording and the VR
% RecToVR_correction = screenTimes(1) - es.screenTimes2(2);
% screenTimes = screenTimes - RecToVR_correction;
% nevclose;
% procSpec.t = procSpec.t - RecToVR_correction./(1000*nevSamplingRateInKHZ);

%% Downsample the es rate

es.origLFPtime  = procSpec.t;
if sampleUP
    
    % upsampling the Hilbert data
    es.bands.theta = procSpec.theta.frange;
    
    % theta
    es.theta.A.hill    = interp1(timestamps.ChnA', procSpec.theta.hilbertA.hill', es.sampleTimes);
    es.theta.B.hill    = interp1(timestamps.ChnB', procSpec.theta.hilbertB.hill', es.sampleTimes);
    
else
    
    % down-sampling the VR data
    es.traj          = interp1q(es.sampleTimes, es.traj,      procSpec.t');
    es.trajspeed     = interp1q(es.sampleTimes, es.trajspeed,      procSpec.t');
    es.ballspeed     = interp1q(es.sampleTimes, es.ballspeed,      procSpec.t');
    es.distTrav      = interp1q(es.sampleTimes, es.distTrav,      procSpec.t');
    es.totDistTrav   = interp1q(es.sampleTimes, es.totDistTrav,      procSpec.t');
    es.trajPercent   = interp1q(es.sampleTimes, es.trajPercent,      procSpec.t');
    es.smthBallSpd   = interp1q(es.sampleTimes, es.smthBallSpd,      procSpec.t');
    es.smthTrajSpd   = interp1q(es.sampleTimes, es.smthTrajSpd,      procSpec.t');
    es.contrast      = interp1q(es.sampleTimes, es.contrast,      procSpec.t');
    es.start         = interp1q(es.sampleTimes, es.start,      procSpec.t');
    es.blanks        = interp1q(es.sampleTimes, es.blanks,      procSpec.t');
    es.active        = interp1q(es.sampleTimes, es.active,      procSpec.t');
    es.rewardPos     = interp1q(es.sampleTimes, es.rewardPos,      procSpec.t');
    es.outcome       = interp1q(es.sampleTimes, es.outcome,      procSpec.t');
    es.roomLength    = interp1q(es.sampleTimes, es.roomLength,      procSpec.t');
    
    
    lickBins = [];
    lickThing          = zeros(size(es.active));
    for ilick = find(es.lick>0 & ~isnan(es.lick))'
        lickTime = es.sampleTimes(ilick);
        newLickTime = max(find((procSpec.t - lickTime)<0));
        lickBins = [lickBins newLickTime];
        lickThing(newLickTime) = lickThing(newLickTime) + 1;
    end
    es.lick = lickThing;
    
    rew = NaN*ones(size(es.active));
    for irew = find(~isnan(es.reward))'
        rewTime = es.sampleTimes(irew);
        newRewTime = max(find((procSpec.t - rewTime)<0));
        rew(newRewTime) = es.reward(irew);
    end
    es.reward = rew;
    
    %         es.lick          = interp1q(es.sampleTimes, es.lick,      procSpec.t');
    %         es.reward        = interp1q(es.sampleTimes, es.reward,      procSpec.t');
    %
    es.trialID       = round(interp1q(es.sampleTimes, es.trialID,      procSpec.t'));
    es.gain          = interp1q(es.sampleTimes, es.gain,      procSpec.t');
    es.iexp          = round(interp1q(es.sampleTimes, es.iexp,      procSpec.t'));
    
    es.sampleRate    = 1./movingwin(2);
    es.origSampleTimes = es.sampleTimes;
    es.sampleTimes   = procSpec.t';
    
    % downsampling the Hilbert data
    es.bands.theta = procSpec.theta.frange;
    % theta
    es.A.theta.hill    = interp1q(timestamps.ChnA', procSpec.theta.hilbertA.hill', procSpec.t');
    es.B.theta.hill    = interp1q(timestamps.ChnB', procSpec.theta.hilbertB.hill', procSpec.t');
end


clear global pepNEV

