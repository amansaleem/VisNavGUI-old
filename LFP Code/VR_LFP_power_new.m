function [es procSpec] = VR_LFP_power_new...
    (animal, iseries, iexp, chan2consider, beh, movingwin, internalRef)
% Calculate the power spectrum at the different frequencies and sample at
% the  es rate

global pepNEV
SetDefaultDirs

if nargin<7
    internalRef = [];
end

addpath(genpath('..\Behaviour analysis\'))
addpath(genpath('..\General\'))
addpath(genpath('..\Spike Analysis\'))

if nargin<5
    beh = 0;
end

if nargin<6
    movingwin=[3 1];
end
% beh = 0;
sampleUP = 0;

exptInfo.animal  = animal;
exptInfo.iseries = iseries;
exptInfo.iexp    = iexp;

%% Load the behaviour data
if beh
    try
        [~, ~ , es] = VRWheelLoad(animal, iseries, iexp);
    catch
        display('Trying Ball load')
        addpath('\\zserver\Code\MouseRoom\BallTools\VRanalysis');
        es = getVRspikeTimes(animal, iseries, iexp);
    end
    runSpd = es.ballspeed(~isnan(es.ballspeed));
    time   = es.sampleTimes(~isnan(es.ballspeed));
    runSmthWin = 500;
    smthRunSpd = smthInTime(runSpd,60,runSmthWin);
end

%% Load the data file
expName = [exptInfo.animal '_' num2str(exptInfo.iseries) '_' num2str(exptInfo.iexp) '.ns5'];
ns5file = ['\\ZSERVER\Data\Cerebus\' exptInfo.animal filesep num2str(exptInfo.iseries) filesep num2str(exptInfo.iexp) filesep expName]

% try
    [~,exptInfo.SamplingRateInKHZ,exptInfo.nchan] = nsopen2(ns5file);
% catch
%     exptInfo.SamplingRateInKHZ = 30000;
%     exptInfo.SamplingRateInKHZ,exptInfo.nchan = 32;
% end

[chan_list] = MichiganGetLayout(animal, iseries);

numChn = length(chan2consider);

ChnC = chan_list(end);

exptInfo.downsample = 1000;
exptInfo.ChnC = chan_list(end);

ichn = 1;
exptInfo.Chn(ichn).Chn = chan2consider(ichn);
exptInfo.Chn(ichn).ChnLoc = chan_list(chan2consider(ichn));

%% Load the particular channels
timestamps.fulllength = length(pepNEV.ns.Data.data(chan_list(exptInfo.Chn(1).Chn),:));
timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;
timestamps.ChnsampInt = 1./exptInfo.downsample;
timestamps.Chnorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
timestamps.Chn = 0:timestamps.ChnsampInt:...
    (timestamps.ChnsampInt*(timestamps.fulllength*(exptInfo.downsample/timestamps.samplingRateOrig)));

    if ~isempty(internalRef)
        exptInfo.Chn(ichn).data = decimate(double(pepNEV.ns.Data.data(exptInfo.Chn(ichn).ChnLoc,:) - pepNEV.ns.Data.data(chan_list(internalRef),:)), timestamps.samplingRateOrig/exptInfo.downsample);
    else
        exptInfo.Chn(ichn).data = decimate(double(pepNEV.ns.Data.data(exptInfo.Chn(ichn).ChnLoc,:)), timestamps.samplingRateOrig/exptInfo.downsample);
    end

% Getting run speed from the last channel
exptInfo.runSmthWin = 300;
ChnC = zeros(1, size(pepNEV.ns.Data.data,2));
maxC = max(pepNEV.ns.Data.data(exptInfo.ChnC,:));
minC = min(pepNEV.ns.Data.data(exptInfo.ChnC,:));
ChnC((abs(diff(pepNEV.ns.Data.data(exptInfo.ChnC,:)-minC)./(maxC-minC)))>0.5) = 1;
ChnC_down = sum(reshape(ChnC(1:end-rem(size(pepNEV.ns.Data.data,2),250)),250,[]),1);

timestamps.ChnCsampInt = 1./exptInfo.downsample;
timestamps.ChnCorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
timestamps.ChnC = 0:timestamps.ChnsampInt:...
    (timestamps.ChnCsampInt*(timestamps.fulllength*(exptInfo.downsample/timestamps.samplingRateOrig)));
if length(timestamps.ChnC)>length(ChnC)
    timestamps.ChnC = timestamps.ChnC(1:length(ChnC));
end
timestamps.ChnCsampInt_down = 1./(timestamps.samplingRateOrig./250);
timestamps.ChnC_down = 0:timestamps.ChnCsampInt_down:...
    (timestamps.ChnCsampInt_down*(timestamps.fulllength*(1/250)));
if length(timestamps.ChnC_down)>length(ChnC_down)
    timestamps.ChnC_down = timestamps.ChnC_down(1:length(ChnC_down));
end
exptInfo.runSpd = ChnC_down;
exptInfo.smthRunSpd = smthInTime(exptInfo.runSpd, 1./timestamps.ChnCsampInt_down, exptInfo.runSmthWin);
exptInfo.time = timestamps.ChnC_down;

%% Get the spectrogram

params.Fs=exptInfo.downsample;
exptInfo.range_low = 0.5;
exptInfo.range_high = 250;
params.fpass=[exptInfo.range_low exptInfo.range_high];
params.tapers=[5 9];

[procSpec.SA,procSpec.t,procSpec.f]=...
            mtspecgramc(exptInfo.Chn(ichn).data',movingwin,params);
params.movingwin = movingwin;

%% Load the screen times on the recording
if beh
    nevSamplingRateInKHZ = 30;
    nevDir = ['\\ZSERVER\Data\Cerebus\' exptInfo.animal filesep num2str(exptInfo.iseries) filesep num2str(exptInfo.iexp) filesep];
    period = nevSamplingRateInKHZ*1000 / nevSamplingRateInKHZ;
    nevFileName = [nevDir animal '_' num2str(iseries) '_' num2str(exptInfo.iexp) '.nev'];
    nevopen(nevFileName);
    
    screenTimes  = pepNEV.sync.timestamps;
    
    screenTimes(find(diff(screenTimes)<100) + 1) = [];
    nevclose;
    procSpec.t = procSpec.t - screenTimes(1)./(1000*nevSamplingRateInKHZ);
end

%% Downsample the es rate

es.freq         = procSpec.f;
es.origLFPtime  = procSpec.t;
es.params = params;
if beh
    
    if sampleUP
        
        es.powA     = interp1q(procSpec.t', procSpec.SA, es.sampleTimes); %procSpec.SA;
        
    else
        es.powA     = procSpec.SA;
        
        % down-sampling the VR data
        t = procSpec.t;
        es.t = t;
        es.freq = procSpec.f;
        es.traj          = interp1q(es.sampleTimes, es.traj,      t');
        es.trajspeed     = interp1q(es.sampleTimes, es.trajspeed,      t');
        es.ballspeed     = interp1q(es.sampleTimes, es.ballspeed,      t');
        es.distTrav      = interp1q(es.sampleTimes, es.distTrav,      t');
        es.totDistTrav   = interp1q(es.sampleTimes, es.totDistTrav,      t');
        es.trajPercent   = interp1q(es.sampleTimes, es.trajPercent,      t');
        es.smthBallSpd   = interp1q(es.sampleTimes, es.smthBallSpd,      t');
        es.smthTrajSpd   = interp1q(es.sampleTimes, es.smthTrajSpd,      t');
        es.contrast      = interp1q(es.sampleTimes, es.contrast,      t');
        es.start         = interp1q(es.sampleTimes, es.start,      t');
        es.blanks        = interp1q(es.sampleTimes, es.blanks,      t');
        es.active        = interp1q(es.sampleTimes, es.active,      t');
        es.rewardPos     = interp1q(es.sampleTimes, es.rewardPos,      t');
        es.outcome       = interp1q(es.sampleTimes, es.outcome,      t');
        es.roomLength    = interp1q(es.sampleTimes, es.roomLength,      t');
        es.lick          = interp1q(es.sampleTimes, es.lick,      t');
        es.reward        = interp1q(es.sampleTimes, es.reward,      t');
        es.trialID       = round(interp1q(es.sampleTimes, es.trialID,      t'));
        es.gain          = interp1q(es.sampleTimes, es.gain,      t');
        es.iexp          = round(interp1q(es.sampleTimes, es.iexp,      t'));
        
        es.sampleRate    = 1./movingwin(2);
        es.origSampleTimes = es.sampleTimes;
        es.sampleTimes   = t';
        
    end
else
    es.powA     = squeeze(procSpec.SA);
    es.freq     = squeeze(procSpec.f);
    es.t        = squeeze(procSpec.t);
    
    es.smthRunSpd = interp1q(exptInfo.time',exptInfo.smthRunSpd', procSpec.t');
    
    clear procSpec
    procSpec.powA       = permute(squeeze(nanmean(es.powA,3)),[3 1 2]);
    procSpec.f          = es.freq;
    procSpec.t          = es.t; 
end
clear global pepNEV