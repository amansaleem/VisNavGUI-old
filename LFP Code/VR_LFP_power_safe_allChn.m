function [es procSpec] = VR_LFP_power_safe_allChn...
    (animal, iseries, iexp, chan2consider, beh, movingwin, internalRef)
% Calculate the power spectrum at the different frequencies and sample at
% the  es rate

global pepNEV

if nargin<7
    internalRef = [];
end

addpath(genpath('..\Behaviour analysis\'))
addpath(genpath('..\General\'))
addpath(genpath('..\Spike Analysis\'))

if nargin<4
    getChn = 1;
else
    getChn = 0;
end

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

[~,exptInfo.SamplingRateInKHZ,exptInfo.nchan] = nsopen2(ns5file);

[chan_list] = MichiganGetLayout(animal, iseries);

if getChn
    numChn = length(numChn);
else
    numChn = length(chan2consider);
end

ChnC = chan_list(end);

exptInfo.downsample = 1000;
exptInfo.ChnC = chan_list(end);

for ichn = 1:numChn
   exptInfo.Chn(ichn).Chn = chan2consider(ichn);
   exptInfo.Chn(ichn).ChnLoc = chan_list(chan2consider(ichn));
end

%% Load the particular channels
timestamps.fulllength = length(pepNEV.ns.Data.data(chan_list(exptInfo.Chn(1).Chn),:));
timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;
timestamps.ChnsampInt = 1./exptInfo.downsample;
timestamps.Chnorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
timestamps.Chn = 0:timestamps.ChnsampInt:...
    (timestamps.ChnsampInt*(timestamps.fulllength*(exptInfo.downsample/timestamps.samplingRateOrig)));

for ichn = 1:numChn
    if ~isempty(internalRef)
        exptInfo.Chn(ichn).data = decimate(double(pepNEV.ns.Data.data(exptInfo.Chn(ichn).ChnLoc,:) - pepNEV.ns.Data.data(chan_list(internalRef),:)), timestamps.samplingRateOrig/exptInfo.downsample);
    else
        exptInfo.Chn(ichn).data = decimate(double(pepNEV.ns.Data.data(exptInfo.Chn(ichn).ChnLoc,:)), timestamps.samplingRateOrig/exptInfo.downsample);
    end
    %     exptInfo.Chn(ichn).data = MKnotch(exptInfo.Chn(ichn).data, 4, exptInfo.downsample, 50);
end

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

for chn1 = 1%:numChn
    for chn2 = 1:numChn
        [procSpec.C(chn1,chn2,:,:),procSpec.phi(chn1,chn2,:,:),procSpec.SAB(chn1,chn2,:,:),procSpec.SA(chn1,chn2,:,:),procSpec.SB(chn1,chn2,:,:),procSpec.t(chn1,chn2,:),procSpec.f(chn1,chn2,:)]=...
            cohgramc(exptInfo.Chn(chn1).data',exptInfo.Chn(chn2).data',movingwin,params);
    end
end
params.movingwin = movingwin;
if numChn<3
    procSpec.SA = squeeze(procSpec.SA(1:chn1,1,:,:));
    procSpec.SB = squeeze(procSpec.SB(1,1:chn2,:,:));
else
    procSpec.SA = squeeze(procSpec.SA);
    procSpec.SB = squeeze(procSpec.SB);
end

% for chn1 = 1%:numChn
%     for chn2 = 1:numChn
%         procSpec.C(chn2,chn1,:,:) = procSpec.C(chn1,chn2,:,:);
%         procSpec.phi(chn2,chn1,:,:) = procSpec.phi(chn1,chn2,:,:);
%         procSpec.SAB(chn2,chn1,:,:) = procSpec.SAB(chn1,chn2,:,:);
%         procSpec.SA(chn2,chn1,:,:) = procSpec.SA(chn1,chn2,:,:);
%         procSpec.SB(chn2,chn1,:,:) = procSpec.SB(chn1,chn2,:,:);
%         procSpec.f(chn2,chn1,:,:) = procSpec.f(chn1,chn2,:);
%         procSpec.t(chn2,chn1,:,:) = procSpec.t(chn1,chn2,:);
%     end
% end
% %% Hilbert info
% procSpec.delta.frange  = [1 4.5];
% procSpec.delta.hilbertA = LFPsignals(procSpec.ChnA', 1000, procSpec.delta.frange(1), procSpec.delta.frange(2));
% procSpec.delta.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.delta.frange(1), procSpec.delta.frange(2));
% procSpec.delta.hilbertC = LFPsignals(procSpec.ChnC', 1000, procSpec.delta.frange(1), procSpec.delta.frange(2));
% 
% procSpec.theta.frange  = [4.5 9];
% procSpec.theta.hilbertA = LFPsignals(procSpec.ChnA', 1000, procSpec.theta.frange(1), procSpec.theta.frange(2));
% procSpec.theta.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.theta.frange(1), procSpec.theta.frange(2));
% procSpec.theta.hilbertC = LFPsignals(procSpec.ChnC', 1000, procSpec.theta.frange(1), procSpec.theta.frange(2));
% 
% procSpec.beta.frange  = [10 30];
% procSpec.beta.hilbertA = LFPsignals(procSpec.ChnA', 1000, procSpec.beta.frange(1), procSpec.beta.frange(2));
% procSpec.beta.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.beta.frange(1), procSpec.beta.frange(2));
% 
% procSpec.gamma.frange  = [52 95];
% procSpec.gamma.hilbertA = LFPsignals(procSpec.ChnA', 1000, procSpec.gamma.frange(1), procSpec.gamma.frange(2));
% procSpec.gamma.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.gamma.frange(1), procSpec.gamma.frange(2));
% 
% procSpec.gamma_narrow.frange  = [58 76];
% procSpec.gamma_narrow.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.gamma_narrow.frange(1), procSpec.gamma_narrow.frange(2));

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
        
        es.coherence= interp1q(procSpec.t', procSpec.C, es.sampleTimes);
        es.cohPhi   = interp1q(procSpec.t', procSpec.phi, es.sampleTimes); %procSpec.phi;
        es.powA     = interp1q(procSpec.t', procSpec.SA, es.sampleTimes); %procSpec.SA;
        es.powB     = interp1q(procSpec.t', procSpec.SB, es.sampleTimes); %procSpec.SB;
        es.powAB    = interp1q(procSpec.t', procSpec.SAB, es.sampleTimes); %procSpec.SAB;
        
        % downsampling the Hilbert data
        es.bands.delta = procSpec.delta.frange;
        es.bands.theta = procSpec.theta.frange;
        es.bands.beta = procSpec.beta.frange;
        es.bands.gamma = procSpec.gamma.frange;
        es.bands.gamma_narrow = procSpec.gamma_narrow.frange;
        
        % delta
        es.A.delta.hill    = interp1q(timestamps.ChnA', procSpec.delta.hilbertA.hill', es.sampleTimes);
        es.B.delta.hill    = interp1q(timestamps.ChnB', procSpec.delta.hilbertB.hill', es.sampleTimes);
        es.C.delta.hill    = interp1q(timestamps.ChnC_down', procSpec.delta.hilbertC.hill', es.sampleTimes);
        % theta
        es.A.theta.hill    = interp1q(timestamps.ChnA', procSpec.theta.hilbertA.hill', es.sampleTimes);
        es.B.theta.hill    = interp1q(timestamps.ChnB', procSpec.theta.hilbertB.hill', es.sampleTimes);
        es.C.theta.hill    = interp1q(timestamps.ChnC_down', procSpec.theta.hilbertC.hill', es.sampleTimes);
        % beta
        es.A.beta.hill    = interp1q(timestamps.ChnA', procSpec.beta.hilbertA.hill', es.sampleTimes);
        es.B.beta.hill    = interp1q(timestamps.ChnB', procSpec.beta.hilbertB.hill', es.sampleTimes);
        % gamma
        es.A.gamma.hill    = interp1q(timestamps.ChnA', procSpec.gamma.hilbertA.hill', es.sampleTimes);
        es.B.gamma.hill    = interp1q(timestamps.ChnB', procSpec.gamma.hilbertB.hill', es.sampleTimes);
        % gamma_narrow
        es.B.gamma_narrow.hill    = interp1q(timestamps.ChnB', procSpec.gamma_narrow.hilbertB.hill', es.sampleTimes);
        
    else
        es.coherence= procSpec.C;
        es.cohPhi   = procSpec.phi;
        es.powA     = procSpec.SA;
        es.powB     = procSpec.SB;
        es.powAB    = procSpec.SAB;
        
        % down-sampling the VR data
        t = squeeze(procSpec.t(1,1,:))';
        es.t = t;
        es.freq = squeeze(procSpec.f(1,1,:));
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
        
%         % downsampling the Hilbert data
%         es.bands.delta = procSpec.delta.frange;
%         es.bands.theta = procSpec.theta.frange;
%         es.bands.beta = procSpec.beta.frange;
%         es.bands.gamma = procSpec.gamma.frange;
%         es.bands.gamma_narrow = procSpec.gamma_narrow.frange;
%         
%         % delta
%         es.A.delta.hill    = interp1q(timestamps.ChnA', procSpec.delta.hilbertA.hill', procSpec.t');
%         es.B.delta.hill    = interp1q(timestamps.ChnB', procSpec.delta.hilbertB.hill', procSpec.t');
%         es.C.delta.hill    = interp1q(timestamps.ChnC_down', procSpec.delta.hilbertC.hill', procSpec.t');
%         % theta
%         es.A.theta.hill    = interp1q(timestamps.ChnA', procSpec.theta.hilbertA.hill', procSpec.t');
%         es.B.theta.hill    = interp1q(timestamps.ChnB', procSpec.theta.hilbertB.hill', procSpec.t');
%         es.C.theta.hill    = interp1q(timestamps.ChnC_down', procSpec.theta.hilbertC.hill', procSpec.t');
%         % beta
%         es.A.beta.hill    = interp1q(timestamps.ChnA', procSpec.beta.hilbertA.hill', procSpec.t');
%         es.B.beta.hill    = interp1q(timestamps.ChnB', procSpec.beta.hilbertB.hill', procSpec.t');
%         % gamma
%         es.A.gamma.hill    = interp1q(timestamps.ChnA', procSpec.gamma.hilbertA.hill', procSpec.t');
%         es.B.gamma.hill    = interp1q(timestamps.ChnB', procSpec.gamma.hilbertB.hill', procSpec.t');
%         % gamma_narrow
%         es.B.gamma_narrow.hill    = interp1q(timestamps.ChnB', procSpec.gamma_narrow.hilbertB.hill', procSpec.t');
        
    end
else
    es.coherence= squeeze(procSpec.C);
    es.cohPhi   = squeeze(procSpec.phi);
    es.powA     = squeeze(procSpec.SA);
    es.powB     = squeeze(procSpec.SB);
    es.powAB    = squeeze(procSpec.SAB);
    es.freq     = squeeze(procSpec.f(1,1,:));
    es.t        = squeeze(procSpec.t(1,1,:));
    
    es.smthRunSpd = interp1q(exptInfo.time',exptInfo.smthRunSpd', squeeze(procSpec.t(1,1,:)));
    
    clear procSpec
    procSpec.coherence  = permute(squeeze(nanmean(es.coherence,3)),[3 1 2]);
    procSpec.phi        = permute(squeeze(nanmean(es.cohPhi,3)),[3 1 2]);
    procSpec.powA       = permute(squeeze(nanmean(es.powA,3)),[3 1 2]);
    procSpec.powB       = permute(squeeze(nanmean(es.powB,3)),[3 1 2]);
    procSpec.powAB      = permute(squeeze(nanmean(es.powAB,3)),[3 1 2]);
    procSpec.f          = es.freq;
    procSpec.t          = es.t;
    
end

clear global pepNEV

