function [es procSpec] = VR_LFP_power_safe(animal, iseries, iexp, beh, ChnA, ChnB, sampleUP)
% Calculate the power spectrum at the different frequencies and sample at
% the  es rate

global pepNEV

addpath(genpath('..\Behaviour analysis\'))
addpath(genpath('..\General\'))
addpath(genpath('..\Spike Analysis\'))

if nargin<4
    ChnA = 20;
    ChnB = 40;
    beh  = 0;
elseif nargin<5
    ChnB = 40;
    ChnA = 20;
elseif nargin<6
    ChnB = 40;
end
if nargin<7
    sampleUP = 0;
    display(['Sample up is set to: ' num2str(sampleUP) '  continue in 1s']);
end

exptInfo.animal  = animal;
exptInfo.iseries = iseries;
exptInfo.iexp    = iexp;

%% Load the behaviour data
if beh
    try
        [~, ~ , es] = VRWheelLoad(animal, iseries, iexp);
        try
%             es = VRLoadMultipleExpts(animal, iseries, iexp,'SPIKES');
        catch
            display('No loadable spikes')
        end
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

% chan_list = [];
% for ichan = 1:exptInfo.nchan
%     chan_list = [chan_list {num2str(ichan)}];
% end
ChnC = chan_list(end);

exptInfo.downsample = 1000;

exptInfo.ChnA = chan_list(ChnA);
exptInfo.ChnB = chan_list(ChnB);
exptInfo.ChnC = chan_list(ChnC);


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

% if ~beh
timestamps.fulllength = length(pepNEV.ns.Data.data(exptInfo.ChnC,:));
timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;

exptInfo.runSmthWin = 300;

ChnC = zeros(1, size(pepNEV.ns.Data.data,2));
maxC = max(pepNEV.ns.Data.data(exptInfo.ChnC,:));
minC = min(pepNEV.ns.Data.data(exptInfo.ChnC,:));
ChnC((abs(diff(pepNEV.ns.Data.data(exptInfo.ChnC,:)-minC)./(maxC-minC)))>0.5) = 1;
ChnC_down = sum(reshape(ChnC(1:end-rem(size(pepNEV.ns.Data.data,2),250)),250,[]),1);

timestamps.ChnCsampInt = 1./exptInfo.downsample;
timestamps.ChnCorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
timestamps.ChnC = 0:timestamps.ChnBsampInt:...
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
% end
%% Get the spectrogram

params.Fs=exptInfo.downsample;
exptInfo.range_low = 0.5;
exptInfo.range_high = 100;
params.fpass=[exptInfo.range_low exptInfo.range_high];
params.tapers=[5 9];
movingwin=[3 1];%[1 0.2]; [0.5 0.2];% [window_size window_shift] 
display(['Note: window size= ' num2str(movingwin(1)) 's  Window slide: ' num2str(movingwin(2))]);

procSpec.ChnA = ChnA';
procSpec.ChnB = ChnB';
procSpec.ChnC = ChnC_down';

[procSpec.C,procSpec.phi,procSpec.SAB,procSpec.SA,procSpec.SB,procSpec.t,procSpec.f]=cohgramc(procSpec.ChnA,procSpec.ChnB,movingwin,params);

%% Hilbert info
procSpec.delta.frange  = [1 4.5];
procSpec.delta.hilbertA = LFPsignals(procSpec.ChnA', 1000, procSpec.delta.frange(1), procSpec.delta.frange(2));
procSpec.delta.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.delta.frange(1), procSpec.delta.frange(2));
procSpec.delta.hilbertC = LFPsignals(procSpec.ChnC', 1000, procSpec.delta.frange(1), procSpec.delta.frange(2));

procSpec.theta.frange  = [4.5 9];
procSpec.theta.hilbertA = LFPsignals(procSpec.ChnA', 1000, procSpec.theta.frange(1), procSpec.theta.frange(2));
procSpec.theta.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.theta.frange(1), procSpec.theta.frange(2));
procSpec.theta.hilbertC = LFPsignals(procSpec.ChnC', 1000, procSpec.theta.frange(1), procSpec.theta.frange(2));

procSpec.beta.frange  = [10 30];
procSpec.beta.hilbertA = LFPsignals(procSpec.ChnA', 1000, procSpec.beta.frange(1), procSpec.beta.frange(2));
procSpec.beta.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.beta.frange(1), procSpec.beta.frange(2));

procSpec.gamma.frange  = [52 95];
procSpec.gamma.hilbertA = LFPsignals(procSpec.ChnA', 1000, procSpec.gamma.frange(1), procSpec.gamma.frange(2));
procSpec.gamma.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.gamma.frange(1), procSpec.gamma.frange(2));

procSpec.gamma_narrow.frange  = [58 76];
procSpec.gamma_narrow.hilbertB = LFPsignals(procSpec.ChnB', 1000, procSpec.gamma_narrow.frange(1), procSpec.gamma_narrow.frange(2));

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
if beh
    
    if sampleUP
        
        es.coherence= interp1q(procSpec.t', procSpec.C, es.sampleTimes);
        es.cohPhi   = interp1q(procSpec.t', procSpec.phi, es.sampleTimes); %procSpec.phi;
        es.powA     = interp1q(procSpec.t', procSpec.SA, es.sampleTimes); %procSpec.SA;
        es.powB     = interp1q(procSpec.t', procSpec.SB, es.sampleTimes); %procSpec.SB;
        es.powAB    = interp1q(procSpec.t', procSpec.SAB, es.sampleTimes); %procSpec.SAB;
        
        es.LFP_A    = procSpec.ChnA';
        es.LFP_B    = procSpec.ChnB';
        
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
%         es.LFP_A    = procSpec.ChnA';
%         es.LFP_B    = procSpec.ChnB';
        
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
        es.bands.delta = procSpec.delta.frange;
        es.bands.theta = procSpec.theta.frange;
        es.bands.beta = procSpec.beta.frange;
        es.bands.gamma = procSpec.gamma.frange;
        es.bands.gamma_narrow = procSpec.gamma_narrow.frange;
        
        % delta
        es.A.delta.hill    = interp1q(timestamps.ChnA', procSpec.delta.hilbertA.hill', procSpec.t');
        es.B.delta.hill    = interp1q(timestamps.ChnB', procSpec.delta.hilbertB.hill', procSpec.t');
        es.C.delta.hill    = interp1q(timestamps.ChnC_down', procSpec.delta.hilbertC.hill', procSpec.t');
        % theta
        es.A.theta.hill    = interp1q(timestamps.ChnA', procSpec.theta.hilbertA.hill', procSpec.t');
        es.B.theta.hill    = interp1q(timestamps.ChnB', procSpec.theta.hilbertB.hill', procSpec.t');
        es.C.theta.hill    = interp1q(timestamps.ChnC_down', procSpec.theta.hilbertC.hill', procSpec.t');
        % beta
        es.A.beta.hill    = interp1q(timestamps.ChnA', procSpec.beta.hilbertA.hill', procSpec.t');
        es.B.beta.hill    = interp1q(timestamps.ChnB', procSpec.beta.hilbertB.hill', procSpec.t');
        % gamma
        es.A.gamma.hill    = interp1q(timestamps.ChnA', procSpec.gamma.hilbertA.hill', procSpec.t');
        es.B.gamma.hill    = interp1q(timestamps.ChnB', procSpec.gamma.hilbertB.hill', procSpec.t');
        % gamma_narrow
        es.B.gamma_narrow.hill    = interp1q(timestamps.ChnB', procSpec.gamma_narrow.hilbertB.hill', procSpec.t');
        
    end
else
    es.coherence= procSpec.C;
    es.cohPhi   = procSpec.phi;
    es.powA     = procSpec.SA;
    es.powB     = procSpec.SB;
    es.powAB    = procSpec.SAB;
    
    es.smthRunSpd = interp1q(exptInfo.time',exptInfo.smthRunSpd', procSpec.t');
    
end

clear global pepNEV

