function [es spikeTimes lfp] = stLFP_oldData(animal, iseries, iexp, ChnA, beh)
% Calculate the power spectrum at the different frequencies and sample at
% the  es rate

global pepNEV

SetDirs

addpath(genpath('..\Behaviour analysis\'))
addpath(genpath('..\General\'))
addpath(genpath('..\Spike Analysis\'))

if nargin<4
    ChnA = 20;
    beh  = 1;
elseif nargin<5
    beh = 1;
end

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
        spikeTimes = UnitLoad([DIRS.spikes filesep 'Klustered'], animal, iseries, iexp);
        es = getVRspikeTimes(animal, iseries, iexp);
        spikeTimes = UnitLoad([DIRS.spikes filesep 'Klustered'], animal, iseries, iexp);
    end
    runSpd = es.ballspeed;%(~isnan(es.ballspeed));
    time   = es.sampleTimes;%(~isnan(es.ballspeed));
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
% ChnC = chan_list(end);

exptInfo.downsample = 1000;

exptInfo.ChnA = ChnA;
% exptInfo.ChnB = ChnB;
% exptInfo.ChnC = ChnC;


%% Load the particular channels
timestamps.fulllength = length(pepNEV.ns.Data.data(exptInfo.ChnA,:));
timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;

lfp.data = decimate(double(pepNEV.ns.Data.data(exptInfo.ChnA,:)), timestamps.samplingRateOrig/exptInfo.downsample);
timestamps.ChnAsampInt = 1./exptInfo.downsample;
timestamps.ChnAorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
lfp.timestamps = 0:timestamps.ChnAsampInt:...
    (timestamps.ChnAsampInt*(timestamps.fulllength*(exptInfo.downsample/timestamps.samplingRateOrig)));

% ChnB = decimate(double(pepNEV.ns.Data.data(exptInfo.ChnB,:)), timestamps.samplingRateOrig/exptInfo.downsample);
% timestamps.ChnBsampInt = 1./exptInfo.downsample;
% timestamps.ChnBorig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
% timestamps.ChnB = 0:timestamps.ChnBsampInt:...
%     (timestamps.ChnBsampInt*(timestamps.fulllength*(exptInfo.downsample/timestamps.samplingRateOrig)));

% %% Get the spectrogram
% 
% params.Fs=exptInfo.downsample;
% exptInfo.range_low = 0.5;
% exptInfo.range_high = 250;
% params.fpass=[exptInfo.range_low exptInfo.range_high];
% params.tapers=[5 9];
% movingwin=[3 1];
% procSpec.ChnA = ChnA';
% procSpec.ChnB = ChnB';
% 
% % [procSpec.C,procSpec.phi,procSpec.SAB,procSpec.SA,procSpec.SB,procSpec.t,procSpec.f]=cohgramc(procSpec.ChnA,procSpec.ChnB,movingwin,params);
% 
% %% Hilbert info
% procSpec.hilbertA = LFPsignals(procSpec.ChnA', 1000, 4.5, 9);
% procSpec.hilbertB = LFPsignals(procSpec.ChnB', 1000, 4.5, 9);

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
%     procSpec.t = procSpec.t - screenTimes(1)./(1000*nevSamplingRateInKHZ);
    lfp.timestamps = lfp.timestamps - screenTimes(1)./(1000*nevSamplingRateInKHZ);
end

% %% Downsample the es rate
% 
% es.freq         = procSpec.f;
% es.origLFPtime  = procSpec.t;
% if beh
%     %     es.coherence= interp1q(procSpec.t', procSpec.C, es.sampleTimes);
%     %     es.cohPhi   = interp1q(procSpec.t', procSpec.phi, es.sampleTimes); %procSpec.phi;
%     %     es.powA     = interp1q(procSpec.t', procSpec.SA, es.sampleTimes); %procSpec.SA;
%     %     es.powB     = interp1q(procSpec.t', procSpec.SB, es.sampleTimes); %procSpec.SB;
%     %     es.powAB    = interp1q(procSpec.t', procSpec.SAB, es.sampleTimes); %procSpec.SAB;
%     
%     es.coherence= procSpec.C;
%     es.cohPhi   = procSpec.phi;
%     es.powA     = procSpec.SA;
%     es.powB     = procSpec.SB;
%     es.powAB    = procSpec.SAB;
%     
%     % down-sampling the VR data
%     es.traj          = interp1q(es.sampleTimes, es.traj,      procSpec.t');
%     es.trajspeed     = interp1q(es.sampleTimes, es.trajspeed,      procSpec.t');
%     es.ballspeed     = interp1q(es.sampleTimes, es.ballspeed,      procSpec.t');
% %     es.distTrav      = interp1q(es.sampleTimes, es.distTrav,      procSpec.t');
% %     es.totDistTrav   = interp1q(es.sampleTimes, es.totDistTrav,      procSpec.t');
% %     es.trajPercent   = interp1q(es.sampleTimes, es.trajPercent,      procSpec.t');
% %     es.smthBallSpd   = interp1q(es.sampleTimes, es.smthBallSpd,      procSpec.t');
% %     es.smthTrajSpd   = interp1q(es.sampleTimes, es.smthTrajSpd,      procSpec.t');
% %     es.contrast      = interp1q(es.sampleTimes, es.contrast,      procSpec.t');
% %     es.start         = interp1q(es.sampleTimes, es.start,      procSpec.t');
% %     es.blanks        = interp1q(es.sampleTimes, es.blanks,      procSpec.t');
% %     es.active        = interp1q(es.sampleTimes, es.active,      procSpec.t');
% %     es.rewardPos     = interp1q(es.sampleTimes, es.rewardPos,      procSpec.t');
% %     es.outcome       = interp1q(es.sampleTimes, es.outcome,      procSpec.t');
% %     es.roomLength    = interp1q(es.sampleTimes, es.roomLength,      procSpec.t');
% %     es.lick          = interp1q(es.sampleTimes, es.lick,      procSpec.t');
% %     es.reward        = interp1q(es.sampleTimes, es.reward,      procSpec.t');
% %     es.trialID       = round(interp1q(es.sampleTimes, es.trialID,      procSpec.t'));
% %     es.gain          = interp1q(es.sampleTimes, es.gain,      procSpec.t');
% %     es.iexp          = round(interp1q(es.sampleTimes, es.iexp,      procSpec.t'));
%     
%     es.sampleRate    = 1./movingwin(2);
%     es.origSampleTimes = es.sampleTimes;
%     es.sampleTimes   = procSpec.t';
%     
%     % downsampling the Hilbert data
% %     es.hilbertA.phase   = interp1q(timestamps.ChnA', procSpec.hilbertA.phase', procSpec.t');
% %     es.hilbertA.hill    = interp1q(timestamps.ChnA', procSpec.hilbertA.hill', procSpec.t');
% %     es.hilbertA.instFreq= interp1q(timestamps.ChnA', procSpec.hilbertA.instFreq', procSpec.t');
% %     es.hilbertA.power   = interp1q(timestamps.ChnA', procSpec.hilbertA.power', procSpec.t');
% %     es.hilbertA.band    = procSpec.hilbertA.band;
% %     
% %     es.hilbertB.phase   = interp1q(timestamps.ChnB', procSpec.hilbertB.phase', procSpec.t');
% %     es.hilbertB.hill    = interp1q(timestamps.ChnB', procSpec.hilbertB.hill', procSpec.t');
% %     es.hilbertB.instFreq= interp1q(timestamps.ChnB', procSpec.hilbertB.instFreq', procSpec.t');
% %     es.hilbertB.power   = interp1q(timestamps.ChnB', procSpec.hilbertB.power', procSpec.t');
% %     es.hilbertB.band    = procSpec.hilbertB.band;
% else
%     es.coherence= procSpec.C;
%     es.cohPhi   = procSpec.phi;
%     es.powA     = procSpec.SA;
%     es.powB     = procSpec.SB;
%     es.powAB    = procSpec.SAB;
% end

clear global pepNEV

