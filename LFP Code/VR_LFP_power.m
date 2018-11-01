function es = VR_LFP_power(animal, iseries, iexp, beh, run_all, ChnA, ChnB)
% Calculate the power spectrum at the different frequencies and sample at
% the  es rate

global pepNEV

addpath(genpath('..\Behaviour analysis\'))
addpath(genpath('..\General\'))
addpath(genpath('..\Spike Analysis\'))

if nargin<4
    ChnA = 20;
    ChnB = 40;
    run_all = 0;
    beh  = 0;
elseif nargin<5
    run_all = 0;
    ChnB = 40;
    ChnA = 20;
elseif nargin<6
    ChnB = 40;
    ChnA = 20;
elseif nargin<7
    ChnB = 40;
end

exptInfo.animal  = animal;
exptInfo.iseries = iseries;
exptInfo.iexp    = iexp;

%% Load the behaviour data
if beh
    [~, ~ , es] = VRWheelLoad(animal, iseries, iexp);
    runSpd = es.ballspeed(~isnan(es.ballspeed));
    time   = es.sampleTimes(~isnan(es.ballspeed));
    runSmthWin = 500;
    smthRunSpd = smthInTime(runSpd,60,runSmthWin);
end
%% Load the data file
expName = [exptInfo.animal '_' num2str(exptInfo.iseries) '_' num2str(exptInfo.iexp) '.ns5'];
ns5file = ['\\ZSERVER\Data\Cerebus\' exptInfo.animal filesep num2str(exptInfo.iseries) filesep num2str(exptInfo.iexp) filesep expName]

[~,exptInfo.SamplingRateInKHZ,exptInfo.nchan] = nsopen2(ns5file);

[chan_list, ~, num_shanks] = MichiganGetLayout(animal, iseries);

% chan_list = [];
% for ichan = 1:exptInfo.nchan
%     chan_list = [chan_list {num2str(ichan)}];
% end
ChnC = chan_list(end);

exptInfo.downsample = 1000;

exptInfo.ChnA = ChnA;
exptInfo.ChnB = ChnB;
exptInfo.ChnC = ChnC;

%% Load the particular channels
timestamps.fulllength = length(pepNEV.ns.Data.data(exptInfo.ChnA,:));
timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;

ichn = 1;
num_tet = num_shanks*2;

for itet = 1:num_tet
    exptInfo.tet(itet).chans = chan_list(ichn:(ichn+3));
    ichn = ichn + 4;
    ic = 1;
    t1= tic;
    exptInfo.tet(itet).signal = decimate(double(pepNEV.ns.Data.data(exptInfo.tet(itet).chans(ic),:)), timestamps.samplingRateOrig/exptInfo.downsample);
    for ic = 2:4
        exptInfo.tet(itet).signal = exptInfo.tet(itet).signal + ...
            decimate(double(pepNEV.ns.Data.data(exptInfo.tet(itet).chans(ic),:)), timestamps.samplingRateOrig/exptInfo.downsample);
    end
    exptInfo.tet(itet).signal = exptInfo.tet(itet).signal./4;
    timeTaken = toc(t1);
    display(['Tetrode ' num2str(itet) ' took ' num2str(timeTaken)]);
end

for idipole = 1:(num_tet/2)
    exptInfo.dipole(idipole).tets = [2*idipole-1 2*idipole];
    exptInfo.dipole(idipole).signal = exptInfo.tet(exptInfo.dipole(idipole).tets(1)).signal - ...
        exptInfo.tet(exptInfo.dipole(idipole).tets(2)).signal;
end

timestamps.sampInt = 1./exptInfo.downsample;
timestamps.orig = 0:(1./timestamps.samplingRateOrig):(timestamps.fulllength./timestamps.samplingRateOrig);
timestamps.final = 0:timestamps.sampInt:...
    (timestamps.sampInt*(timestamps.fulllength*(exptInfo.downsample/timestamps.samplingRateOrig)));

ChnA = decimate(double(pepNEV.ns.Data.data(exptInfo.ChnA,:)), timestamps.samplingRateOrig/exptInfo.downsample);
ChnB = decimate(double(pepNEV.ns.Data.data(exptInfo.ChnB,:)), timestamps.samplingRateOrig/exptInfo.downsample);
%% Get the spectrograms

params.Fs=exptInfo.downsample;
exptInfo.range_low = 0.5;
exptInfo.range_high = 250;
params.fpass=[exptInfo.range_low exptInfo.range_high];
params.tapers=[5 9];
movingwin=[3 1];

procSpec.ChnA = ChnA';
procSpec.ChnB = ChnB';

[procSpec.C,procSpec.phi,procSpec.SAB,procSpec.SA,procSpec.SB,procSpec.t,procSpec.f]=cohgramc(procSpec.ChnA,procSpec.ChnB,movingwin,params);

for ipair = 1:4
    [procSpec.pairs(ipair).C,procSpec.pairs(ipair).phi,procSpec.pairs(ipair).SAB,procSpec.pairs(ipair).SA,procSpec.pairs(ipair).SB,~,~] = ...
        cohgramc( exptInfo.dipole(ipair).signal',...
        exptInfo.dipole(ipair+4).signal', movingwin,params);
end
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
    %     es.coherence= interp1q(procSpec.t', procSpec.C, es.sampleTimes);
    %     es.cohPhi   = interp1q(procSpec.t', procSpec.phi, es.sampleTimes); %procSpec.phi;
    %     es.powA     = interp1q(procSpec.t', procSpec.SA, es.sampleTimes); %procSpec.SA;
    %     es.powB     = interp1q(procSpec.t', procSpec.SB, es.sampleTimes); %procSpec.SB;
    %     es.powAB    = interp1q(procSpec.t', procSpec.SAB, es.sampleTimes); %procSpec.SAB;
    
    es.coherence= procSpec.C;
    es.cohPhi   = procSpec.phi;
    es.powA     = procSpec.SA;
    es.powB     = procSpec.SB;
    es.powAB    = procSpec.SAB;
    
    es.pairs    = procSpec.pairs;
    % down-sampling the VR data
    es.traj          = interp1q(es.sampleTimes, es.traj,      procSpec.t');
    es.trajspeed     = interp1q(es.sampleTimes, es.trajspeed,      procSpec.t');
    es.ballspeed     = interp1q(es.sampleTimes, es.ballspeed,      procSpec.t');
    es.distTrav      = interp1q(es.sampleTimes, es.distTrav,      procSpec.t');
    es.totDistTrav   = interp1q(es.sampleTimes, es.totDistTrav,      procSpec.t');
    es.trajPercent   = interp1q(es.sampleTimes, es.trajPercent,      procSpec.t');
    es.smthBallSpd   = interp1q(es.sampleTimes, es.smthBallSpd,      procSpec.t');
    es.contrast      = interp1q(es.sampleTimes, es.contrast,      procSpec.t');
    es.start         = interp1q(es.sampleTimes, es.start,      procSpec.t');
    es.blanks        = interp1q(es.sampleTimes, es.blanks,      procSpec.t');
    es.active        = interp1q(es.sampleTimes, es.active,      procSpec.t');
    es.rewardPos     = interp1q(es.sampleTimes, es.rewardPos,      procSpec.t');
    es.outcome       = interp1q(es.sampleTimes, es.outcome,      procSpec.t');
    es.roomLength    = interp1q(es.sampleTimes, es.roomLength,      procSpec.t');
    es.lick          = interp1q(es.sampleTimes, es.lick,      procSpec.t');
    es.reward        = interp1q(es.sampleTimes, es.reward,      procSpec.t');
    es.trialID       = round(interp1q(es.sampleTimes, es.trialID,      procSpec.t'));
    es.gain          = interp1q(es.sampleTimes, es.gain,      procSpec.t');
    es.iexp          = round(interp1q(es.sampleTimes, es.iexp,      procSpec.t'));
    
    es.sampleRate    = 1./movingwin(2);
    es.origSampleTimes = es.sampleTimes;
    es.sampleTimes   = procSpec.t';
else
    es.coherence= procSpec.C;
    es.cohPhi   = procSpec.phi;
    es.powA     = procSpec.SA;
    es.powB     = procSpec.SB;
    es.powAB    = procSpec.SAB;
end

clear global pepNEV

