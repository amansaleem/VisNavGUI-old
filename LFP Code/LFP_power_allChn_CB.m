function [es procSpec] = LFP_power_allChn_CB(animal, iseries, iexp, chan2consider, beh)
% Calculate the power spectrum at the different frequencies and sample at
% the  es rate

global pepNEV

SetDefaultDirs

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

beh = 0;

exptInfo.animal  = animal;
exptInfo.iseries = iseries;
exptInfo.iexp    = iexp;

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

exptInfo.downsample = 1000;

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
    exptInfo.Chn(ichn).data = decimate(double(pepNEV.ns.Data.data(exptInfo.Chn(ichn).ChnLoc,:)), timestamps.samplingRateOrig/exptInfo.downsample);
    %     exptInfo.Chn(ichn).data = MKnotch(exptInfo.Chn(ichn).data, 4, exptInfo.downsample, 50);
end

%% Get the spectrogram

params.Fs=exptInfo.downsample;
exptInfo.range_low = 0.5;
exptInfo.range_high = 250;
params.fpass=[exptInfo.range_low exptInfo.range_high];
params.tapers=[5 9];
movingwin=[3 1];

for chn1 = 1%:numChn
    for chn2 = 1:numChn
        [procSpec.C(chn1,chn2,:,:),procSpec.phi(chn1,chn2,:,:),procSpec.SAB(chn1,chn2,:,:),procSpec.SA(chn1,chn2,:,:),procSpec.SB(chn1,chn2,:,:),procSpec.t(chn1,chn2,:),procSpec.f(chn1,chn2,:)]=...
            cohgramc(exptInfo.Chn(chn1).data',exptInfo.Chn(chn2).data',movingwin,params);
    end
end
for chn1 = 1%:numChn
    for chn2 = 1:numChn
        procSpec.C(chn2,chn1,:,:) = procSpec.C(chn1,chn2,:,:);
        procSpec.phi(chn2,chn1,:,:) = procSpec.phi(chn1,chn2,:,:);
        procSpec.SAB(chn2,chn1,:,:) = procSpec.SAB(chn1,chn2,:,:);
        procSpec.SA(chn2,chn1,:,:) = procSpec.SA(chn1,chn2,:,:);
        procSpec.SB(chn2,chn1,:,:) = procSpec.SB(chn1,chn2,:,:);
        procSpec.f(chn2,chn1,:,:) = procSpec.f(chn1,chn2,:);
        procSpec.t(chn2,chn1,:,:) = procSpec.t(chn1,chn2,:);
    end
end

%% Downsample the es rate

es.freq         = procSpec.f;
es.origLFPtime  = procSpec.t;

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

clear global pepNEV

