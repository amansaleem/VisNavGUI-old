function [exptInfo, procSpec, timestamps] = LFP_loadFile(animal, iseries, iexp, ChnA, ChnB)

global pepNEV

if nargin<4
    ChnA = 40; % V1 electrodes should be 33:64
    ChnB = 55;
end

exptInfo.animal  = animal;
exptInfo.iseries = iseries;
exptInfo.iexp    = iexp;

expName = [exptInfo.animal '_' num2str(exptInfo.iseries) '_' num2str(exptInfo.iexp) '.ns5'];
ns5file = ['\\ZSERVER\Data\Cerebus\' exptInfo.animal filesep num2str(exptInfo.iseries) filesep num2str(exptInfo.iexp) filesep expName];

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

%% Get the spectrogram

params.Fs=exptInfo.downsample;
exptInfo.range_low = 0.5;
exptInfo.range_high = 100;
params.fpass=[exptInfo.range_low exptInfo.range_high];
params.tapers=[5 9];
movingwin=[3 1];%[0.5 0.2];% [window_size window_shift] 
procSpec.ChnA = ChnA';
procSpec.ChnB = ChnB';
procSpec.ChnC = ChnC_down';

% This is calling chronux
[procSpec.C,procSpec.phi,procSpec.SAB,procSpec.SA,procSpec.SB,procSpec.t,procSpec.f]=cohgramc(procSpec.ChnA,procSpec.ChnB,movingwin,params);
