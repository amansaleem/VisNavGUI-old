function [es, procSpec] = LFP_power_only...
    (data, chan2consider, samplingRate, movingwin)
% Calculate the power spectrum at the different frequencies and sample at
% the  es rate

if nargin<3
    samplingRate = 1250;
end
if nargin<4
    movingwin=[3 1];
end

data_2consider = data(chan2consider,:);

% %% Load the particular channels
%% Get the spectrogram

params.Fs = samplingRate;
exptInfo.range_low = 0.5;
exptInfo.range_high = 250;
params.fpass=[exptInfo.range_low exptInfo.range_high];
params.tapers=[5 9];

[procSpec.SA,procSpec.t,procSpec.f]=...
            mtspecgramc(data_2consider',movingwin,params);
params.movingwin = movingwin;

%% Downsample the es rate
es.freq     = procSpec.f;
es.t        = procSpec.t;
es.params   = params;
es.powA     = procSpec.SA;
es.sampleRate       = 1./movingwin(2);