function [screenTimes, stimOn, stimOff] = loadStimTimes(animal, iseries, iexp)

SetDefaultDirs
global pepNEV
nevDir = fullfile(DIRS.Cerebus,animal,num2str(iseries));

KlustersDir = fullfile(DIRS.multichanspikes,animal,num2str(iseries));

% nevSamplingRateInKHZ = 30;
period = 1;%nevSamplingRateInKHZ*1000 / samplingRate;

nevFileName = [nevDir filesep num2str(iexp) filesep animal '_' num2str(iseries) '_' num2str(iexp) '.nev'];
nevopen(nevFileName);

screenTimes  = pepNEV.sync.timestamps;


screenTimes(find(diff(screenTimes)<100) + 1) = [];

stimOn = screenTimes(1:2:end);
stimOff = screenTimes(2:2:end);

nevclose;

