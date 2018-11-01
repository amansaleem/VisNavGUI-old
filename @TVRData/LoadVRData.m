function obj = LoadVRData(obj, shank_list, suffix, speed_th, nthetaphsbins, SmthTimeWindow, samplerate)
SetDirs;
if nargin < 2
    shank_list = [];
    suffix = [];
end
if nargin < 4
    speed_th = 5;
end
if nargin < 5
    nthetaphsbins = 0;
end
if nargin < 6
    SmthTimeWindow = 150;
end
if nargin < 7
    samplerate = 60;
end
obj.SpeedThresh = speed_th;
obj.SmthTimeWindow = SmthTimeWindow;
if isdir([DIRS.data2p filesep obj.animal filesep num2str(obj.iseries)])
    obj.data.es = VRLoadMultipleExpts(obj.animal, obj.iseries, obj.exptList, '2PDATA', shank_list, [], SmthTimeWindow, samplerate);
    obj.data.ephys = false;
    obj.data.twophoton = true;
elseif isdir([DIRS.multichanspikes filesep obj.animal filesep num2str(obj.iseries)])
    obj.data.es = VRLoadMultipleExpts(obj.animal, obj.iseries, obj.exptList, 'SPIKES', shank_list, suffix, SmthTimeWindow, samplerate);
    obj.data.ephys = true;
    obj.data.twophoton = false;
else
    obj.data.es = VRLoadMultipleExpts(obj.animal, obj.iseries, obj.exptList, 'BEHAV_ONLY', shank_list, suffix, SmthTimeWindow, samplerate);
    obj.data.ephys = false;
    obj.data.twophoton = false;
end

%if nthetaphsbins > 0, theta phase info must have been saved before by 
%preprocessing LFP with preprocessLFP.m
chnum = 34;%ch from which to measure theta phase. 
obj.data.es = VRbinbythetaphase(obj.data.es,nthetaphsbins,chnum);

% obj.data.es = getESDataSubset(obj.data.es, 'smthBallSpd', speed_th, []);
    
obj.CalculateSubsets();
end