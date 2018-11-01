function createScreenTimeFiles2p(animal, iseries, expt_list)

SetDefaultDirs
TimelineDir = fullfile(DIRS.expInfo,animal,num2str(iseries));

TwophotonDir = fullfile(DIRS.data2p,animal,num2str(iseries));

photodiodeID = 2;
neuralFrames = 4;
TimelineSamplingrate = 1000;

for fileID = 1:length(expt_list)
    TimelineFileName = [TimelineDir filesep num2str(expt_list(fileID)) filesep num2str(iseries) '_' num2str(expt_list(fileID)) '_' animal '_Timeline.mat'];
    S = load(TimelineFileName);
    
    photodiode = medfilt1(S.Timeline.rawDAQData(:,photodiodeID),5);
    maxVoltage = max(photodiode);
    screenTimes  = find(diff(sign(photodiode - maxVoltage/2))~=0);
    screenTimes = screenTimes(:)';
    
    screenTimes(find(diff(screenTimes)<10*TimelineSamplingrate/1000) + 1) = [];
    screenFileName = [TwophotonDir filesep animal '_' num2str(iseries) '_' num2str(expt_list(fileID)) '_screenTimes.mat'];
    
%     screenTimes(1) = [];

    neuralFrameTimes = find(diff(S.Timeline.rawDAQData(:,neuralFrames))>0);
    neuralFrameTimes = neuralFrameTimes(:)';
    
    save(screenFileName, 'screenTimes', 'neuralFrameTimes');
    display(['Done with experiment ' num2str(expt_list(fileID))]);
end