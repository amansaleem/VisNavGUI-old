function createScreenTimeFiles(animal, iseries, expt_list)

SetDefaultDirs
global pepNEV
nevDir = fullfile(DIRS.Cerebus,animal,num2str(iseries));

KlustersDir = fullfile(DIRS.multichanspikes,animal,num2str(iseries));

% nevSamplingRateInKHZ = 30;
period = 1;%nevSamplingRateInKHZ*1000 / samplingRate;

for ns5ID = 1:length(expt_list)
    nevFileName = [nevDir filesep num2str(expt_list(ns5ID)) filesep animal '_' num2str(iseries) '_' num2str(expt_list(ns5ID)) '.nev'];
    openNEV(nevFileName);
    nevopen(nevFileName);
    
    screenTimes  = pepNEV.sync.timestamps;
    
    
    screenTimes(find(diff(screenTimes)<100) + 1) = [];
    %                 screenTimes(1) = [];
    
    %             % open .ns* file and load data
    %             nsFileName = [nevDir filesep num2str(expt_list(ns5ID)) filesep animal '_' num2str(iseries) '_' num2str(expt_list(ns5ID)) '.ns5'];
    %             [~,SamplingRateInKHZ,nchan] = nsopen(nsFileName);
    screenFileName = [KlustersDir filesep animal '_' num2str(iseries) '_' num2str(expt_list(ns5ID)) '_screenTimes.mat'];
    
    nevclose;
    
    screenTimes(1) = [];
    save(screenFileName, 'screenTimes');
    display(['Done with experiment ' num2str(expt_list(ns5ID))]);
end
