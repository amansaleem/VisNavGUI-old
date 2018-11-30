%% Tomaso Muzzu - UCL - 10/04/2018



for i = 1:length(FoldersList)
    
    FileList = dir([FoldersList{i} filesep '101_ADC*.continuous']);
    clear data timestamps info
    
    MpepFlag = 1;
    [data, timestamps, ACInfo] = getEPhysAnalogSignals( MpepFlag, FoldersList, FileList);

    ACInfo = getSyncPulseSignal(data,timestamps,ACInfo);
    
    ACInfo = getPhotoDiodeSignal(data,timestamps,ACInfo, MpepFlag);
    
    ACInfo = getSpeedSignal(data,timestamps,ACInfo);
    
end