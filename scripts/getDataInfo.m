function [dataInfo]= getDataInfo(animal)
% AA 2010-07-22 created

SetDirs;

% get folder names with date format
folders = sprintf('%s%s%s',[DIRS.ball filesep],[animal filesep],'*-*-*');


animalDir = sprintf('%s%s%s',[DIRS.ball filesep],animal);
% create a data structure that contains data file names for each day
% number of sessions and number of trials in each session
f = dir(animalDir);
f(1) = [];
f(1) = [];
days = struct;
if ismac
    temp_start = 2;
else
    temp_start = 1;
end
idxTemp = 0;
for i = temp_start:numel(f)
    s= sprintf('%s%s',[animalDir filesep],f(i).name);
    
    files = what(s);
    if strcmp(f(i).name, 'EXP.mat') | strfind(f(i).name, '_log')>1;
        continue
    end
    idxTemp = idxTemp + 1;
    days(idxTemp).date = f(i).name;
    days(idxTemp).datafiles = files.mat;
    days(idxTemp).animalDir = animalDir;
    
    starti = regexpi(days(idxTemp).datafiles, 'session_[0-9]+_trial');
    endi = regexpi(days(idxTemp).datafiles, '_trial');
    sessionIds = [];
    trialIds = [];
    trials=[];
    for j=1:numel(days(idxTemp).datafiles)
        sessionIds(j) = str2num(days(idxTemp).datafiles{j}(starti{j}+8:endi{j}-1));
        trialIds(j) = sscanf((days(idxTemp).datafiles{j}((endi{j}+7):(length(days(idxTemp).datafiles{j})-4))),'%d');
    end
    
    days(i).sessions = unique(sessionIds);
    for k=1: numel(days(i).sessions)
        trials(k) = max(trialIds(find(sessionIds==(days(i).sessions(k)))));
    end
    days(i).trials = trials;
    
end

dataInfo = days;