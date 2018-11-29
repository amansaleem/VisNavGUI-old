% function DIRS = SetDefaultDirs
global DIRS serverName server2Name server3Name

[~,hostName] = system('hostname');
hostName = hostName(1:end-1);

switch hostName
    case 'saleem12'
        serverName = 'Y:\Saleem Lab\';
    case 'saleem05'
        serverName = 'Y:\Saleem Lab\';
    case 'solomon02'
        serverName = 'Y:\';
    case 'saleem01'
        serverName = 'S:\';
    case 'saleem11'
        %serverName = 'Z:\';
        serverName = 'Z:\Data\NewFileServer\Data\Subjects'; % temporary testing of new server format
end

if strcmp(hostName,'saleem11')
    
    animal= 'M180910_MMc_sms';
    iseries='1030';
    iexp='';
    suffix='V1';
    
    DIRS.animals        = serverName;
    DIRS.ball           = fullfile(serverName,animal,'VRBehaviour'); % This is where the VR data is
    DIRS.multichanspikes= fullfile(serverName,animal,'ePhys',iseries); % This is where the kwikFiles, photodiode (VR processed files are)
    DIRS.data           = fullfile(serverName,animal,'BonsaiLog'); % This is where all the processed files are. Also where the p-files are stored

else
    DIRS.ball           = fullfile(serverName,'Data','Behav'); % This is where the VR data is
    DIRS.multichanspikes= fullfile(serverName,'Data','ePhys'); % This is where the kwikFiles, photodiode (VR processed files are)
    DIRS.data           = fullfile(serverName,'Data','processed'); % This is where all the processed files are. Also where the p-files are stored

    DIRS.xfiles         = fullfile(serverName,'Code','SLVS','xfiles'); % stimulus running files.

    DIRS.stimInfo       = fullfile(serverName,'Data','stimInfo'); % Not really necessary
    DIRS.spikes         = fullfile(serverName,'Data','Spikes'); % seems redundant
    DIRS.expInfo        = fullfile(serverName,'Data','expInfo');

end

% DIRS.camera         = fullfile(serverName,'Data','Intrinsic');
% DIRS.Intrinsic      = fullfile(serverName,'Data','Intrinsic'); % Piperico changed to zserver
% DIRS.EyeCamera      = fullfile(serverName,'Data','EyeCamera'); % Piperico changed to zserver 
% DIRS.EyeTrack       = fullfile(serverName,'Data','EyeTrack');
% DIRS.michigan       = fullfile(serverName,'Data','michigan');
% DIRS.Cerebus        = fullfile(serverName,'Data','Cerebus');
% DIRS.behavior       = fullfile(serverName,'Data','behavior');
% DIRS.mouselogs      = fullfile(serverName,'Data','logs','mouse','behavior');

% DIRS.Stacks         = fullfile(server3Name,'Data','Stacks'); % Piperico changed to zserver3
