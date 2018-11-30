%% Tomaso Muzzu - UCL - 10/04/2018

%% photodiode signal: save original signal downsampled at 1kHz, timestamps of stimulus onsets and offsets
function [stimON, stimOFF, ACInfo] = getPhotoDiodeSignal(animal, iseries, iexp)
    
    global DIRS
    DIRname  = [DIRS.data filesep animal filesep num2str(iseries) filesep num2str(iexp)];
    % load the mat file containing stim ON and OFF events recorded with OE
    
    if (nargout<3) & exist([DIRname filesep 'AC_Info.mat'],'file')
        load([DIRname filesep 'AC_Info.mat'],'stimON','stimOFF');
        % else load the analog channels now, find when the timestamps of the VS
        % onsets and offsets
    else
        display('Analogue channels have not been converted yet or some error reading the saved file. Please select folder containing relevant OE recording for recovering ACs info.');
        ACInfo = getEPhysAnalogSignals;
        % data(:,1) = sync pulse signal
        % data(:,2) = photodiode signal
        % data(:,3) = signal A from rotary encoder
        % data(:,4) = signal B from rotary encoder
        
        % Plot the PD signal --added 21/08/18 MM 
        figure; plot(ACInfo.Data(:,2))
        
        % resample and normalise signal
        [temp_pks,temp_locs] = findpeaks(ACInfo.Data(:,2)/max(ACInfo.Data(:,2)));
        pks = temp_pks(temp_pks>0.5);
        locs = temp_locs(temp_pks>0.5);
        clear temp_pks temp_locs
        % find first peaks signalling the onset of the white square
        PksTimeDifference = diff(locs);
        temp_VSOnsets = find(PksTimeDifference>900);
        VSOnsetsIndecesON = [locs(1); locs(temp_VSOnsets+1)]; % save indexes of stimulus onsets
        % find last peaks signalling the offset of the white square
        VSOnsetsIndecesOFF = [locs(temp_VSOnsets-1); locs(end)]; % save indexes of stimulus offsets
        if ~isempty(find(diff(VSOnsetsIndecesON)>3*ACInfo.SamplingRateOE)) || ~isempty(find(diff(VSOnsetsIndecesOFF)>3*ACInfo.SamplingRateOE))
            fprintf('Error during photodiode signal extraction. Please check it manually \n');
            return
        end
        ACInfo.stimON = ACInfo.Timestamps(VSOnsetsIndecesON)-min(ACInfo.Timestamps);
        ACInfo.stimOFF = ACInfo.Timestamps(VSOnsetsIndecesOFF)-min(ACInfo.Timestamps);
        stimON = ACInfo.stimON;
        stimOFF = ACInfo.stimOFF;
        % save Onsets and Offsets into relevant folder
        save([DIRname filesep 'AC_Info.mat'],'-struct', 'ACInfo');
        
    end
    
    % verify detection of stimulus presentation
%     figure
%     plot(ACInfo.Timestamps-min(ACInfo.Timestamps),ACInfo.Data(:,2))
%     hold on
%     plot(ACInfo.Timestamps(VSOnsetsIndecesON)-min(ACInfo.Timestamps),ACInfo.Data(VSOnsetsIndecesON,2),'s')   
%     hold on
%     plot(ACInfo.Timestamps(VSOnsetsIndecesOFF)-min(ACInfo.Timestamps),ACInfo.Data(VSOnsetsIndecesOFF,2),'s')   
%     
end    