%% Tomaso Muzzu - UCL - 10/04/2018

%% syncpulse signal: save original signal downsampled at 1kHz, distance between two consecutive pulses falling-to-rising edge
function ACInfo = getSyncPulseSignal(animal, iseries, iexp)
    
    global DIRS
    DIRname  = [DIRS.data filesep animal filesep num2str(iseries) filesep num2str(iexp)];
    
    if (nargout<3) && exist([DIRname filesep 'AC_Info.mat'],'file')
        load([DIRname filesep 'AC_Info.mat'],'SyncPulse');
        % else load the analog channels now, find when the timestamps of the VS
        % onsets and offsets
    else
        display('Analogue channels have not been converted yet or some error reading the saved file. Please select the folder containing relevant OE recording for recovering ACs info.');
        ACInfo = getEPhysAnalogSignals;
        ACInfo.SamplingFactor = ACInfo.SamplingRateOE/1000;
        if isfield(ACInfo,'SessionStarts')
            for i = 1:length(ACInfo.SessionStarts)+1
                if i == 1
                    ACInfo.TimestampsDownsampled{i} = linspace(min(ACInfo.Timestamps(1:ACInfo.SessionStarts(i)-1))-min(ACInfo.Timestamps(1:ACInfo.SessionStarts(i)-1)),...
                                                                max(ACInfo.Timestamps(1:ACInfo.SessionStarts(i)-1))-min(ACInfo.Timestamps(1:ACInfo.SessionStarts(i)-1)),...
                                                                round(length(ACInfo.Timestamps(1:ACInfo.SessionStarts(i)-1))/ACInfo.SamplingFactor))'; %resample at 1kHz
                    ACInfo.SyncPulse{i} = interp1qr(ACInfo.Timestamps(1:ACInfo.SessionStarts(i)-1)-min(ACInfo.Timestamps(1:ACInfo.SessionStarts(i)-1)),...
                                                    ACInfo.Data(1:ACInfo.SessionStarts(i)-1,1),...
                                                    ACInfo.TimestampsDownsampled{i});
                elseif i == length(ACInfo.SessionStarts)+1
                    ACInfo.TimestampsDownsampled{i} = linspace(min(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):end))-min(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):end)),...
                                                                max(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):end))-min(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):end)),...
                                                                round(length(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):end))/ACInfo.SamplingFactor))'; %resample at 1kHz
                    ACInfo.SyncPulse{i} = interp1qr(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):end)-min(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):end)),...
                                                    ACInfo.Data(ACInfo.SessionStarts(i-1):end,1),...
                                                    ACInfo.TimestampsDownsampled{i});
                else 
                    ACInfo.TimestampsDownsampled{i} = linspace(min(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):ACInfo.SessionStarts(i)))-min(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):ACInfo.SessionStarts(i))),...
                                                                max(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):ACInfo.SessionStarts(i)))-min(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):ACInfo.SessionStarts(i))),...
                                                                round(length(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):ACInfo.SessionStarts(i)))/ACInfo.SamplingFactor))'; %resample at 1kHz
                    ACInfo.SyncPulse{i} = interp1qr(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):end)-min(ACInfo.Timestamps(ACInfo.SessionStarts(i-1):ACInfo.SessionStarts(i))),...
                                                    ACInfo.Data(ACInfo.SessionStarts(i-1):ACInfo.SessionStarts(i),1),...
                                                    ACInfo.TimestampsDownsampled{i});    
                end
            end
        else
            ACInfo.TimestampsDownsampled = linspace(min(ACInfo.Timestamps)-min(ACInfo.Timestamps),...
                                                       max(ACInfo.Timestamps)-min(ACInfo.Timestamps),...
                                                       round(length(ACInfo.Timestamps)/ACInfo.SamplingFactor))'; %resample at 1kHz
            ACInfo.SyncPulse = interp1qr(ACInfo.Timestamps-min(ACInfo.Timestamps),...
                                            ACInfo.Data(:,1),...
                                            ACInfo.TimestampsDownsampled);
        end
        
%         if iscell(ACInfo.SyncPulse)
%             for i = 1:length(ACInfo.SyncPulse)
%                 SPDiff = diff(ACInfo.SyncPulse{i});
%                 % find rising edges
%                 [temp_ppks,temp_plocs] = findpeaks(SPDiff);
%                 ppks = temp_ppks(temp_ppks>2); % select real edges [Volts]
%                 plocs = temp_plocs(temp_ppks>2);
%                 clear temp_ppks temp_plocs
%                 % find falling edges
%                 [temp_npks,temp_nlocs] = findpeaks(-SPDiff);
%                 npks = temp_npks(temp_npks>2); % select real edges [Volts]
%                 nlocs = temp_nlocs(temp_npks>2);
%                 clear temp_npks temp_nlocs
%                  % verify quality of pks finding
% %                         figure
% %                         plot(ACInfo.SyncPulse{i},'k')
% %                         hold on
% %                         plot(SPDiff,'r')
% %                         hold on
% %                         plot(plocs,ppks,'o')
% %                         hold on
% %                         plot(nlocs,-npks,'o')
%                 % find if the first edge is rising or falling
%                 if min(plocs)<min(nlocs)
%                     % compute the difference between the second rising edge and the
%                     % first falling edge
%                     SyncPulseIntervals(:,1) = plocs(2:end)-nlocs(1:length(plocs(2:end)));
%                     SyncPulseIntervals(:,2) = plocs(2:end);
%                     % check if any of the values is bigger than 3.5 seconds, i.e. there
%                     % is an error in the steps above and stops the program
%                     if ~isempty(find(SyncPulseIntervals(:,1)>3500))
%                         fprintf('Error computing sync pulse. Please check it manually \n')
%                         return
%                     end
%                 else
%                     % compute the difference between the first rising edge and the
%                     % first falling edge
%                     SyncPulseIntervals(:,1) = plocs(1:end)-nlocs(1:length(plocs(1:end)));
%                     SyncPulseIntervals(:,2) = plocs(1:end);
%                     if ~isempty(find(SyncPulseIntervals(:,1)>3500))
%                         fprintf('Error computing sync pulse. Please check it manually \n');
%                         return
%                     end
%                 end
%                 clear nlocs plocs ppks nks
%                 % construct an array as long as the timeseries (1kHz) with the
%                 % inter-pulse intervals time
%                 temp_SyncPulseTimeSeries = zeros(length(ACInfo.TimestampsDownsampled{i}),1);
%                 for sp_i=1:length(SyncPulseIntervals(:,1))
%                     if sp_i==1
%                         temp_SyncPulseTimeSeries(1:SyncPulseIntervals(sp_i,2)) = SyncPulseIntervals(sp_i,1);
%                     else
%                         temp_SyncPulseTimeSeries(SyncPulseIntervals(sp_i-1,2)+1:SyncPulseIntervals(sp_i,2)) = SyncPulseIntervals(sp_i,1);
%                     end
%                 end
%                 ACInfo.SyncPulseTimeSeries{i} = temp_SyncPulseTimeSeries;
%                 clear SyncPulseIntervals temp_SyncPulseTimeSeries plocs nlocs
%                 % verify timeseries
% %                 figure
% %                 plot(ACInfo.TimestampsDownsampled{i},ACInfo.SyncPulseTimeSeries{i}/max(ACInfo.SyncPulseTimeSeries{i}))
% %                 hold on
% %                 plot(ACInfo.TimestampsDownsampled{i},ACInfo.SyncPulse{i}/max(ACInfo.SyncPulse{i}))
%             end
%         else
%             SPDiff = diff(ACInfo.SyncPulse);
%             % find rising edges
%             [temp_ppks,temp_plocs] = findpeaks(SPDiff);
%             ppks = temp_ppks(temp_ppks>2); % select real edges [Volts]
%             plocs = temp_plocs(temp_ppks>2);
%             clear temp_ppks temp_plocs
%             % find falling edges
%             [temp_npks,temp_nlocs] = findpeaks(-SPDiff);
%             npks = temp_npks(temp_npks>2); % select real edges [Volts]
%             nlocs = temp_nlocs(temp_npks>2);
%             clear temp_npks temp_nlocs
%             % verify quality of pks finding
%             %             figure
%             %             plot(ACInfo.SyncPulse,'k')
%             %             hold on
%             %             plot(SPDiff,'r')
%             %             hold on
%             %             plot(plocs,ppks,'o')
%             %             hold on
%             %             plot(nlocs,-npks,'o')
%             % find if the first edge is rising or falling
%             if min(plocs)<min(nlocs)
%                 % compute the difference between the second rising edge and the
%                 % first falling edge
%                 SyncPulseIntervals(:,1) = plocs(2:end)-nlocs(1:length(plocs(2:end)));
%                 SyncPulseIntervals(:,2) = plocs(2:end);
%                 % check if any of the values is bigger than 3.5 seconds, i.e. there
%                 % is an error in the steps above and stops the program
%                 if ~isempty(find(SyncPulseIntervals(:,1)>3500))
%                     fprintf('Error computing sync pulse. Please check it manually \n')
%                     return
%                 end
%             else
%                 % compute the difference between the first rising edge and the
%                 % first falling edge
%                 SyncPulseIntervals(:,1) = plocs(1:end)-nlocs(1:length(plocs(1:end)));
%                 SyncPulseIntervals(:,2) = plocs(1:end);
%                 if ~isempty(find(SyncPulseIntervals(:,1)>3500))
%                     fprintf('Error computing sync pulse. Please check it manually \n');
%                     return
%                 end
%             end
%             clear nlocs plocs ppks nks
%             % construct an array as long as the timeseries (1kHz) with the
%             % inter-pulse intervals time
%             temp_SyncPulseTimeSeries = zeros(length(ACInfo.TimestampsDownsampled),1);
%             for sp_i=1:length(SyncPulseIntervals(:,1))
%                 if sp_i==1
%                     temp_SyncPulseTimeSeries(1:SyncPulseIntervals(sp_i,2)) = SyncPulseIntervals(sp_i,1);
%                 else
%                     temp_SyncPulseTimeSeries(SyncPulseIntervals(sp_i-1,2)+1:SyncPulseIntervals(sp_i,2)) = SyncPulseIntervals(sp_i,1);
%                 end
%             end
%             ACInfo.SyncPulseTimeSeries = temp_SyncPulseTimeSeries;
%             % verify timeseries
%             %             figure
%             %             plot(ACInfo.SyncPulseTimeSeries)
%         end
        
    end

end