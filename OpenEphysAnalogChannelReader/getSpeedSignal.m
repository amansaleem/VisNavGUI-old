%% Tomaso Muzzu - UCL - 10/04/2018

%% speed signal: save original signal downsampled at 1kHz, position in cm and speed in cm/s at 1kHz
function temp_ACInfo = getSpeedSignal(animal, iseries, iexp)
    
    global DIRS
    DIRname  = [DIRS.data filesep animal filesep num2str(iseries) filesep num2str(iexp)];
    % load the mat file containing the speed signalf
    
if (nargout<3) & exist([DIRname filesep 'AC_Info.mat'],'file')
        load([DIRname filesep 'AC_Info.mat']); 
        %%
        %%
        %% TO BE CONTINUED ........................
        %%
        %%
        % else load the analog channels now, find when the timestamps of the VS
        % onsets and offsets
    else
        
    temp_EncoderSignal1kHz(:,1) = resample(data(:,3),1,30);
    temp_EncoderSignal1kHz(:,2) = resample(data(:,4),1,30);
    % digitalise/normalise the signal
    temp_EncoderSignal1kHz(temp_EncoderSignal1kHz(:,1)>2.5,1) = 5; % originally it is a digital signal so it is 0-5V
    temp_EncoderSignal1kHz(temp_EncoderSignal1kHz(:,1)<=2.5,1) = 0; 
    temp_EncoderSignal1kHz(temp_EncoderSignal1kHz(:,2)>2.5,2) = 5; 
    temp_EncoderSignal1kHz(temp_EncoderSignal1kHz(:,2)<=2.5,2) = 0; 
    % save it
    temp_ACInfo.EncoderSignal1kHz = temp_EncoderSignal1kHz;
    % count numbers of rising edges 
    temp_DiffEncoderSignal1kHz(:,1) = diff(temp_ACInfo.EncoderSignal1kHz(:,1));
    temp_DiffEncoderSignal1kHz(:,2) = diff(temp_ACInfo.EncoderSignal1kHz(:,2));
    % find rising edges
    [temp_ppks,temp_plocs] = findpeaks( temp_DiffEncoderSignal1kHz(:,1));
    ppks_A = temp_ppks(temp_ppks>2.5); % select real edges [Volts]
    plocs_A = temp_plocs(temp_ppks>2.5);
    [temp_ppks,temp_plocs] = findpeaks( temp_DiffEncoderSignal1kHz(:,2));
    ppks_B = temp_ppks(temp_ppks>2.5); % select real edges [Volts]
    plocs_B = temp_plocs(temp_ppks>2.5);
    clear temp_ppks temp_plocs
    % find falling edges
    [temp_npks,temp_nlocs] = findpeaks(- temp_DiffEncoderSignal1kHz(:,1));
    npks_A = temp_npks(temp_npks>2.5); % select real edges [Volts]
    nlocs_A = temp_nlocs(temp_npks>2.5);
    [temp_npks,temp_nlocs] = findpeaks(- temp_DiffEncoderSignal1kHz(:,2));
    npks_B = temp_npks(temp_npks>2.5); % select real edges [Volts]
    nlocs_B = temp_nlocs(temp_npks>2.5);
    clear temp_npks temp_nlocs
    % 
    [v,ind] = min([min(plocs_A) min(nlocs_A) min(plocs_B) min(nlocs_B)]);
    temp_RFEvents = zeros(length(temp_ACInfo.EncoderSignal1kHz),1);
    temp_RFEvents(plocs_A) = 1; % if A is rising set value to 1
    temp_RFEvents(nlocs_A) = -1; % if A is falling set value to -1
    temp_RFEvents(plocs_B) = 2; % if B is rising set value to 2
    temp_RFEvents(nlocs_B) = -2; % if B is falling set value to -2
    for RF_i = 1:length(temp_RFEvents) 
        switch temp_RFEvents(RF_i)
            case 1
                if temp_ACInfo.EncoderSignal1kHz(RF_i,2)==0
                    temp_Position(RF_i)=1;
                elseif temp_ACInfo.EncoderSignal1kHz(RF_i,2)==5
                    temp_Position(RF_i)=-1;
                end
            case -1
                if temp_ACInfo.EncoderSignal1kHz(RF_i,2)==5
                    temp_Position(RF_i)=1;
                elseif temp_ACInfo.EncoderSignal1kHz(RF_i,2)==0
                    temp_Position(RF_i)=-1;
                end
            case 2
                if temp_ACInfo.EncoderSignal1kHz(RF_i,1)==5
                    temp_Position(RF_i)=1;
                elseif temp_ACInfo.EncoderSignal1kHz(RF_i,1)==0
                    temp_Position(RF_i)=-1;
                end
            case -2
                if temp_ACInfo.EncoderSignal1kHz(RF_i,1)==0
                    temp_Position(RF_i)=1;
                elseif temp_ACInfo.EncoderSignal1kHz(RF_i,1)==5
                    temp_Position(RF_i)=-1;
                end
            otherwise
                temp_Position(RF_i) = 0;
        end
    end         
    
    % verify speed signal
%     figure
%     axis([ACInfo.timestampsDownsampled(1) ACInfo.timestampsDownsampled(end) -10 10])
%     plot(ACInfo.timestampsDownsampled,temp_EncoderSignal1kHz(:,1))
%     hold on
%     plot(ACInfo.timestampsDownsampled,ACInfo.EncoderSignal1kHz(:,1),'k')
%     hold on
%     plot(ACInfo.timestampsDownsampled(2:end),temp_DiffEncoderSignal1kHz(:,1),'r')
%     hold on
%     plot(ACInfo.timestampsDownsampled(plocs_B), ppks_B,'ro')
%     hold on
%     plot(ACInfo.timestampsDownsampled(nlocs_B), npks_B,'ko')
%     hold on
%     plot(ACInfo.timestampsDownsampled,ACInfo.EncoderSignal1kHz(:,2),'c')
%     ylim([-6 6])
%     hold on
%     plot(ACInfo.timestampsDownsampled,temp_RFEvents,'g')
%     hold on 
%     plot(ACInfo.timestampsDownsampled,temp_Position,'y')
    
    temp_PositionSum = -cumsum(temp_Position); % signal A and B are inverted here.        
    % convert units into cm
    ResolEncoder = 1024;
    WheelRadius = 10;
    framerate = ((max(temp_ACInfo.timestampsDownsampled)-min(temp_ACInfo.timestampsDownsampled))/length(temp_ACInfo.timestampsDownsampled))^(-1);
    WheelCircum = WheelRadius*2*pi;
    temp_PositionSum = temp_PositionSum*(WheelCircum/ResolEncoder); % [cm]
    temp_ACInfo.PositionINcm = temp_PositionSum'; 
    % compute speed
    temp_speed = -temp_Position.*(WheelCircum/ResolEncoder).*framerate;
    % smooth speed signal to computer speed
    sigma = 50; % standard deviation in number of samples, in this case is 50ms
    Width = sigma*3; % convert size from seconds to number of samples
    x_g = linspace(-Width/2, Width/2, Width);
    gaussFilter = exp(-x_g.^2/(2*sigma^2));
    gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
    temp_speedSmoothed = conv(temp_speed, gaussFilter_, 'same');    
    temp_ACInfo.SpeedSpeedSmoothed = [temp_speed; temp_speedSmoothed]';
    temp_ACInfo.SamplingRateSignalsOut = framerate;
    % verify that position increases in general as the mice move forward mainly 
%     figure
%     plot(ACInfo.timestampsDownsampled,temp_PositionSum)
%     hold on
%     %plot(ACInfo.timestampsDownsampled,temp_speed)
%     hold on
%     plot(ACInfo.timestampsDownsampled,temp_speedSmoothed,'k')  
end