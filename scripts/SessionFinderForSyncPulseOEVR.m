

function [VR_SessionOI ind_M] = SessionFinderForSyncPulseOEVR(es, AC_Info)

    % OE sync pulse signal
    OE_SyncPulse = AC_Info.SyncPulse; 
    OE_SyncPulse(OE_SyncPulse<1) = 0;
    OE_SyncPulse(OE_SyncPulse>1) = 1;
    OE_TimeStamps = AC_Info.TimestampsDownsampled;
    OE_SamplingRate = AC_Info.SamplingRateOE/AC_Info.SamplingFactor;
    %% VR sync pulse signal
    VR_SyncPulse = es.SyncPulse;
    VR_TimeStamps = es.sampleTimes;
    VR_SamplingRate = str2num(es.sampleRate(1:2));
    
%     figure
%     plot(TimeStamps_OE,SyncPulse_OE)
%     ylabel('OE - 120 Hz'); xlabel('seconds');
%     figure
%     plot(TimeStamps_VR,SyncPulse_VR)
%     ylabel('VR - 60 Hz'); xlabel('seconds');
    
    % select a chunk of signal from the middle, long as ComparisonTime
    ComparisonTime = 20; % seconds ( not the number of samples!!! )
    ChunkStart = 0.5;
    clear OE_chunk OE_chunkTime
    OE_chunk = OE_SyncPulse(round(ChunkStart*length(OE_SyncPulse)):round(ChunkStart*length(OE_SyncPulse)+OE_SamplingRate*ComparisonTime));
    OE_chunkTime = OE_TimeStamps(round(ChunkStart*length(OE_TimeStamps)):round(ChunkStart*length(OE_TimeStamps)+OE_SamplingRate*ComparisonTime)) - ...
                   min(OE_TimeStamps(round(ChunkStart*length(OE_TimeStamps)):round(ChunkStart*length(OE_TimeStamps)+OE_SamplingRate*ComparisonTime)));
    figure
    plot(OE_chunkTime,OE_chunk)
    
    % downsample OE_chunk to 60Hz
    OE_chunk_DS = interp1(OE_chunkTime,...
                            OE_chunk, ...
                            1/VR_SamplingRate:1/VR_SamplingRate:ComparisonTime);
    OE_chunkTime_DS = interp1(OE_chunkTime,...
                            OE_chunkTime, ...
                            1/VR_SamplingRate:1/VR_SamplingRate:ComparisonTime);
    % scan the entire sync pulse signal recorded from the VR and find the
    % best match.
    clear k
    for i = 1:length(VR_SyncPulse)-(VR_SamplingRate*ComparisonTime)-1 % add 200 samples of buffering for the end
%         [v,indm] = min(abs(TimeStamps_VR(i+1:end)-(TimeStamps_VR(i+1)+ComparisonTime)));
%         VR_ChunkTime = TimeStamps_VR(i+1:i+indm)-min(TimeStamps_VR(i+1:i+indm));
%         VR_SyncPulse_res = interp1qr(VR_ChunkTime,...
%                                    SyncPulse_VR(i+1:i+indm), ...
%                                    (1/SamplingRate_VR:1/SamplingRate_VR:ComparisonTime)');
%         R = corrcoef(VR_SyncPulse_res,OE_chunk_DS); % 500 is duration of pulse
        R = corrcoef(VR_SyncPulse(i+1:i+ComparisonTime*VR_SamplingRate),OE_chunk_DS); % 500 is duration of pulse
        k(i) = R(1,2);
        clear R;
    end
    figure; plot(k,'.')  
    [val_M, ind_M] = max(k);

    plot(OE_chunkTime_DS,OE_chunk_DS,'k')
    hold on
    plot(VR_TimeStamps(ind_M:ind_M+ComparisonTime*VR_SamplingRate)-min(VR_TimeStamps(ind_M:ind_M+ComparisonTime*VR_SamplingRate)),...
        VR_SyncPulse(ind_M:ind_M+ComparisonTime*VR_SamplingRate),'r.')
    
    % Determine which VR session was recorded from
    VR_SessionOI = es.iexp(ind_M)
        
end







