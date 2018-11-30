function [es, spikeRate, VRdata, VRdata_o, chans, spontRate, spikeTimes] = ...
    getVRspikes(animal,iseries,iexp, ...
    discreet_steps, flag_load, flag_spont, flag_spkrate, igroup,iaddinfo,loadTimes,SmthTimeWindow,samplerate,FloadVRfile)

%% initialization
global DIRS
global THRES

% SetDirs
SetDefaultDirs2018

if nargin<5
    flag_load = 1;
end

if nargin<4
    discreet_steps = 100;
end

if nargin<6
    flag_spont = 1;
end

if nargin<7
    flag_spkrate = 0;
end

if nargin<10
    loadTimes = 0;
end

if nargin<11
    SmthTimeWindow = 150;
end

if nargin<12
    samplerate = 60;
end

if nargin<13
    FloadVRfile = true;
end

if numel(iaddinfo) == 1 && numel(igroup) > 1
    iaddinfotemp = iaddinfo{1};
    iaddinfo = cell(1,numel(igroup));
    for group_idx = 1:numel(igroup)
        iaddinfo{group_idx} = iaddinfotemp;
    end
end


fname = [animal '_' num2str(iseries) '_' num2str(iexp)];
%dDIRname = [DIRS.multichanspikes filesep animal filesep num2str(iseries)];
dDIRname = [DIRS.multichanspikes filesep num2str(iseries)];
%dDIRname2 = [DIRS.data2 filesep animal filesep num2str(iseries)];

if flag_load
    if exist([dDIRname filesep fname '_contspikes' '.mat'],'file')
        try
            load([dDIRname filesep fname '_contspikes','.mat'],'es');
%         if exist('es','var')
            if size(es.trajpos,2)==discreet_steps
                run_rest = 0;
            else
                display('Old structure has different position sampling')
                run_rest = 1;
            end
        catch %else
            display('No even sampling')
            run_rest = 1;
        end
    else
        display('No existing VR file');
        run_rest = 1;
    end
else
    run_rest = 1;
end



if run_rest
    
    load([dDIRname filesep fname '_screenTimes']);
    
    if FloadVRfile
        [VRdata, VRdata_o, es_test] = VRWheelLoad(animal, iseries, iexp, SmthTimeWindow);
    else
        display('won''t load VR file');
        VRdata = [];
        VRdata_o = [];
        es_test = [];
        es = [];
    end
    
    if nargin<9
        inp = inputdlg('Enter the group to be loaded:','0');
        igroup = str2num(inp{1});
    end
    %for loading kwik file clusters
    idx = 1;
    for igroup_idx = 1:length(igroup)
        if strcmp(iaddinfo{igroup_idx},'V1')
            SpikeSorterType = 'klusta';%'kilosort';%
        else
            SpikeSorterType = 'klusta';%'kilosort';%
        end
        try
            switch SpikeSorterType
                case 'klusta'
                    temp = getKwikSpikes(animal, iseries, iexp, igroup(igroup_idx),iaddinfo{igroup_idx});
                    tempNumCells = length(temp);
                    for k = 1:tempNumCells
                        temp(k).ProbeIDs = iaddinfo{igroup_idx};
                    end
                    if idx==1
                        chans=temp;
                        numCells = tempNumCells;
                    else
                        for icell = 1:tempNumCells
                            chans(numCells+1) = temp(icell);
                            numCells = numCells + 1;
                        end
                    end
                case 'kilosort'
                    temp = getKiloSortSpikes(animal, iseries, iexp, igroup(igroup_idx),iaddinfo{igroup_idx});
                    tempNumCells = length(temp);
                    for k = 1:tempNumCells
                        temp(k).ProbeIDs = iaddinfo{igroup_idx};
                    end
                    if idx==1
                        chans=temp;
                        numCells = tempNumCells;
                    else
                        for icell = 1:tempNumCells
                            chans(numCells+1) = temp(icell);
                            numCells = numCells + 1;
                        end
                    end
                case 'spikeconsensus'
                    temp = getConsensusSpikes(animal, iseries, iexp, igroup(igroup_idx),'CA1');
                    tempNumCells = length(temp);
                    numCells = tempNumCells;
                    chans=temp;
            end
            idx = idx + 1;
        catch
            warning(['no spikes sorted on shank ' num2str(igroup(igroup_idx))]);
        end
    end
    isolDist = zeros(1,length(chans));
    
    screenTimes = screenTimes./30000;
    %  Correcting and matching the recording and the VR
    if ~isempty(es_test)
        RecToVR_correction = screenTimes(1) - es_test.screenTimes2(2);
        screenTimes = screenTimes - RecToVR_correction;
    else
        es.screenTimes = screenTimes;
        es.screenTimes2 = screenTimes;
        RecToVR_correction = screenTimes(1);
        screenTimes = screenTimes - RecToVR_correction;
    end
    
    %% To get the even sampled parameters
    evenSampleTime = (1/samplerate):(1/samplerate):max(screenTimes);
    if ~isempty(VRdata)
        es = VREvenSample(VRdata, evenSampleTime, screenTimes, samplerate);
    else
        es.sampleTimes = evenSampleTime;
    end
    es.iexp = zeros(size(es.sampleTimes));
    es.iexp(:) = iexp;
    es.series = ones(size(es.sampleTimes));
    
    es.isolDist = isolDist;
    
    numCells = length(chans);
    
    if exist([dDIRname filesep fname '_LFP3.mat']) && FloadVRfile
        try
        nChCA1 = 34;
        LFPfile = matfile([dDIRname filesep fname '_LFP3.mat']);
        es.LFPphase = mod(double(LFPfile.LFPthetaPhase(nChCA1,:))*2,360)';
        es.LFPpower = LFPfile.LFPthetaPower(nChCA1,:)';
        es.LFPtheta = LFPfile.LFPtheta(nChCA1,:)';
        SpkthetaPhase = LFPfile.SpkthetaPhase;
        if numel(SpkthetaPhase) == length(chans)
            for ichan = 1:length(chans)
                chans(ichan).spikephase = SpkthetaPhase{ichan}(:,nChCA1);
            end
        else
            warning('theta phase of spikes won''t be loaded')
        end
        catch
            warning('CA1 LFP must be preprocessed again')
        end
    else
        warning('theta phase of spikes won''t be loaded')
    end
    if exist([dDIRname filesep fname '_V1_LFP3.mat'])%([dDIRname2 filesep fname '_V1_LFP.mat'])
%         try
%         LFPfile = matfile([dDIRname filesep fname '_V1_LFP2.mat']);
%         es.LFPphaseV1 = mod(double(LFPfile.LFPthetaPhase)*2,360)';
%         es.LFPpowerV1 = LFPfile.LFPthetaPower';
%         es.LFPthetaV1 = LFPfile.LFPtheta';
%         SpkthetaPhase = LFPfile.SpkthetaPhase;
%         if numel(SpkthetaPhase) == length(chans)
%             for ichan = 1:length(chans)
%                 chans(ichan).spikephaseV1 = SpkthetaPhase{ichan};
%             end
%         else
%             warning('theta phase of spikes won''t be loaded')
%         end
%         catch
%             warning('V1 LFP must be preprocessed again')
%         end
    end  
    
    % getting the spike trains
    numCells = length(chans);
    es.spikeTrain = zeros(length(es.sampleTimes), numCells);
    if isfield(chans(1),'spikephase')
        es.spikePhase = NaN(length(es.sampleTimes), numCells);
    end
%     if isfield(chans(ichan),'spikephaseV1')
%         es.spikePhaseV1 = NaN(length(es.sampleTimes), numCells);
%     end
    for ichan = 1:length(chans)
        chans(ichan).spiketimes = chans(ichan).spiketimes - RecToVR_correction;
        
%         if isfield(chans(ichan),'spikephase')
%             chans(ichan).spikephase(chans(ichan).spiketimes<0) = [];
%         end

%         
%         if isfield(chans(ichan),'spikephaseV1')
%             chans(ichan).spikephaseV1(chans(ichan).spiketimes<0) = [];
%         end
        chans(ichan).spiketimes(chans(ichan).spiketimes<0) = [];
        
%         display(['Processing cell: ' num2str(ichan)]);
        % finding the bin for each spike
        st = ceil(chans(ichan).spiketimes./(1/samplerate));
        numExcessSpikes = sum(st>length(es.sampleTimes));
        st(st>length(es.sampleTimes)) = [];
        st(st==0) = [];
        
        if numExcessSpikes>0
%             if isfield(chans(ichan),'spikephase')
%                 chans(ichan).spikephase(end-numExcessSpikes+1:end) = [];
%             end
%             if isfield(chans(ichan),'spikephaseV1')
%                 chans(ichan).spikephaseV1(end-numExcessSpikes+1:end) = [];
%             end
            chans(ichan).spiketimes(end-numExcessSpikes+1:end) = [];
        end
        
        [ust, idx] = unique(st);
        es.spikeTrain(st,ichan) = 1;
        if isfield(chans(ichan),'spikephase')
            es.spikePhase(ust,ichan) = chans(ichan).spikephase(min(end,idx));
        end
%         if isfield(chans(ichan),'spikephaseV1')
%             es.spikePhaseV1(ust,ichan) = chans(ichan).spikephaseV1(idx);
%         end
        
        repeats = st(diff(st)==0);
        repspknum = find(diff(st)==0);
        for irep = 1:length(repeats)
            es.spikeTrain(repeats(irep),ichan) = es.spikeTrain(repeats(irep),ichan) + 1;
            if isfield(chans(ichan),'spikephase')
                es.spikePhase(repeats(irep),ichan) = mod(360/(2*pi)*circ_mean(2*pi/360*[es.spikePhase(repeats(irep),ichan);chans(ichan).spikephase(repspknum(irep)+1)]),360);
            end
%             if isfield(chans(ichan),'spikephaseV1')
%                 es.spikePhaseV1(repeats(irep),ichan) = mod(360/(2*pi)*circ_mean(2*pi/360*[es.spikePhaseV1(repeats(irep),ichan);es.spikePhaseV1(repeats(irep)+1,ichan)]),360);
%             end
        end
        es.spikeIDs{ichan} = chans(ichan).id;
        es.chanIDs{ichan} = chans(ichan).ichan;
        es.bestchan{ichan} = chans(ichan).bestchan;
        es.waveform{ichan} = chans(ichan).waveform;
        es.ProbeIDs{ichan} = chans(ichan).ProbeIDs;
        if loadTimes
            es.spikeTimes{ichan} = chans(ichan).spiketimes;
        end
        clear st trash st_unique order repeats irep
    end
    es.mua = sum(es.spikeTrain,2);
    %% Checking the the sync between the display and recording
    if length(es.screenTimes) ~= length(es.screenTimes2)-1
        disp(['!!!!!WARNING!!!!!!: ' num2str(animal) '_' num2str(iseries) '_' num2str(iexp) ' trajectory has ' num2str(length(es.screenTimes) - length(es.screenTimes2)+1) ' more frames than screen refreshes']);
    end
    
    
    %% get all spikeTimes
    if flag_spkrate
        for ichan = 1:length(chans)
            if isempty(chans(ichan).spiketimes)
                spikeTimes(ichan).t = [];
                spikeRate(ichan).t = [];
            else
                spikeTimes(ichan).t = chans(ichan).spiketimes;
                for iscreen = 2:length(screenTimes)
                    nspike = sum((spikeTimes(ichan).t > screenTimes(iscreen-1)) & (spikeTimes(ichan).t <= screenTimes(iscreen)));
                    spikeRate(ichan).t(iscreen-1) = nspike/(screenTimes(iscreen)-screenTimes(iscreen-1));
                end
            end            
        end
    else
        spikeRate = [];
    end
    %% the spont rate for the period before the start of the VR
    if flag_spont
        spontRate = zeros(1,length(chans));
        for cellID = 1:length(chans)
            numSpks(cellID) = sum(chans(cellID).spiketimes<screenTimes(1));
            
        end
        spontRate = numSpks./screenTimes(1);
        clear numSpks
    else
        spontRate = [];
    end
    
    %% save them to expt
%     if exist([dDIRname filesep fname '_contspikes' '.mat'],'file')
%         display('Appending the new structures to the file!');
%         if flag_spkrate
%             save([dDIRname filesep fname '_contspikes'],'spikeRate','spikeTimes','es', '-append');
%         else
%             save([dDIRname filesep fname '_contspikes'],'es', '-append');
%         end
%     else
%         display('Saving the file!');
%         display('Problem with append, saving without append!')
%         
%         if flag_spkrate
%             save([dDIRname filesep fname '_contspikes'],'spikeRate','VRdata','VRdata_o','chans','spikeTimes','es');
%         else
%             save([dDIRname filesep fname '_contspikes'],'VRdata','VRdata_o','chans','es');
%         end
%     end
%     save([dDIRname filesep fname '_es'],'es');
end

