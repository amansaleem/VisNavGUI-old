function [es, VRdata, VRdata_o, chans] = getVR2pdata(animal,iseries,iexp,irecexp, ...
    discreet_steps, flag_load, igroup, Nplane)

%% initialization
global DIRS
global THRES

SetDirs


if nargin~=6
    flag_load = 1;
end

if nargin<5
    discreet_steps = 100;
end

if nargin<7
    flag_spont = 1;
end

if nargin<8
    flag_spkrate = 0;
end

if nargin < 10
    Nplane = 4;
end


fname = [animal '_' num2str(iseries) '_' num2str(iexp)];
dDIRname = [DIRS.data2p filesep animal filesep num2str(iseries)];

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
        display('No existing file');
        run_rest = 1;
    end
else
    run_rest = 1;
end



if run_rest
    
    load([dDIRname filesep fname '_screenTimes']);
    
    [VRdata, VRdata_o, es_test] = VRWheelLoad(animal, iseries, iexp);   
    
    screenTimes = screenTimes./1000;
    neuralFrameTimes = neuralFrameTimes./1000;
    %  Correcting and matching the recording and the VR
    RecToVR_correction = screenTimes(1) - es_test.screenTimes2(2);
    screenTimes = screenTimes - RecToVR_correction;
    neuralFrameTimes = neuralFrameTimes - RecToVR_correction;
    
    %% To get the even sampled parameters
    SampleRate = 60;
    evenSampleTime = (1/SampleRate):(1/SampleRate):max(screenTimes);
   
    es = VREvenSample(VRdata, evenSampleTime, screenTimes);
    es.iexp = zeros(size(es.sampleTimes));
    es.iexp(:) = iexp;
    es.series = ones(size(es.sampleTimes));
    
    es.isolDist = [];
    
    TimelineDir = fullfile(DIRS.expInfo,animal,num2str(iseries));
    
    NeuralRate = 1/nanmean(diff(neuralFrameTimes(2:Nplane:end)));
    highFq = 0.95*0.5*NeuralRate;
    lowFq = 0.05;
    
    es.spikeTrain = [];
    icell = 0;
    FgetSpikes = true;
    for igroup_idx = 1:length(igroup)
        if exist([dDIRname filesep 'F_' animal '_' num2str(iseries) '_plane' num2str(igroup(igroup_idx)) '_proc.mat'],'file')
            load([dDIRname filesep 'F_' animal '_' num2str(iseries) '_plane' num2str(igroup(igroup_idx)) '_proc.mat']);
            if isfield(dat,'sp') && FgetSpikes
                disp('Using deconvolved spikes');
                raw2pdata = dat.sp{irecexp};
            else
                if isfield(dat.stat,'neuropilCoefficient')
                    disp('Using raw fluorescence (with neuropil correction)');
                    raw2pdata = dat.Fcell{irecexp} - (dat.stat.neuropilCoefficient.*dat.FcellNeu{irecexp});
                else
                    disp('Using raw fluorescence (no neuropil correction)');
                    raw2pdata = dat.Fcell{irecexp};
                end
            end

            Fiscell = [dat.stat(:).iscell]>0;
            raw2pdata = raw2pdata(Fiscell,:);
            numCells = size(raw2pdata,1);

            if numel(neuralFrameTimes)~= 4*size(raw2pdata,2)
                warning(['neuralFrameTimes has ' num2str(numel(neuralFrameTimes)) ' frames ~= ' num2str(4*size(raw2pdata,2)) ' data points']);
                neuralFrameTimes = neuralFrameTimes(1:min(floor(end/4)*4,size(raw2pdata,2)*4));
                raw2pdata = raw2pdata(:,1:min(end,numel(neuralFrameTimes)/4));
            end
            if isempty(es.spikeTrain)
                es.spikeTrain = zeros(size(es.sampleTimes,1),1);
            end

            for n = 1:numCells
                icell = icell + 1;
                labelType = 'good';
                es.spikeTrain(:,icell) = interp1(neuralFrameTimes((igroup_idx):Nplane:end), raw2pdata(n,:), es.sampleTimes, 'linear');
                es.spikeTrain(isnan(es.spikeTrain(:,icell)),icell) = 0;
                if isfield(dat,'sp') && FgetSpikes
                    [pks,locs] = findpeaks(es.spikeTrain(:,icell));
                    es.spikeTrain(:,icell) = 0*es.spikeTrain(:,icell);
                    es.spikeTrain(locs,icell) = 1;
                else
                    es.spikeTrain(:,icell) = es.spikeTrain(:,icell) - min(es.spikeTrain(:,icell));
                end
                es.spikeIDs{icell} = ['plane#' num2str(igroup) '_' labelType '_cluster' num2str(icell)];
                es.chanIDs{icell} = igroup(igroup_idx);
                es.ProbeIDs{icell} = '1';
            end
        end
    end
    es.mua = sum(es.spikeTrain,2);
    
    %% Checking the the sync between the display and recording
    if length(es.screenTimes) ~= length(es.screenTimes2)-1
        disp(['!!!!!WARNING!!!!!!: ' num2str(animal) '_' num2str(iseries) '_' num2str(iexp) ' trajectory has ' num2str(length(es.screenTimes) - length(es.screenTimes2)+1) ' more frames than screen refreshes']);
    end
    
%     %% save them to expt
%     if exist([dDIRname filesep fname '_contspikes' '.mat'],'file')
%         display('Appending the new structures to the file!');
% %         save([dDIRname filesep fname '_contspikes'],'es', '-append');
%     else
%         display('Saving the file!');
%         display('Problem with append, saving without append!')
%         save([dDIRname filesep fname '_contspikes'],'VRdata','VRdata_o','es');
%     end
%     save([dDIRname filesep fname '_es'],'es');
end

