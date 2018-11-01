function res = preprocessLFP_VS(animal, iseries, iexplist, Probe)
% dwonsample the LFP to 1kHz and save it to a mat file
SetDirs;
global pepNEV
global DIRS    

addpath(genpath('..\Behaviour analysis\'))
addpath(genpath('..\General\'))
addpath(genpath('..\Spike Analysis\'))
if nargin < 3
    iexplist = [];
end
if nargin < 4
    Probe = 'CA1';
end
[chan_list] = MichiganGetLayout(animal, iseries);
switch Probe
    case 'V1'
        igroup = 0;
        el_pos = 2e-6:2e-6:(32*2e-6);
        if numel(chan_list)>32+5
            chan_list = chan_list((32+1):(32+32));
        end
        nbCh = 32;
    case'CA1'
        igroup = 0:7;
        chan_list = chan_list(1:32);
        nbCh = 32;
end
for s = 1:numel(iseries)
    if isempty(iexplist)
        infoAll = getDataInfo(animal);
        whichdate = num2str(iseries(s));
        iexplist = [];
        for i=1:numel(infoAll)
            if strcmpi(infoAll(i).date,whichdate) == 1
                iexplist = infoAll(i).sessions;
                break;
            end
        end
        if isempty(iexplist)
            keyboard
        end
    end
    exptInfo.animal  = animal;
    exptInfo.iseries = iseries(s);
    
    fqrange = [0.5 5];
    for iprotocol = 1:10
        VSexpName = [exptInfo.animal '_' num2str(exptInfo.iseries) '_' num2str(iprotocol) '.ns5'];
        ns5file = [DIRS.Cerebus filesep exptInfo.animal filesep num2str(exptInfo.iseries) filesep num2str(iprotocol) filesep VSexpName];
        
        
        exptInfo.downsample = 1000;
        
        if exist(ns5file)
            disp(ns5file)
            [~,exptInfo.SamplingRateInKHZ,exptInfo.nchan] = nsopen(ns5file);
            Timestamps.fulllength = length(pepNEV.ns.Data.data(chan_list(1),:));
            Timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;
            try
                [stimTimes,p] = createStimMatrix(exptInfo.animal, exptInfo.iseries, iprotocol);
                stimTimes.on = floor((stimTimes.on/Timestamps.samplingRateOrig)*exptInfo.downsample);
                stimTimes.off = floor((stimTimes.off/Timestamps.samplingRateOrig)*exptInfo.downsample);
                numStim = size(stimTimes.on,1);
                numReps = size(stimTimes.on,2);
                if numReps>=10
                    nbinPreStim = 500;
                    nbinPostStim = 1000;
                    if strfind(p.description,'9 y1')
                        res.StimType{iprotocol} = 'ypos';
                    elseif strfind(p.description,'10 x1')
                        res.StimType{iprotocol} = 'xpos';
                    elseif strfind(p.description,'12 ori')
                        res.StimType{iprotocol} = 'ori';
                    else
                        res.StimType{iprotocol} = 'other';
                    end
                    [~,exptInfo.SamplingRateInKHZ,exptInfo.nchan] = nsopen(ns5file);
                    res.LFPfilt_VS{iprotocol} = zeros(numStim,nbinPreStim + nbinPostStim + 1,nbCh,numReps);
                    LFPfilt = [];
                    for ch = 1:nbCh
                        %% Load the particular channels
                        fprintf([num2str(ch) ' ']);
                        Timestamps.fulllength = length(pepNEV.ns.Data.data(chan_list(ch),:));
                        Timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;
                        
                        % decimate(double(pepNEV.ns.Data.data(exptInfo.Chn(ichn).ChnLoc,:)), Timestamps.samplingRateOrig/exptInfo.downsample);
                        Chntemp = decimate(double(pepNEV.ns.Data.data(chan_list(ch),:)), Timestamps.samplingRateOrig/exptInfo.downsample);
                        Timestamps.ChnsampInt = 1./exptInfo.downsample;
                        Timestamps.Chnorig = 0:(1./Timestamps.samplingRateOrig):(Timestamps.fulllength./Timestamps.samplingRateOrig);
                        Timestamps.Chn = 0:Timestamps.ChnsampInt:...
                            (Timestamps.ChnsampInt*(Timestamps.fulllength*(exptInfo.downsample/Timestamps.samplingRateOrig)));
                        if numel(Timestamps.Chn) > numel(Chntemp)
                            Timestamps.Chn = Timestamps.Chn(1:numel(Chntemp));
                        end
                        LFPfilt.filtSig = LFPfilter(Chntemp, fqrange(1), fqrange(2), 1000);%LFPsignals(Chntemp, 1000, fqrange(1), fqrange(2));
                        procSpec.t = Timestamps.Chn;
                        
                        for istim = 1:numStim
                            nstim = 0;
                            for irep = 1:numReps
                                startTime = stimTimes.on(istim,irep);
                                endTime   = stimTimes.off(istim,irep);
                                res.LFPfilt_VS{iprotocol}(istim,:,ch,irep) = squeeze(res.LFPfilt_VS{iprotocol}(istim,:,ch,irep)) + LFPfilt.filtSig((startTime - nbinPreStim):(startTime + nbinPostStim));
                                nstim = nstim + 1;
                            end
                            %                     res(s).LFPfilt_VS{iprotocol}(istim,:,ch) = res(s).LFPfilt_VS{iprotocol}(istim,:,ch)/nstim;
                        end
                    end
                    res.LFPfilt_VS{iprotocol} = res.LFPfilt_VS{iprotocol} - repmat(median(res.LFPfilt_VS{iprotocol},3),[1 1 32 1]);
                    
                    res.LFPfilt_VS{iprotocol} = squeeze(mean(mean(res.LFPfilt_VS{iprotocol}(~ismember(1:p.nstim,p.blankstims),:,:,:),4),1))';
                    if strcmp(Probe,'V1')
                        res.CSDfilt_VS{iprotocol} = LFP2iCSD(res.LFPfilt_VS{iprotocol},el_pos);
                    else
                        res.CSDfilt_VS{iprotocol} = [];
                    end
                else
                    res.LFPfilt_VS{iprotocol} = [];
                    res.CSDfilt_VS{iprotocol} = [];
                    res.StimType{iprotocol} = [];
                end
%                 res.CSDfilt_VS{iprotocol} = zeros(size(res.LFPfilt_VS{iprotocol}));
%                 if strcmp(Probe,'V1')
%                     for istim = 1:numStim
%                         res.CSDfilt_VS{iprotocol}(istim,:,:) = squeeze(LFP2iCSD(squeeze(res.LFPfilt_VS{iprotocol}(istim,:,:))', el_pos))';
%                     end
%                 end
            catch
                disp(['pb with ' ns5file])
                res.LFPfilt_VS{iprotocol} = [];
                res.CSDfilt_VS{iprotocol} = [];
                res.StimType{iprotocol} = [];
            end
        else
            res.LFPfilt_VS{iprotocol} = [];
            res.CSDfilt_VS{iprotocol} = [];
            res.StimType{iprotocol} = [];
        end
    end
end

clear global pepNEV


end

function [spk,Chans] = getspiketrain(chans,binsize, nbins)
numCells = length(chans);
spk.spikeTrain = zeros(nbins, numCells);
for ichan = 1:length(chans)
    st = ceil(chans(ichan).spiketimes./binsize);
    numExcessSpikes = sum(st>nbins);
    st(st>nbins) = [];
    st(st==0) = [];
    if numExcessSpikes>0
        chans(ichan).spiketimes(end-numExcessSpikes+1:end) = [];
        %             chans(ichan).spikephase(end-numExcessSpikes+1:end) = [];
    end
    chans(ichan).spiketimes = ceil(chans(ichan).spiketimes./binsize)*binsize;
    spk.spikeTrain(st,ichan) = 1;
    
    repeats = st(diff(st)==0);
    for irep = 1:length(repeats)
        spk.spikeTrain(repeats(irep),ichan) = spk.spikeTrain(repeats(irep),ichan) + 1;
        %need to switch to circular coordinates to average the spike phase
        %             es.spikePhase(repeats(irep),ichan) = es.spikePhase(repeats(irep),ichan) + es.spikePhase(repeats(irep)+1,ichan);
    end
    spk.spikeIDs{ichan} = chans(ichan).id;
    spk.chanIDs{ichan} = chans(ichan).ichan;
    spk.bestchan{ichan} = chans(ichan).bestchan;
    spk.ProbeIDs{ichan} = chans(ichan).ProbeIDs;
    spk.spikeTimes{ichan} = chans(ichan).spiketimes;
    Chans(ichan).spiketimes = chans(ichan).spiketimes(:);
end
end

function cellinfo = getCellInfo(res)
spdth = 5;
FCircularMaze = true;

cellinfo.field = [];

tidx = res.gain == 0.5 & res.outcome == 2 & ~res.blanks & res.smthBallSpd > spdth & ~isnan(res.smthBallSpd) & res.traj <= 100;
varX = res.traj(:);
dx = 1;
nXbins = 100;
samplerate = 1000;
Xbinsize = 4;
nbXbinsmth = round(1/(Xbinsize/100));

nCells = size(res.spikeTrain,2);
cellinfo.field = zeros(nXbins,nCells);
cellinfo.fieldmax = zeros(nCells,1);

for icell = 1:nCells
    [map,~] = fast1Dmap(varX(tidx & ~isnan(res.spikeTrain(:,icell))), res.spikeTrain(tidx & ~isnan(res.spikeTrain(:,icell)),icell), dx, samplerate,nbXbinsmth,FCircularMaze);
    cellinfo.field(:,icell) = map(:);
    [~, imax] = max(map);
    cellinfo.fieldmax(icell) = imax;
end
end

