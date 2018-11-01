function res = preprocessLFP(animal, iseries, iexplist, Probe, F_VS)
% dwonsample the LFP to 1kHz and save it to a mat file
SetDirs;
global pepNEV
global DIRS    

FcomputeCSD = false;
if nargin < 5
    F_VS = false;
end

addpath(genpath('..\Behaviour analysis\'))
addpath(genpath('..\General\'))
addpath(genpath('..\Spike Analysis\'))
if nargin < 3
    iexplist = [];
end
if nargin < 4
    Probe = 'CA1';
end
switch Probe
    case 'V1'
        igroup = 0;
        suffix = {'V1'};
    case'CA1'
        if strcmp(animal,'M160114C_BALL') && iseries == 323
            igroup = [0:15];
            suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
        else
            igroup = [0:7];
            suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
        end
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
    
    bestCh = [];
    res(s).LFPfilt = [];
    res(s).LFPthetaGlobal = [];
    res(s).CSDfilt = [];
    res(s).chans = [];
    res(s).spikeTrain = [];
    res(s).cellinfo = [];
    res(s).mua = [];
    res(s).LFPtimestamps = [];
    res(s).Behtimestamps = [];
    res(s).mua_dw = [];
    res(s).traj = [];
    res(s).smthBallSpd = [];
    res(s).outcome = [];
    res(s).gain = [];
    res(s).blanks = [];
        
    for i = 1:numel(iexplist)
        %try
        iexp = iexplist(i);
        exptInfo.iexp    = iexp;
        
        %% Load the data file
        expName = [exptInfo.animal '_' num2str(exptInfo.iseries) '_' num2str(exptInfo.iexp) '.ns5'];
        ns5file = [DIRS.Cerebus filesep exptInfo.animal filesep num2str(exptInfo.iseries) filesep num2str(exptInfo.iexp) filesep expName];
        if exist(ns5file,'file')
            disp(ns5file);
            %     es = getVRspikes(animal,iseries,iexplist(1),1,100,1,0,0:7,'CA1',true);
            fname = [animal '_' num2str(iseries(s)) '_' num2str(iexp)];
            dDIRname = [DIRS.multichanspikes filesep animal filesep num2str(iseries(s))];
            try
                load([dDIRname filesep fname '_screenTimes']);
                Fskip = false;
            catch
                warning(['no screentimes for ' num2str(exptInfo.iexp)]);
                Fskip = true;
            end            
            
            if ~Fskip
                [~,exptInfo.SamplingRateInKHZ,exptInfo.nchan] = nsopen(ns5file);%nsopen2(ns5file);
                
                [chan_list] = MichiganGetLayout(animal, iseries(s));
                
                exptInfo.downsample = 1000;
                nbCh_probe = [32 32];
                
                if strcmp(Probe,'V1')
                    nbCh = nbCh_probe(2);
                    chan_list = chan_list((nbCh_probe(1)+1):(nbCh_probe(1)+nbCh_probe(2)));
                    chgroup = (1:numel(chan_list))';
                    Espacing=20e-6;
                    el_pos=(Espacing:Espacing:numel(chan_list)*Espacing)';
                    el_num = (1:nbCh)';
                elseif strcmp(Probe,'CA1')
                    nbCh = nbCh_probe(1);
                    chan_list = chan_list(1:nbCh_probe(1));
                    Espacing = 17e-6;
                    gap = 500e-6;
                    Pmap = [Espacing;2*Espacing;2*Espacing;3*Espacing;gap + Espacing;gap + 2*Espacing;gap + 2*Espacing;gap + 3*Espacing];
                    el_pos = repmat(Pmap,[1 4]);
                    el_num = [(1:8)' (9:16)' (17:24)' (25:32)'];
                end
                
                fqrange  = [0.2 100];%[0.1 exptInfo.downsample/2*(1-0.05)];%[4 9];%
                thetafrange = [6 9];
                fprintf('Dowsampling Channel #: ');
                
                FloadVRfile = false;
                [es, spikeRate, VRdata, VRdata_o, chans] = getVRspikes(animal,iseries(s),iexp,1,100,1,0,igroup,suffix,false,50,60,FloadVRfile);
                
                if strcmp(Probe,'CA1')
                    res(s).CA1chref = 4*(mode(cell2mat(es.chanIDs))+1-1) + 1;
                end
                
                nevSamplingRateInKHZ = 30;
                LFPfilt_all = [];
                LFPshifted = [];
                corr_LFPmua = [];
                
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
                    
                    LFPfilt = LFPsignals(Chntemp, 1000, fqrange(1), fqrange(2));
                    LFPfilt.filtSig = Chntemp;
                    %         LFPphase(ch,:) = uint8(round(mod(360+180/pi*theta.phase,360)/2));
                    
                    procSpec.t = Timestamps.Chn;
                    procSpec.t = procSpec.t - screenTimes(1)./(1000*nevSamplingRateInKHZ);
                    Timestamps.Chn = Timestamps.Chn - screenTimes(1)./(exptInfo.downsample*nevSamplingRateInKHZ);
                    
                    
                    if isempty(LFPfilt_all)
                        LFPfilt_all = zeros(numel(LFPfilt.filtSig),nbCh);
                        %                     ThetaPhs_all = zeros(numel(LFPfilt.filtSig),nbCh);
                    end
                    LFPfilt_all(:,ch) = LFPfilt.filtSig(:);%interp1(Timestamps.Chn', theta.hill', es.sampleTimes);
                    
                    theta = LFPsignals(LFPfilt.filtSig, 1000, thetafrange(1), thetafrange(2));
                    
                    %                 thetaphs0 = round(mod(360 + 180/pi*angle(theta.hill'),360));
                    %                 thetapeak_idx0 = find([0;abs(diff(mod(thetaphs0,360)))>180]);
                    %                 shorttheta_lim = 90;%ms
                    %                 shorttheta = NaN;
                    %                 while ~isempty(shorttheta)
                    %                     thetaperiod0 = diff(thetapeak_idx0);
                    %                     shorttheta = find(thetaperiod0<shorttheta_lim);
                    %                     thetapeak_idx0(shorttheta + 1) = [];
                    %                 end
                    %                 thetaperiod0 = diff(thetapeak_idx0);
                    %                 thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
                    %
                    %                 thetaphase = zeros(size(thetaphs0));
                    %                 for tt = 1:numel(thetapeak_idx0)-1
                    %                     thetaphase(thetapeak_idx0(tt):thetapeak_idx0(tt+1)) = 360*(0:(thetapeak_idx0(tt+1)-thetapeak_idx0(tt)))/thetaperiod0(tt);
                    %                 end
                    %                 ThetaPhs_all(:,ch) = round(thetaphase);
                    
                    if isempty(LFPshifted)
                        LFPshifted = zeros(numel(real(theta.hill)),nbCh);
                    end
                    LFPshifted(:,ch) = real(theta.hill);
                end
                res(s).LFPfilt = cat(1,res(s).LFPfilt,LFPfilt_all);
                
                if strcmp(Probe,'CA1')
                    spkratemax = 10;
                    spkrate = sum(es.spikeTrain,1)/size(es.spikeTrain,1)*(1/nanmean(diff(es.sampleTimes)));
                    [N,~] = histcounts(cell2mat(es.chanIDs(spkrate < spkratemax)),0:8);
                    [~,pyrtet] = sort(N,'descend');
                    ch_pyrlayer = [];
                    for chtet = 1:4
                        ch_pyrlayer = [ch_pyrlayer (pyrtet(1:4)-1)*4 + chtet];
                    end
                    muapyr = sum(es.spikeTrain(:,spkrate < 10),2);
                    mua = interp1(es.sampleTimes, muapyr, Timestamps.Chn');
                    mua = mua - nanmean(mua);
                    maxlag = 80;
                    for ch = 1:nbCh
                        corr_LFPmua(ch,:) = xcorr(LFPshifted(~isnan(mua),ch) - nanmean(LFPshifted(~isnan(mua),ch)),mua(~isnan(mua)),maxlag,'coeff');%xcorr(real(theta),es.mua,maxlag,'coeff');
                        tshiftmax = find((corr_LFPmua(ch,:)) == max((corr_LFPmua(ch,:)))) - (maxlag+1);
                        tshiftmin = find((corr_LFPmua(ch,:)) == min((corr_LFPmua(ch,:)))) - (maxlag+1);
                        if abs(tshiftmax) < abs(tshiftmin)
                            tshift = tshiftmax;
                            signX = -1;
                        else
                            tshift = tshiftmin;
                            signX = +1;
                        end
                        %                 tshift = find((corr) == min((corr))) - (maxlag+1);
                        LFPshifted(:,ch) = signX * circshift(LFPshifted(:,ch),-tshift);%circshift(real(theta),-tshift);
                        tshiftLFP(ch) = tshift;
                    end
                    
                    theta = LFPsignals(mean(LFPshifted(:,ch_pyrlayer),2)', 1000, thetafrange(1), thetafrange(2));
                    theta = theta.hill';
                    res(s).LFPthetaGlobal = cat(1,res(s).LFPthetaGlobal,real(theta));
                end
                
                if strcmp(Probe,'V1') && FcomputeCSD
                    CSDfilt_All = zeros(size(LFPfilt_all));
                    Ntimes = size(CSDfilt_All,1);
                    Nwin = 100000;
                    Nseg = floor(Ntimes/Nwin)+1;
                    
                    for ishank = 1:size(el_pos,2)
                        for tk = 1:Nseg
                            CSDfilt_All(((tk-1)*Nwin + 1):min(tk*Nwin,Ntimes),el_num(:,ishank)) = (LFP2iCSD(LFPfilt_all(((tk-1)*Nwin + 1):min(tk*Nwin,Ntimes),el_num(:,ishank))', el_pos(:,ishank)))';
                        end
                    end
                    res(s).CSDfilt = cat(1,res(s).CSDfilt,CSDfilt_All);
                    
                    %                 for ch = 1:size(CSDfilt_All,2)
                    %                     theta = LFPsignals(CSDfilt_All(:,ch)', 1000, thetafrange(1), thetafrange(2));
                    %
                    %                     thetaphs0 = round(mod(360 + 180/pi*angle(theta.hill'),360));
                    %                     thetapeak_idx0 = find([0;abs(diff(mod(thetaphs0,360)))>180]);
                    %                     shorttheta_lim = 90;%ms
                    %                     shorttheta = NaN;
                    %                     while ~isempty(shorttheta)
                    %                         thetaperiod0 = diff(thetapeak_idx0);
                    %                         shorttheta = find(thetaperiod0<shorttheta_lim);
                    %                         thetapeak_idx0(shorttheta + 1) = [];
                    %                     end
                    %                     thetaperiod0 = diff(thetapeak_idx0);
                    %                     thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
                    %
                    %                     thetaphase = zeros(size(thetaphs0));
                    %                     for tt = 1:numel(thetapeak_idx0)-1
                    %                         thetaphase(thetapeak_idx0(tt):thetapeak_idx0(tt+1)) = 360*(0:(thetapeak_idx0(tt+1)-thetapeak_idx0(tt)))/thetaperiod0(tt);
                    %                     end
                    %
                    %                     ThetaPhs_all(:,ch) = round(thetaphase);
                    %                 end
                end
                LFPfilt_all = [];
                CSDfilt_All = [];
                
                res(s).spkrate = sum(es.spikeTrain,1)/size(es.spikeTrain,1)*(1/nanmean(diff(es.sampleTimes)));
                
                nbins = sum(Timestamps.Chn>0);
                spkstart = find(Timestamps.Chn>0,1,'first');
                idxoffset = numel(res(s).traj);
                binsize = 1/exptInfo.downsample;
                if isempty(res(s).chans)
                    res(s).chans = chans;
                    for ichan = 1:length(chans)
                        res(s).chans(ichan).spiketimes = [];
                        %                     res(s).chans(ichan).thetaphase = [];
                        res(s).Fgoodunit(ichan) = ~isempty(strfind(chans(ichan).id{5},'Good'));
                    end
                end
                for ichan = 1:length(chans)
                    st = ceil(chans(ichan).spiketimes./binsize);
                    numExcessSpikes = sum(st>nbins);
                    st(st>nbins) = [];
                    st(st==0) = [];
                    if numExcessSpikes>0
                        chans(ichan).spiketimes(end-numExcessSpikes+1:end) = [];
                    end
                    chans(ichan).spiketimes = ceil(chans(ichan).spiketimes./binsize)+spkstart;%*binsize;
                    %                 chans(ichan).thetaphase = ThetaPhs_all(chans(ichan).spiketimes,:);
                    try
                        res(s).chans(ichan).spiketimes = cat(1,res(s).chans(ichan).spiketimes,chans(ichan).spiketimes(:)+idxoffset);
                        %                 res(s).chans(ichan).thetaphase = cat(1,res(s).chans(ichan).thetaphase,chans(ichan).thetaphase);
                    catch
                        keyboard
                    end
                end
                
                %             [spk,chans] = getspiketrain(chans, 1/exptInfo.downsample,sum(Timestamps.Chn>0));
                %             spk.spikeTrain(ceil(es.sampleTimes(end)*exptInfo.downsample):end,:) = [];
                %             spknpts = size(spk.spikeTrain,1);
                %             spkncell = size(spk.spikeTrain,2);
                %             goodunits = false(1,spkncell);
                %             spktrain = [];
                %             for icell = 1:spkncell
                %                 goodunits(icell) = ~isempty(strfind(spk.spikeIDs{icell}{5},'Good'));
                %                 celltrain= interp1(1/exptInfo.downsample:1/exptInfo.downsample:(spknpts*1/exptInfo.downsample), spk.spikeTrain(:,icell), Timestamps.Chn','nearest');
                %                 if isempty(spktrain)
                %                     spktrain = zeros(numel(celltrain),spkncell);
                %                 end
                %                 spktrain(:, icell) = celltrain;
                %             end
                %             try
                %             mua = sum(spktrain(:,goodunits), 2);
                %             res(s).mua = cat(1, res(s).mua, mua(:));
                %             catch
                %                 keyboard
                %             end
                %
                %             res(s).spikeTrain = cat(1,res(s).spikeTrain,spktrain(:,goodunits));
                %             res(s).bestchan = spk.bestchan(goodunits);
                %             spktrain = [];
                
                %             res(s).nunits = hist(cell2mat(spk.bestchan(goodunits)),1:32);
                %
                %             mua_dw = interp1(es.sampleTimes, es.mua, Timestamps.Chn');
                %             res(s).mua_dw = cat(1,res(s).mua_dw,mua_dw(:));
                if isfield(es,'trajPercent')
                    traj = interp1(es.sampleTimes, es.trajPercent, Timestamps.Chn','nearest');
                    res(s).traj = cat(1,res(s).traj,traj(:));
                    
                    smthBallSpd = interp1(es.sampleTimes, es.smthBallSpd, Timestamps.Chn','nearest');
                    res(s).smthBallSpd = cat(1,res(s).smthBallSpd,smthBallSpd(:));
                    
                    outcome = interp1(es.sampleTimes, es.outcome, Timestamps.Chn','nearest');
                    res(s).outcome = cat(1,res(s).outcome,outcome(:));
                    
                    gain = interp1(es.sampleTimes, es.gain, Timestamps.Chn','nearest');
                    res(s).gain = cat(1,res(s).gain,gain(:));
                    
                    blanks = interp1(es.sampleTimes, single(es.blanks), Timestamps.Chn','nearest');
                    res(s).blanks = cat(1,res(s).blanks,blanks(:) == 1);
                    
                    if isempty(res(s).LFPtimestamps)
                        res(s).LFPtimestamps = Timestamps.Chn(:);
                        res(s).Behtimestamps = es.sampleTimes(:);
                    else
                        res(s).LFPtimestamps = cat(1, res(s).LFPtimestamps, res(s).LFPtimestamps(end) + Timestamps.Chn(:));
                        res(s).Behtimestamps = cat(1, res(s).Behtimestamps, res(s).Behtimestamps(end) + es.sampleTimes(:));
                    end
                else
                    res(s).traj = cat(1,res(s).traj,NaN(numel(Timestamps.Chn),1));
                end
                
                fprintf('\r');
            end
        end
    end
%     for icell = 1:size(res(s).spikeTrain,2)
%         res(s).chans(icell).spiketimes = find(res(s).spikeTrain(:,icell) >= 0.5)*1/exptInfo.downsample;
%     end
%     res(s).cellinfo = getCellInfo(res(s));
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

