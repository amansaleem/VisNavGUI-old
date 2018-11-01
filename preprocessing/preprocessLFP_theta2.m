function preprocessLFP_theta2(animal, iseries, iexplist, Probe)
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
    chref = [];
    for i = 1:numel(iexplist)
%         try
        iexp = iexplist(i);                
        exptInfo.animal  = animal;
        exptInfo.iseries = iseries(s);
        exptInfo.iexp    = iexp;
        
        %% Load the data file
        expName = [exptInfo.animal '_' num2str(exptInfo.iseries) '_' num2str(exptInfo.iexp) '.ns5'];
        ns5file = [DIRS.Cerebus filesep exptInfo.animal filesep num2str(exptInfo.iseries) filesep num2str(exptInfo.iexp) filesep expName]
        
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
            nbChCA1 = 32;
            nbChV1 = 32;
            nbCh = nbChCA1;%numel(chan_list);
            
            if strcmp(Probe,'V1')
                chan_list = chan_list(nbChCA1+1:nbChCA1+32);
            end
            
            LFPthetaPower = [];
            LFPthetaPhase = [];
            LFPtheta = [];
            thetafrange  = [6 9];%[4 9];%
            nevSamplingRateInKHZ = 30;
            fprintf('Dowsampling Channel #: ');
            
            Timestamps.fulllength = length(pepNEV.ns.Data.data(chan_list(1),:));
            Timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;
            Timestamps.ChnsampInt = 1./exptInfo.downsample;
            Timestamps.Chnorig = 0:(1./Timestamps.samplingRateOrig):(Timestamps.fulllength./Timestamps.samplingRateOrig);
            Timestamps.Chn = 0:Timestamps.ChnsampInt:...
                (Timestamps.ChnsampInt*(Timestamps.fulllength*(exptInfo.downsample/Timestamps.samplingRateOrig)));
            procSpec.t = Timestamps.Chn;
            procSpec.t = procSpec.t - screenTimes(1)./(1000*nevSamplingRateInKHZ);
            Timestamps.Chn = Timestamps.Chn - screenTimes(1)./(exptInfo.downsample*nevSamplingRateInKHZ);
            

            if strcmp(animal,'M160114C_BALL') && iseries(s) == 323
                shanknum = [0:15 0];
                suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
            else
                shanknum = [0:7 0];
                suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
            end
            [es, spikeRate, VRdata, VRdata_o, chans] = getVRspikes(animal,iseries(s),iexp,1,100,1,0,shanknum,suffix);
            
            probeID = zeros(1,numel(es.ProbeIDs));
            for icell = 1:numel(es.ProbeIDs)
                if strcmp(es.ProbeIDs{icell},'CA1')
                    probeID(icell) = 1;
                elseif strcmp(es.ProbeIDs{icell},'V1')
                    probeID(icell) = 2;
                end
            end
            if strcmp(Probe,'CA1')
                muacellidx = probeID == 1;
            elseif strcmp(Probe,'V1')
                muacellidx = probeID == 2;
            end
            
            nbins = sum(Timestamps.Chn>0);
            spkstart = find(Timestamps.Chn>0,1,'first');
            binsize = 1/exptInfo.downsample;
            for ichan = 1:length(chans)
%                 st = ceil(chans(ichan).spiketimes./binsize);
%                 numExcessSpikes = sum(st>nbins);
%                 st(st>nbins) = [];
%                 st(st==0) = [];
%                 if numExcessSpikes>0
%                     chans(ichan).spiketimes(end-numExcessSpikes+1:end) = [];
%                 end
                chans(ichan).spiketimes = ceil(chans(ichan).spiketimes./binsize)+spkstart;
                chans(ichan).thetaphase = NaN(numel(chans(ichan).spiketimes),nbCh+2);
            end
            
            LFPshifted = [];
            corr_LFPmua = [];
            
            for ch = 1:nbCh
                %% Load the particular channels
                fprintf([num2str(ch) ' ']);
%                 Timestamps.fulllength = length(pepNEV.ns.Data.data(chan_list(ch),:));
%                 Timestamps.samplingRateOrig = exptInfo.SamplingRateInKHZ * 1000;
%                 
                % decimate(double(pepNEV.ns.Data.data(exptInfo.Chn(ichn).ChnLoc,:)), Timestamps.samplingRateOrig/exptInfo.downsample);
                Chntemp = decimate(double(pepNEV.ns.Data.data(chan_list(ch),:)), Timestamps.samplingRateOrig/exptInfo.downsample);
%                 Timestamps.ChnsampInt = 1./exptInfo.downsample;
%                 Timestamps.Chnorig = 0:(1./Timestamps.samplingRateOrig):(Timestamps.fulllength./Timestamps.samplingRateOrig);
%                 Timestamps.Chn = 0:Timestamps.ChnsampInt:...
%                     (Timestamps.ChnsampInt*(Timestamps.fulllength*(exptInfo.downsample/Timestamps.samplingRateOrig)));
                if numel(Timestamps.Chn) > numel(Chntemp)
                    Timestamps.Chn = Timestamps.Chn(1:numel(Chntemp));
                end
                
                theta = LFPsignals(Chntemp, 1000, thetafrange(1), thetafrange(2));
                %         LFPphase(ch,:) = uint8(round(mod(360+180/pi*theta.phase,360)/2));
                
%                 procSpec.t = Timestamps.Chn;
%                 procSpec.t = procSpec.t - screenTimes(1)./(1000*nevSamplingRateInKHZ);
%                 Timestamps.Chn = Timestamps.Chn - screenTimes(1)./(exptInfo.downsample*nevSamplingRateInKHZ);
                try
                    theta = theta.hill';%interp1(Timestamps.Chn', theta.hill', es.sampleTimes);
%                     thetaphase = angle(theta);
                    thetaphs0 = round(mod(360 + 180/pi*angle(theta),360));
                    thetapeak_idx0 = find([0;abs(diff(mod(thetaphs0,360)))>180]);
                    shorttheta_lim = 90;%ms
                    shorttheta = NaN;
                    while ~isempty(shorttheta)
                        thetaperiod0 = diff(thetapeak_idx0);
                        shorttheta = find(thetaperiod0<shorttheta_lim);
                        thetapeak_idx0(shorttheta + 1) = [];
                    end
                    thetaperiod0 = diff(thetapeak_idx0);
                    thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
                    
                    thetaphase = zeros(size(thetaphs0));
                    for tt = 1:numel(thetapeak_idx0)-1
                        thetaphase(thetapeak_idx0(tt):thetapeak_idx0(tt+1)) = 360*(0:(thetapeak_idx0(tt+1)-thetapeak_idx0(tt)))/thetaperiod0(tt);
                    end
                    thetaphase = round(thetaphase(:));
                    
                    for ichan = 1:length(chans)
                        chans(ichan).thetaphase(:,ch) = thetaphase(chans(ichan).spiketimes);
                    end
                catch
                    keyboard
                end
                theta_unwraped = unwrap(thetaphase/360*2*pi);%unwrap(thetaphase);
                theta_resample = interp1(Timestamps.Chn', theta_unwraped', es.sampleTimes);
                LFPthetaPhase(ch,:) = uint8(round(mod(360/(2*pi)*theta_resample,360)/2));%uint8(round(mod(360-180/pi*angle(theta),360)/2));
                
                LFPthetaPower(ch,:) = interp1(Timestamps.Chn', abs(theta)', es.sampleTimes);
                LFPtheta(ch,:) = interp1(Timestamps.Chn', real(theta)', es.sampleTimes);
                
                if isempty(LFPshifted)
                    LFPshifted = zeros(numel(real(theta)),nbCh);
                end
                LFPshifted(:,ch) = real(theta);
            end
            
            
            if strcmp(Probe,'CA1')
                spkratemax = 10;
                spkrate = sum(es.spikeTrain,1)/size(es.spikeTrain,1)*(1/nanmean(diff(es.sampleTimes)));
                [N,~] = histcounts(cell2mat(es.chanIDs(muacellidx & spkrate < spkratemax)),0:8);
                [~,pyrtet] = sort(N,'descend');
                ch_pyrlayer = [];
                for chtet = 1:4
                    ch_pyrlayer = [ch_pyrlayer (pyrtet(1:4)-1)*4 + chtet];
                end
                muapyr = sum(es.spikeTrain(:,muacellidx & spkrate < spkratemax),2);
                mua = interp1(es.sampleTimes, muapyr, Timestamps.Chn');
                mua = mua - nanmean(mua);
                maxlag = 80;
                %instead of realiging according to mua, we realign on the
                %phase of the theta from the tetrode with a large
                %number of units that is the closest to the peak of theta
                %firing in the mua
                if isempty(chref)
                    for ch = 1:nbCh
                        corr_LFPmua2(ch,:) = xcorr(LFPshifted(~isnan(mua),ch) - nanmean(LFPshifted(~isnan(mua),ch)),mua(~isnan(mua)),maxlag,'coeff');
                        tshiftmaxLFPmua(ch) = find((corr_LFPmua2(ch,:)) == max((corr_LFPmua2(ch,:)))) - (maxlag+1);
                        tshiftminLFPmua(ch) = find((corr_LFPmua2(ch,:)) == min((corr_LFPmua2(ch,:)))) - (maxlag+1);
                    end
                    if abs(tshiftmaxLFPmua(ch_pyrlayer(1))) <= abs(tshiftminLFPmua(ch_pyrlayer(1)))
                        [~,ichref] = min(abs(tshiftmaxLFPmua(ch_pyrlayer)));
                    else
                        [~,ichref] = min(abs(tshiftminLFPmua(ch_pyrlayer)));
                    end
                    chref = ch_pyrlayer(ichref);
                end
                disp(['Reference channel: ' num2str(chref)]);
                
                for ch = 1:nbCh
                    corr_LFPmua(ch,:) = xcorr(LFPshifted(:,ch) - nanmean(LFPshifted(:,ch)),LFPshifted(:,chref) - nanmean(LFPshifted(:,chref)),maxlag,'coeff');%xcorr(LFPshifted(~isnan(mua),ch) - nanmean(LFPshifted(~isnan(mua),ch)),mua(~isnan(mua)),maxlag,'coeff');
                    tshiftmax = find((corr_LFPmua(ch,:)) == max((corr_LFPmua(ch,:)))) - (maxlag+1);
                    tshiftmin = find((corr_LFPmua(ch,:)) == min((corr_LFPmua(ch,:)))) - (maxlag+1);
                    if abs(tshiftmax) < abs(tshiftmin)
                        tshift = tshiftmax;
                        signX = +1;%-1;
                    else
                        tshift = tshiftmin;
                        signX = -1;%+1;
                    end
                    %                 tshift = find((corr) == min((corr))) - (maxlag+1);
                    LFPshifted(:,ch) = signX * circshift(LFPshifted(:,ch),-tshift);%circshift(real(theta),-tshift);
                    tshiftLFP(ch) = tshift;
                end
                corr_LFPmua2All = xcorr(nanmean(LFPshifted(~isnan(mua),ch_pyrlayer),2) - nanmean(nanmean(LFPshifted(~isnan(mua),ch_pyrlayer),2),1),mua(~isnan(mua)),maxlag,'coeff');
                tshiftmax = find((corr_LFPmua2All) == max((corr_LFPmua2All))) - (maxlag+1);
                tshiftmin = find((corr_LFPmua2All) == min((corr_LFPmua2All))) - (maxlag+1);
                if abs(tshiftmax) < abs(tshiftmin)
                    tshift = tshiftmax;
                    signX = -1;
                else
                    tshift = tshiftmin;
                    signX = +1;
                end
                LFPshifted = signX*LFPshifted;
            else
                ch_pyrlayer = 1:nbChV1;
            end
            
%             maxlag = 80;
%             corr2 = [];
%             if isempty(bestCh)
%                 for ch = 1:nbCh
%                     corr2(ch,:) = xcorr(mean(LFPshifted(~isnan(mua),ch),2),mua(~isnan(mua)),maxlag,'coeff');
%                 end
%                 [~, bestCh] = max(max(abs(corr2(:,floor(maxlag/2):(maxlag + floor(maxlag/2)))),[],2));
%             end
            
            %theta from mua
            theta = LFPsignals(sum(es.spikeTrain(:,muacellidx),2)', 60, thetafrange(1), thetafrange(2));  
            theta = theta.hill';
            thetaphs0 = round(mod(360 + 180/pi*angle(theta),360));
            thetapeak_idx0 = find([0;abs(diff(mod(thetaphs0,360)))>180]);
            thetaperiod0 = diff(thetapeak_idx0);
            thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
            
            thetaphase = zeros(size(thetaphs0));
            for tt = 1:numel(thetapeak_idx0)-1
                thetaphase(thetapeak_idx0(tt):thetapeak_idx0(tt+1)) = 360*(0:(thetapeak_idx0(tt+1)-thetapeak_idx0(tt)))/thetaperiod0(tt);
            end
            thetaphase = round(thetaphase(:));
            
            LFPthetaPhase(nbCh+1,:) = uint8(round(mod(thetaphase,360)/2));%uint8(round(mod(360-180/pi*angle(theta),360)/2));
            LFPthetaPower(nbCh+1,:) = abs(theta);
            LFPtheta(nbCh+1,:) = real(theta);
            
            theta = LFPsignals(mean(LFPshifted(:,ch_pyrlayer),2)', 1000, thetafrange(1), thetafrange(2));
            theta = theta.hill';
            thetaphs0 = round(mod(360 + 180/pi*angle(theta),360));
            thetapeak_idx0 = find([0;abs(diff(mod(thetaphs0,360)))>180]);
            shorttheta_lim = 90;%ms
            shorttheta = NaN;
            while ~isempty(shorttheta)
                thetaperiod0 = diff(thetapeak_idx0);
                shorttheta = find(thetaperiod0<shorttheta_lim);
                thetapeak_idx0(shorttheta + 1) = [];
            end
            thetaperiod0 = diff(thetapeak_idx0);
            thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
            
            thetaphase = zeros(size(thetaphs0));
            for tt = 1:numel(thetapeak_idx0)-1
                thetaphase(thetapeak_idx0(tt):thetapeak_idx0(tt+1)) = 360*(0:(thetapeak_idx0(tt+1)-thetapeak_idx0(tt)))/thetaperiod0(tt);
            end
            thetaphase = round(thetaphase(:));
            for ichan = 1:length(chans)
                chans(ichan).thetaphase(:,nbCh+2) = thetaphase(chans(ichan).spiketimes);
            end
            
            theta_unwraped = unwrap(2*pi*thetaphase/360);%unwrap(angle(theta));
            theta_resample = interp1(Timestamps.Chn', theta_unwraped', es.sampleTimes);
            LFPthetaPhase(nbCh+2,:) = uint8(round(mod(360/(2*pi)*theta_resample,360)/2));%uint8(round(mod(360 - 180/pi*theta_resample,360)/2));
            LFPthetaPower(nbCh+2,:) = interp1(Timestamps.Chn', abs(theta)', es.sampleTimes);
            LFPtheta(nbCh+2,:) = interp1(Timestamps.Chn', real(theta)', es.sampleTimes);
            
            CSDthetaPhase = zeros(size(LFPthetaPhase));
            CSDthetaPower = zeros(size(LFPthetaPower));
            CSDtheta = zeros(size(LFPtheta));
            if strcmp(Probe,'V1')
                for ichan = 1:length(chans)
                    chans(ichan).CSDthetaphase = NaN(numel(chans(ichan).spiketimes),nbCh+2);
                end
                Espacing=20e-6;
                el_pos=(Espacing:Espacing:nbChCA1*Espacing)';
                
                CSDthetafull = zeros(size(LFPshifted));
                Ntimes = size(CSDthetafull,1);
                Nwin = 100000;
                Nseg = floor(Ntimes/Nwin)+1;
                for tk = 1:Nseg
                    CSDthetafull(((tk-1)*Nwin + 1):min(tk*Nwin,Ntimes),:) = (LFP2iCSD(LFPshifted(((tk-1)*Nwin + 1):min(tk*Nwin,Ntimes),:)', el_pos))';
                end
                
                for ch = 1:nbChCA1
                    theta = LFPsignals(mean(CSDthetafull(:,ch),2)', 1000, thetafrange(1), thetafrange(2));
                    theta = theta.hill';
                    thetaphs0 = round(mod(360 + 180/pi*angle(theta),360));
                    thetapeak_idx0 = find([0;abs(diff(mod(thetaphs0,360)))>180]);
                    shorttheta_lim = 90;%ms
                    shorttheta = NaN;
                    while ~isempty(shorttheta)
                        thetaperiod0 = diff(thetapeak_idx0);
                        shorttheta = find(thetaperiod0<shorttheta_lim);
                        thetapeak_idx0(shorttheta + 1) = [];
                    end
                    thetaperiod0 = diff(thetapeak_idx0);
                    thetaperiod0 = [thetaperiod0;thetaperiod0(end)];
                    
                    thetaphase = zeros(size(thetaphs0));
                    for tt = 1:numel(thetapeak_idx0)-1
                        thetaphase(thetapeak_idx0(tt):thetapeak_idx0(tt+1)) = 360*(0:(thetapeak_idx0(tt+1)-thetapeak_idx0(tt)))/thetaperiod0(tt);
                    end
                    thetaphase = round(thetaphase(:));
                    for ichan = 1:length(chans)
                        chans(ichan).CSDthetaphase(:,ch) = thetaphase(chans(ichan).spiketimes);
                    end
                    
                    theta_unwraped = unwrap(2*pi*thetaphase/360);%unwrap(angle(theta));
                    theta_resample = interp1(Timestamps.Chn', theta_unwraped', es.sampleTimes);
                    CSDthetaPhase(ch,:) = uint8(round(mod(360/(2*pi)*theta_resample,360)/2));%uint8(round(mod(360 - 180/pi*theta_resample,360)/2));
                    CSDthetaPower(ch,:) = interp1(Timestamps.Chn', abs(theta)', es.sampleTimes);
                    CSDtheta(ch,:) = interp1(Timestamps.Chn', real(theta)', es.sampleTimes);
                end
            end
            
            
            
            fprintf('\r');
            timestamps = Timestamps.Chn;
            
            global DIRS
            fname = [animal '_' num2str(iseries(s)) '_' num2str(iexp)];
            dDIRname = [DIRS.multichanspikes filesep animal filesep num2str(iseries(s))];
            
            SpkthetaPhase = cell(1,length(chans));
            for ichan = 1:length(chans)
                SpkthetaPhase{ichan} = chans(ichan).thetaphase;
            end
            if strcmp(Probe,'V1')
                SpkCSDthetaPhase = cell(1,length(chans));
                for ichan = 1:length(chans)
                    SpkCSDthetaPhase{ichan} = chans(ichan).CSDthetaphase;
                end
            end
            for ichan = 1:length(chans)
                if size(chans(ichan).thetaphase,1) ~= numel(chans(ichan).spiketimes)
                    keyboard
                end
            end
%             dDIRnameSave = [DIRS.data2 filesep animal filesep num2str(iseries(s))];
%             if ~isdir(dDIRnameSave)
%                 mkdir(dDIRnameSave)
%             end
            if strcmp(Probe,'CA1')
                save([dDIRname filesep fname '_LFP3.mat'], 'LFPthetaPhase','LFPthetaPower','LFPtheta','corr_LFPmua','SpkthetaPhase','chref','corr_LFPmua2All','-v7.3');%,'LFP', 'LFPphase', 'theta',  'timestamps','-v7.3');
            else
                save([dDIRname filesep fname '_' Probe '_LFP3.mat'], 'LFPthetaPhase','LFPthetaPower','LFPtheta', 'CSDthetaPhase','CSDthetaPower','CSDtheta','corr_LFPmua','SpkthetaPhase','SpkCSDthetaPhase','-v7.3');%,'LFP', 'LFPphase', 'theta',  'timestamps','-v7.3');
            end
        end
%         catch
%             warning(['problem in ' fname]);
%         end
    end
end

clear global pepNEV

end

% to batch across sessions
% expt = getExperimentList;
% for ianimal = 1:numel(expt)
% for iseries = 1:numel(expt(ianimal).series)
%     if expt(ianimal).CA1{iseries}
%         preprocessLFP_theta2(expt(ianimal).animal, expt(ianimal).series{iseries}, [], 'CA1');
%     end
%     if expt(ianimal).V1{iseries}
%         preprocessLFP_theta2(expt(ianimal).animal, expt(ianimal).series{iseries}, [], 'V1');
%     end
% end
% end

