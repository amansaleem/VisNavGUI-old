function preprocessLFP_theta(animal, iseries, iexplist, Probe)
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
    bestCh = [];
    for i = 1:numel(iexplist)
        try
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
            nbCh = nbChCA1;%numel(chan_list);
            
            if strcmp(Probe,'V1')
                chan_list = chan_list(nbChCA1+1:nbChCA1+32);
            end
            
            LFPthetaPower = [];
            LFPthetaPhase = [];
            LFPtheta = [];
            thetafrange  = [6 9];%[4 9];%
            fprintf('Dowsampling Channel #: ');
            
            es = getVRspikes(animal,iseries(s),iexp,1,100,1,0,0:7,{Probe});
            
            nbins = sum(Timestamps.Chn>0);
            spkstart = find(Timestamps.Chn>0,1,'first');
            binsize = 1/exptInfo.downsample;
            for ichan = 1:length(chans)
                st = ceil(chans(ichan).spiketimes./binsize);
                numExcessSpikes = sum(st>nbins);
                st(st>nbins) = [];
                st(st==0) = [];
                if numExcessSpikes>0
                    chans(ichan).spiketimes(end-numExcessSpikes+1:end) = [];
                end
                chans(ichan).spiketimes = ceil(chans(ichan).spiketimes./binsize)+spkstart;
                chans(ichan).thetaphase = ThetaPhs_all(chans(ichan).spiketimes,:);
            end
            
            nevSamplingRateInKHZ = 30;
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
                
                theta = LFPsignals(Chntemp, 1000, thetafrange(1), thetafrange(2));
                %         LFPphase(ch,:) = uint8(round(mod(360+180/pi*theta.phase,360)/2));
                
                procSpec.t = Timestamps.Chn;
                procSpec.t = procSpec.t - screenTimes(1)./(1000*nevSamplingRateInKHZ);
                Timestamps.Chn = Timestamps.Chn - screenTimes(1)./(exptInfo.downsample*nevSamplingRateInKHZ);
                try
                    theta = theta.hill';%interp1(Timestamps.Chn', theta.hill', es.sampleTimes);
                    mua = interp1(es.sampleTimes, es.mua, Timestamps.Chn');
                catch
                    keyboard
                end
                theta_unwraped = unwrap(angle(theta));
                theta_resample = interp1(Timestamps.Chn', theta_unwraped', es.sampleTimes);
                LFPthetaPhase(ch,:) = uint8(round(mod(360 - 180/pi*theta_resample,360)/2));%uint8(round(mod(360-180/pi*angle(theta),360)/2));
                
                LFPthetaPower(ch,:) = interp1(Timestamps.Chn', abs(theta)', es.sampleTimes);
                LFPtheta(ch,:) = interp1(Timestamps.Chn', real(theta)', es.sampleTimes);
                
                if isempty(LFPshifted)
                    LFPshifted = zeros(numel(real(theta)),nbCh);
                end
                LFPshifted(:,ch) = real(theta);
                if strcmp(Probe,'CA1')
                    maxlag = 80;
                    corr_LFPmua(ch,:) = xcorr(LFPshifted(~isnan(mua),ch),mua(~isnan(mua)),maxlag,'coeff');%xcorr(real(theta),es.mua,maxlag,'coeff');
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
            theta = LFPsignals(es.mua', 60, thetafrange(1), thetafrange(2));  
            theta = theta.hill';
            LFPthetaPhase(nbCh+1,:) = uint8(round(mod(360-180/pi*angle(theta),360)/2));
            LFPthetaPower(nbCh+1,:) = abs(theta);
            LFPtheta(nbCh+1,:) = real(theta);
            
            theta = LFPsignals(mean(LFPshifted(:,1:nbChCA1),2)', 1000, thetafrange(1), thetafrange(2));
            theta = theta.hill';
            theta_unwraped = unwrap(angle(theta));
            theta_resample = interp1(Timestamps.Chn', theta_unwraped', es.sampleTimes);
            LFPthetaPhase(nbCh+2,:) = uint8(round(mod(360 - 180/pi*theta_resample,360)/2));
            LFPthetaPower(nbCh+2,:) = interp1(Timestamps.Chn', abs(theta)', es.sampleTimes);
            LFPtheta(nbCh+2,:) = interp1(Timestamps.Chn', real(theta)', es.sampleTimes);
            
            if strcmp(Probe,'V1')
                Espacing=20e-6;
                el_pos=(Espacing:Espacing:nbChCA1*Espacing)';
                
                CSDtheta = zeros(size(LFPshifted));
                Ntimes = size(CSDtheta,1);
                Nwin = 100000;
                Nseg = floor(Ntimes/Nwin)+1;
                for tk = 1:Nseg
                    CSDtheta(((tk-1)*Nwin + 1):min(tk*Nwin,Ntimes),:) = (LFP2iCSD(LFPshifted(((tk-1)*Nwin + 1):min(tk*Nwin,Ntimes),:)', el_pos))';
                end
                theta = LFPsignals(mean(CSDtheta(:,1:nbChCA1),2)', 1000, thetafrange(1), thetafrange(2));
                theta = theta.hill';
                theta_unwraped = unwrap(angle(theta));
                theta_resample = interp1(Timestamps.Chn', theta_unwraped', es.sampleTimes);
                CSDthetaPhase = uint8(round(mod(360 - 180/pi*theta_resample,360)/2));
                CSDthetaPower = interp1(Timestamps.Chn', abs(theta)', es.sampleTimes);
                CSDtheta = interp1(Timestamps.Chn', real(theta)', es.sampleTimes);
            end
            
            
            
            fprintf('\r');
            timestamps = Timestamps.Chn;
            
            global DIRS
            fname = [animal '_' num2str(iseries(s)) '_' num2str(iexp)];
            dDIRname = [DIRS.multichanspikes filesep animal filesep num2str(iseries(s))];
%             dDIRnameSave = [DIRS.data2 filesep animal filesep num2str(iseries(s))];
%             if ~isdir(dDIRnameSave)
%                 mkdir(dDIRnameSave)
%             end
            if strcmp(Probe,'CA1')
                save([dDIRname filesep fname '_LFP.mat'], 'LFPthetaPhase','LFPthetaPower','LFPtheta','corr_LFPmua','-v7.3');%,'LFP', 'LFPphase', 'theta',  'timestamps','-v7.3');
            else
                save([dDIRname filesep fname '_' Probe '_LFP.mat'], 'LFPthetaPhase','LFPthetaPower','LFPtheta', 'CSDthetaPhase','CSDthetaPower','CSDtheta','corr_LFPmua','-v7.3');%,'LFP', 'LFPphase', 'theta',  'timestamps','-v7.3');
            end
        end
        catch
            warning(['problem in ' fname]);
        end
    end
end

clear global pepNEV

end

% to batch across sessions
% expt = getExperimentList;
% for ianimal = 1:numel(expt)
% for iseries = 1:numel(expt(ianimal).series)
%     if expt(ianimal).CA1{iseries}
%         preprocessLFP_theta(expt(ianimal).animal, expt(ianimal).series{iseries}, [], 'CA1');
%     end
%     if expt(ianimal).V1{iseries}
%         preprocessLFP_theta(expt(ianimal).animal, expt(ianimal).series{iseries}, [], 'V1');
%     end
% end
% end

