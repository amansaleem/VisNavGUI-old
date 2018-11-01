function cellprop = BatchCellInfo
SetDirs;
cellprop = [];


strlistvarname = {'2p data','electrophys data'};
[varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'dataset', 'SelectionMode', 'single', 'InitialValue', 1);
if ok && varnamesel == 1
    batch2p = true;
    expt = getExperimentList2p;
    datadir = 'D:\DATA\batch\2p';%DIRS.data2p;
    
    cellprop.Tsmthwin = 250;%250;%250;%250;%150;%300;%40;%120;%50
    cellprop.Xsmthwin = 2;%1;%
    cellprop.SpeedThreshold = 1;
    cellprop.nthetaphsbins = 0;%1;%
    cellprop.cellstr = 'goodonly';
    
    suffix = ['Twin' num2str(cellprop.Tsmthwin) '_Xwin' num2str(cellprop.Xsmthwin) '_spdth' num2str(cellprop.SpeedThreshold) '_' num2str(cellprop.nthetaphsbins) 'thetabins_' cellprop.cellstr];
    disp(suffix)
    suffix_cellprop = [suffix '_cellProperties'];
    suffix_map1d = [suffix '_maps1d'];
elseif ok
    batch2p = false;
    expt = getExperimentList;
    datadir = 'D:\DATA\batch';%datadir = DIRS.multichanspikes;
    
    cellprop.Tsmthwin = 15;%15;%250;%150;%300;%
    cellprop.Xsmthwin = 4;%4;%1;%
    cellprop.SpeedThreshold = 5;
    cellprop.nthetaphsbins = 0;%1;%
    cellprop.cellstr = 'goodonly';
    suffix = ['Twin' num2str(cellprop.Tsmthwin) '_Xwin' num2str(cellprop.Xsmthwin) '_spdth' num2str(cellprop.SpeedThreshold) '_' num2str(cellprop.nthetaphsbins) 'thetabins_' cellprop.cellstr];
    disp(suffix)
    suffix_cellprop = [suffix '_cellProperties'];
    suffix_map1d = [suffix '_maps1d'];
    suffix_map2dtheta = [suffix '_maps2d_phs'];
    suffix_lfpspkcoherence = [cellprop.cellstr '_LFPSpikeCoherence'];
end

nanimal = numel(expt);

nbXpos = 10; nbYpos = 9; nVStimepts = 61;

for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        dDIRname = [datadir  filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed'];
        if exist([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' suffix_cellprop '.mat'],'file')
            S = load([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' suffix_cellprop '.mat']);
            Smaps1d = load([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' suffix_map1d '.mat']);
            Smaps2dtheta = load([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' suffix_map2dtheta '.mat']);
            Slfpspkcoherence = load([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' suffix_lfpspkcoherence '.mat']);
            if ~isfield(cellprop,'animal')
                cellprop.animal = ianimal*ones(1,S.CellInfo.NumCells);
                cellprop.iseries = iseries*ones(1,S.CellInfo.NumCells);
                cellprop.CellList = S.CellInfo.CellList;
                cellprop.Shank = S.CellInfo.Shank;
                cellprop.MUAcluster = S.CellInfo.MUAcluster;
                cellprop.Goodcluster = S.CellInfo.Goodcluster;
                cellprop.Finterneuron = S.CellInfo.Finterneuron;
                cellprop.Unsortedcluster = S.CellInfo.Unsortedcluster;
                if isfield(S.CellInfo,'rateautocorr')
                    cellprop.Spatialmodulation = S.CellInfo.Spatialmodulation;
                    cellprop.burstindex = S.CellInfo.burstindex;
                    cellprop.minFR = S.CellInfo.minFR;
                    cellprop.rateautocorr = S.CellInfo.rateautocorr;
                    cellprop.min2maxSpkwf = S.CellInfo.min2maxSpkwf;
                    cellprop.spkautocorr = S.CellInfo.spkautocorr;
                end
                cellprop.Probe = S.CellInfo.Probe;
                cellprop.bestchan = S.CellInfo.bestchan;
                if strcmp(expt(ianimal).area{iseries},'V1medial')
                    cellprop.Cellpos2p = 1*ones(1,S.CellInfo.NumCells);
                elseif strcmp(expt(ianimal).area{iseries},'V1lateral')
                    cellprop.Cellpos2p = 2*ones(1,S.CellInfo.NumCells);
                elseif strcmp(expt(ianimal).area{iseries},'PPC')
                    cellprop.Cellpos2p = 3*ones(1,S.CellInfo.NumCells);
                elseif strcmp(expt(ianimal).area{iseries},'AL')
                    cellprop.Cellpos2p = 4*ones(1,S.CellInfo.NumCells);
                else
                    cellprop.Cellpos2p = zeros(1,S.CellInfo.NumCells);
                end
                if ~isempty(S.CellInfo.RFXpos)
                    cellprop.RFXpos = S.CellInfo.RFXpos;
                    cellprop.Xpos = S.CellInfo.Xpos;
                    cellprop.XposZmax = S.CellInfo.XposZmax;
                    cellprop.globalXposrep = S.CellInfo.globalXposrep;
                    [~,xposmax] = max(mean(S.CellInfo.RFXpos(S.CellInfo.XposZmax>2,:),1));
                    cellprop.XposPop = xposmax*ones(1,S.CellInfo.NumCells);

                    cellprop.RFYpos = S.CellInfo.RFYpos;
                    cellprop.Ypos = S.CellInfo.Ypos;
                    cellprop.YposZmax = S.CellInfo.YposZmax;
                    cellprop.globalYposrep = S.CellInfo.globalYposrep;
                    [~,yposmax] = max(mean(S.CellInfo.RFYpos(S.CellInfo.YposZmax>2,:),1));
                    cellprop.YposPop = yposmax*ones(1,S.CellInfo.NumCells);
                    
                    cellprop.VStime = (-20:40)*0.05;
                else
                    cellprop.RFXpos = NaN(S.CellInfo.NumCells,nbXpos);
                    cellprop.Xpos = NaN(1,S.CellInfo.NumCells);
                    cellprop.XposZmax = NaN(1,S.CellInfo.NumCells);
                    cellprop.globalXposrep = NaN(S.CellInfo.NumCells,nVStimepts);
                    cellprop.XposPop = NaN(1,S.CellInfo.NumCells);

                    cellprop.RFYpos = NaN(S.CellInfo.NumCells,nbYpos);
                    cellprop.Ypos = NaN(1,S.CellInfo.NumCells);
                    cellprop.YposZmax = NaN(1,S.CellInfo.NumCells);
                    cellprop.globalYposrep = NaN(S.CellInfo.NumCells,nVStimepts);
                    cellprop.YposPop = NaN(1,S.CellInfo.NumCells);
                    
                    cellprop.VStime = (-20:40)*0.05;
                end
                
                for iprobe = 1:2
                    if ~isempty(S.CellInfo.LFP2Spike_phscorrMUA{iprobe})
                        cellprop.phsfieldMUA{iprobe} = S.CellInfo.phsfieldMUA{iprobe}';
                        cellprop.phsfieldMUASE{iprobe} = S.CellInfo.phsfieldMUASE{iprobe}';
                        cellprop.phsfieldZMUA{iprobe} = S.CellInfo.phsfieldZMUA{iprobe};
                        cellprop.phsfieldPosMUA{iprobe} = S.CellInfo.phsfieldPosMUA{iprobe};
                        cellprop.LFP2Spike_phscorrMUA{iprobe} = S.CellInfo.LFP2Spike_phscorrMUA{iprobe};
                        
                        cellprop.phsfieldMUAnorm{iprobe} = S.CellInfo.phsfieldMUAnorm{iprobe};
                        cellprop.phsfieldPosMUAnorm{iprobe} = S.CellInfo.phsfieldPosMUAnorm{iprobe};
                        cellprop.LFP2Spike_phscorrMUAnorm{iprobe} = S.CellInfo.LFP2Spike_phscorrMUAnorm{iprobe};
                        
                        cellprop.phsfieldMean{iprobe} = S.CellInfo.phsfieldMean{iprobe};
                        cellprop.phsfieldPosMean{iprobe} = S.CellInfo.phsfieldPosMean{iprobe};
                        cellprop.LFP2Spike_phscorrMean{iprobe} = S.CellInfo.LFP2Spike_phscorrMean{iprobe};
                    else
                        cellprop.phsfieldMUA{iprobe} = NaN(size(S.CellInfo.phsfieldMUA{1}'));
                        cellprop.phsfieldMUASE{iprobe} = NaN(size(S.CellInfo.phsfieldMUASE{1}'));
                        cellprop.phsfieldZMUA{iprobe} = NaN(size(S.CellInfo.phsfieldZMUA{1}));
                        cellprop.phsfieldPosMUA{iprobe} = NaN(size(S.CellInfo.phsfieldPosMUA{1}));
                        cellprop.LFP2Spike_phscorrMUA{iprobe} = NaN(size(S.CellInfo.LFP2Spike_phscorrMUA{1}));
                        
                        cellprop.phsfieldMUAnorm{iprobe} = NaN(size(S.CellInfo.phsfieldMUAnorm{1}));
                        cellprop.phsfieldPosMUAnorm{iprobe} = NaN(size(S.CellInfo.phsfieldPosMUAnorm{1}));
                        cellprop.LFP2Spike_phscorrMUAnorm{iprobe} = NaN(size(S.CellInfo.LFP2Spike_phscorrMUAnorm{1}));
                        
                        cellprop.phsfieldMean{iprobe} = NaN(size(S.CellInfo.phsfieldMean{1}));
                        cellprop.phsfieldPosMean{iprobe} = NaN(size(S.CellInfo.phsfieldPosMean{1}));
                        cellprop.LFP2Spike_phscorrMean{iprobe} = NaN(size(S.CellInfo.LFP2Spike_phscorrMean{1}));
                    end
                end
                
                allcontidx = size(S.CellInfo.field,1);
                for g = [2 1 3]
                    cellprop.field{g} = S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldSE{g} = S.CellInfo.fieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldZ{g} = S.CellInfo.fieldZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldPos{g} = S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldPosSE{g} = S.CellInfo.fieldPosSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldCOM{g} = S.CellInfo.fieldCOM{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldXcorr{g} = S.CellInfo.fieldXcorr{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.rate{g} = S.CellInfo.rate{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldMax{g} = S.CellInfo.fieldMax{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldMin{g} = S.CellInfo.fieldMin{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.SSI{g} = S.CellInfo.SSI{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.SpatialInfo{g} = S.CellInfo.SpatialInfo{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.SpatialInfoPerSpike{g} = S.CellInfo.SpatialInfoPerSpike{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
%                     cellprop.thetaPhaseMax{g} = S.CellInfo.PhaseMax{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXPhstheta{g} = S.CellInfo.field2dXPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXPhsthetaSE{g} = S.CellInfo.field2dXPhsthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXthetapos{g} = S.CellInfo.field2dXthetapos{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXthetaposSE{g} = S.CellInfo.field2dXthetaposSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXthetaposNorm{g} = S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXthetaposSENorm{g} = S.CellInfo.field2dXthetaposSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    cellprop.field2dPhstheta{g} = S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dPhsthetaSE{g} = S.CellInfo.field2dPhsthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    cellprop.field2dXcorrtheta{g} = S.CellInfo.field2dXcorrtheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetaSE{g} = S.CellInfo.field2dXcorrthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetamax{g} = S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetamaxSE{g} = S.CellInfo.field2dXcorrthetamaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    cellprop.field2dXcorrthetamax_set1{g} = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                    cellprop.field2dXcorrthetamax_set2{g} = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                    thetafieldset1 = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                    thetafieldset2 = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                    cellprop.field2dXcorrthetamax_set1{g} = thetafieldset1;
                    cellprop.field2dXcorrthetamax_set2{g} = thetafieldset2;
                    thetareliabilityCorr = sum((thetafieldset1-repmat(nanmean(thetafieldset1,2),[1 size(thetafieldset1,2)])).*...
                        (thetafieldset2-repmat(nanmean(thetafieldset2,2),[1 size(thetafieldset2,2)])),2)./...
                        sqrt(sum((thetafieldset1-repmat(nanmean(thetafieldset1,2),[1 size(thetafieldset1,2)])).^2,2).*...
                        sum((thetafieldset2-repmat(nanmean(thetafieldset2,2),[1 size(thetafieldset2,2)])).^2,2));
                    cellprop.thetareliabilityCorr{g} = thetareliabilityCorr;
                    
                    cellprop.phsfieldPos{g} = S.CellInfo.phsfieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldPosSE{g} = S.CellInfo.phsfieldPosSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldCOM{g} = S.CellInfo.phsfieldCOM{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldCOMSE{g} = S.CellInfo.phsfieldCOMSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfield{g} = S.CellInfo.phsfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldSE{g} = S.CellInfo.phsfieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldZ{g} = S.CellInfo.phsfieldZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsmodulation{g} = S.CellInfo.phsmodulation{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    ZXposthetaMax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    ZXposthetaMin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    XposthetaMax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    XposthetaMin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    ZXPhsthetaMax = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    ZXPhsthetaMin = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    PhsthetaMod = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    Phsthetacross = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    phaseXmin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    phaseXmax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    phaseXahead = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    phaseXbehind = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    phaseXcross = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    nphsbins = size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},2);
                    for icell = 1:size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                         Xpostheta = S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
%                         XposthetaSE = S.CellInfo.field2dXthetaposSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                        Xpostheta = S.CellInfo.field2dXcorrthetamaxNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                        XposthetaSE = S.CellInfo.field2dXcorrthetamaxSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                        [~, imax] = max(Xpostheta);
                        ZXposthetaMax(icell) = abs((Xpostheta(imax) - mean(Xpostheta))./XposthetaSE(imax));
                        [~, imin] = min(Xpostheta);
                        ZXposthetaMin(icell) = abs((Xpostheta(imin) - mean(Xpostheta))./XposthetaSE(imin));
                        
                        Xpostheta = S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                        [XposthetaMax(icell), ~] = max(Xpostheta);
                        [XposthetaMin(icell), ~] = min(Xpostheta);
                        
                        Phstheta = S.CellInfo.phsfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                        PhsthetaSE = S.CellInfo.phsfieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                        [~, imax] = max(Phstheta);
                        ZXPhsthetaMax(icell) = abs((Phstheta(imax) - mean(Phstheta))./PhsthetaSE(imax));
                        [~, imin] = min(Phstheta);
                        ZXPhsthetaMin(icell) = abs((Phstheta(imin) - mean(Phstheta))./PhsthetaSE(imin));
%                         Phsthetacross(icell) = (imax-1)*360/nphsbins;
                        
                        Phstheta = S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:)./mean(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:));
                        Phi = (0:(nphsbins-1))/nphsbins*2*pi;
                        Phi_num = sum(Phstheta.*sin(Phi));
                        Phi_den = sum(Phstheta.*cos(Phi));
                        Phsthetacross(icell) = mod(90 - 360*atan2(Phi_den,Phi_num)/(2*pi),360);
                        [PhsthetaMod(icell), ~] = max(Phstheta);
                        
                        Xpostheta = squeeze(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:));
                        imax = getCircularAverage(Xpostheta',0,0.01,0.05);%[~, imax] = max(Xpostheta);
                        phaseXmax(icell) = imax;
                        imin = getCircularAverage(-Xpostheta',0,0.01,0.05);%[~, imin] = min(Xpostheta);
                        phaseXmin(icell) = imin;
                        
                        Xpostheta_ahead = Xpostheta - mean(Xpostheta);
                        Xpostheta_ahead(Xpostheta_ahead>0) = 0;
                        phaseXahead(icell) = getCircularAverage(abs(Xpostheta_ahead)',0,1);
                        Xpostheta_behind = Xpostheta - mean(Xpostheta);
                        Xpostheta_behind(Xpostheta_behind<0) = 0;
                        phaseXbehind(icell) = getCircularAverage(abs(Xpostheta_behind)',0,1);
                        
                        Phi = (0:(nphsbins-1))/nphsbins*2*pi;
                        Phi_num = sum(Xpostheta.*sin(Phi));
                        Phi_den = sum(Xpostheta.*cos(Phi));
                        phaseXcross(icell) = 180 - 360*atan2(Phi_den,Phi_num)/(2*pi);
                    end
                    cellprop.ZXmaxthetapos{g} = ZXposthetaMax;
                    cellprop.ZXminthetapos{g} = ZXposthetaMin;
                    cellprop.Xmaxthetapos{g} = XposthetaMax;
                    cellprop.Xminthetapos{g} = XposthetaMin;
                    cellprop.ZPhsmaxtheta{g} = ZXPhsthetaMax;
                    cellprop.ZPhsmintheta{g} = ZXPhsthetaMin;
                    cellprop.thetaPhaseMax{g} = Phsthetacross;
                    cellprop.PhsthetaMod{g} = PhsthetaMod;
                    cellprop.PhsXmaxthetapos{g} = phaseXmax;
                    cellprop.PhsXminthetapos{g} = phaseXmin;
                    cellprop.PhsXaheadthetapos{g} = phaseXahead;
                    cellprop.PhsXbehindthetapos{g} = phaseXbehind;
                    cellprop.PhsXcrossthetapos{g} = phaseXcross;
                    
%                     ZXposthetaMax = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                     ZXposthetaMin = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                     phaseXmin = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                     phaseXmax = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                     for icell = 1:size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                         Xpostheta = S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
%                         XposthetaSE = S.CellInfo.field2dXthetaposSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
%                         [~, imax] = max(Xpostheta);
%                         phaseXmax(icell) = imax;
%                         ZXposthetaMax(icell) = abs((Xpostheta(imax) - mean(Xpostheta))./XposthetaSE(imax));
%                         [~, imin] = min(Xpostheta);
%                         ZXposthetaMin(icell) = abs((Xpostheta(imin) - mean(Xpostheta))./XposthetaSE(imin));
%                         phaseXmin(icell) = imin;
%                     end
%                     cellprop.ZXmaxthetaposNorm{g} = ZXposthetaMax;
%                     cellprop.ZXminthetaposNorm{g} = ZXposthetaMin;
%                     cellprop.PhsXmaxthetaposNorm{g} = phaseXmax;
%                     cellprop.PhsXminthetaposNorm{g} = phaseXmin;
                    
                    cellprop.PhaseRho{g} = S.CellInfo.PhaseRho{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.PhasePval{g} = S.CellInfo.PhasePval{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.PhaseRayleighPval{g} = S.CellInfo.PhaseRayleighPval{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.PhaseRayleighZ{g} = S.CellInfo.PhaseRayleighZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    cellprop.field_half1{g} = S.CellInfo.field_half1{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field_half2{g} = S.CellInfo.field_half2{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    fieldset1 = squeeze(S.CellInfo.field_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                    fieldset2 = squeeze(S.CellInfo.field_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                    cellprop.field_set1{g} = fieldset1;
                    cellprop.field_set2{g} = fieldset2;
                    reliabilityCorr = sum((fieldset1-repmat(nanmean(fieldset1,2),[1 size(fieldset1,2)])).*...
                                               (fieldset2-repmat(nanmean(fieldset2,2),[1 size(fieldset2,2)])),2)./...
                                               sqrt(sum((fieldset1-repmat(nanmean(fieldset1,2),[1 size(fieldset1,2)])).^2,2).*...
                                               sum((fieldset2-repmat(nanmean(fieldset2,2),[1 size(fieldset2,2)])).^2,2));
                    cellprop.reliabilityCorr{g} = reliabilityCorr;
                    
                    
%                     maxpos = [];
%                     for icell = 1:size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                         map = repmat(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:),[1 3]);
%                         xorig = 0.5:size(map,2)-0.5;
%                         xinterp = 0.1:0.1:size(map,2);
%                         mapinterp = interp1(xorig,map,xinterp,'spline');
%                         mapinterp = mapinterp((numel(mapinterp)/3+1):(2*numel(mapinterp)/3));
%                         [~,imax] = max(mapinterp);
%                         maxpos(icell) = xinterp(imax);
%                     end
%                     cellprop.fieldPos{g} = maxpos(:)';
                    
                    maxcorr{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                    maxcorrstd{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                    if ~isempty(Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model)
                        for icell = 1:numel(Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning)
                            map = Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).meanrespModel;
                            mapbase = Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).meanrespModel;
                            
                            map = map - nanmean(map);
                            map = map./sqrt(sum(map.^2));
                            mapbase = mapbase - nanmean(mapbase);
                            mapbase = mapbase./sqrt(sum(mapbase.^2));
                            fieldXcorr_all = zeros(1,numel(mapbase));
                            xshiftlim = floor(numel(mapbase)/2);
                            ishift = 0;%floor(numel(mapbase)/2)-xshiftlim-1;
                            for xshift = -xshiftlim:xshiftlim-1
                                ishift = ishift + 1;
                                fieldXcorr_all(ishift) = map*circshift(mapbase,xshift)';
                            end
                            
%                             xorig = 0.5:numel(fieldXcorr_all)-0.5;
%                             xinterp = 0:0.1:numel(fieldXcorr_all);
%                             
%                             if sum(isnan(fieldXcorr_all)) == 0
%                                 fieldXcorr_all = interp1(xorig,fieldXcorr_all,xinterp,'spline');
%                             else
%                                 fieldXcorr_all = interp1(xorig,fieldXcorr_all,xinterp);
%                             end
                            
                            maxcorr_all = getCircularAverage(fieldXcorr_all(:),0,0.1,0.05);
%                             [~ , maxcorr_all] = max(fieldXcorr_all);
%                             maxcorr_all = xinterp(maxcorr_all) - xinterp(floor(size(fieldXcorr_all,2)/2) + 1);
                            
                            stdmaxcorr = 0;
                            kfold = size(Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel,1);
                            %SEM computed using Jacknife method
                            for i = 1:kfold
                                map_iter = Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel(i,:);%map;%
                                mapbase_iter = mapbase;%Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel(i,:);%
                                
                                map_iter = map_iter - nanmean(map_iter);
                                map_iter = map_iter./sqrt(sum(map_iter.^2));
                                mapbase_iter = mapbase_iter - nanmean(mapbase_iter);
                                mapbase_iter = mapbase_iter./sqrt(sum(mapbase_iter.^2));
                                fieldXcorr_iter = zeros(1,numel(mapbase_iter));
                                xshiftlim = floor(numel(mapbase_iter)/2);
                                ishift = 0;%floor(numel(mapbase_iter)/2)-xshiftlim-1;
                                for xshift = -xshiftlim:xshiftlim-1
                                    ishift = ishift + 1;
                                    fieldXcorr_iter(ishift) = map_iter*circshift(mapbase_iter,xshift)';
                                end
                                
%                                 if sum(isnan(fieldXcorr_iter)) == 0
%                                     fieldXcorr_iter = interp1(xorig,fieldXcorr_iter,xinterp,'spline');
%                                 else
%                                     fieldXcorr_iter = interp1(xorig,fieldXcorr_iter,xinterp);
%                                 end
                                maxcorr_iter = getCircularAverage(fieldXcorr_iter(:),0,0.1,0.05);
                                numbins = numel(fieldXcorr_iter);
%                                 [~ , maxcorr_iter] = max(fieldXcorr_iter);
%                                 maxcorr_iter = xinterp(maxcorr_iter) - xinterp(floor(size(fieldXcorr_iter,2)/2) + 1);
                                stdmaxcorr = stdmaxcorr + (kfold - 1)/kfold*numbins/(2*pi)*circ_dist(2*pi/numbins*maxcorr_iter,2*pi/numbins*maxcorr_all).^2;
                            end
                            stdmaxcorr = (stdmaxcorr).^0.5;
                            
                            maxcorr{g}(icell) = maxcorr_all;
                            maxcorrstd{g}(icell) = stdmaxcorr;
                        end
                    end
                    cellprop.fieldXcorrMax{g} = maxcorr{g};
                    cellprop.fieldXcorrMaxSE{g} = maxcorrstd{g};
                    
                    if ~isempty(Slfpspkcoherence.resCA1V1)
                        if isempty(Slfpspkcoherence.resCA1V1(2).spk_CohSpecChAll)
                            lfpspkcoherence = Slfpspkcoherence.resCA1V1(1).spk_CohSpecChAll';
                            lfpPhsspkcoherence = Slfpspkcoherence.resCA1V1(1).spk_PhsCohSpecChAll';
                            lfpspkthetaPhscoherence = circ_mean(lfpPhsspkcoherence(:,Slfpspkcoherence.resCA1V1(1).f>=6 & Slfpspkcoherence.resCA1V1(1).f<=9),[],2);
                        else
                            lfpspkcoherence = (cat(2,Slfpspkcoherence.resCA1V1(1).spk_CohSpecChAll,Slfpspkcoherence.resCA1V1(2).spk_CohSpecChAll))';
                            lfpPhsspkcoherence = (cat(2,Slfpspkcoherence.resCA1V1(1).spk_PhsCohSpecChAll,Slfpspkcoherence.resCA1V1(2).spk_PhsCohSpecChAll))';
                            lfpspkthetaPhscoherence = circ_mean(lfpPhsspkcoherence(:,Slfpspkcoherence.resCA1V1(1).f>=6 & Slfpspkcoherence.resCA1V1(1).f<=9),[],2);
                        end
                        cellprop.lfpSpkCoherence{g} = lfpspkcoherence;
                        cellprop.lfpSpkThetaPhsCoherence{g} = lfpspkthetaPhscoherence;
                    end
                end
            else
                cellprop.animal = [cellprop.animal ianimal*ones(1,S.CellInfo.NumCells)];
                cellprop.iseries = [cellprop.iseries iseries*ones(1,S.CellInfo.NumCells)];
                cellprop.CellList = [cellprop.CellList S.CellInfo.CellList];
                cellprop.Shank = [cellprop.Shank S.CellInfo.Shank];
                cellprop.MUAcluster = [cellprop.MUAcluster S.CellInfo.MUAcluster];
                cellprop.Goodcluster = [cellprop.Goodcluster S.CellInfo.Goodcluster];
                cellprop.Finterneuron = [cellprop.Finterneuron S.CellInfo.Finterneuron];
                cellprop.Unsortedcluster = [cellprop.Unsortedcluster S.CellInfo.Unsortedcluster];
                if isfield(S.CellInfo,'rateautocorr')
                    cellprop.Spatialmodulation = [cellprop.Spatialmodulation S.CellInfo.Spatialmodulation];
                    cellprop.burstindex = [cellprop.burstindex S.CellInfo.burstindex];
                    cellprop.minFR = [cellprop.minFR S.CellInfo.minFR];
                    cellprop.rateautocorr = [cellprop.rateautocorr S.CellInfo.rateautocorr];
                    cellprop.min2maxSpkwf = [cellprop.min2maxSpkwf S.CellInfo.min2maxSpkwf];
                    cellprop.spkautocorr = cat(1,cellprop.spkautocorr,S.CellInfo.spkautocorr);
                end
                cellprop.Probe = [cellprop.Probe S.CellInfo.Probe];
                cellprop.bestchan = [cellprop.bestchan S.CellInfo.bestchan];
                if strcmp(expt(ianimal).area{iseries},'V1medial')
                    cellprop.Cellpos2p = [cellprop.Cellpos2p 1*ones(1,S.CellInfo.NumCells)];
                elseif strcmp(expt(ianimal).area{iseries},'V1lateral')
                    cellprop.Cellpos2p = [cellprop.Cellpos2p 2*ones(1,S.CellInfo.NumCells)];
                elseif strcmp(expt(ianimal).area{iseries},'PPC')
                    cellprop.Cellpos2p = [cellprop.Cellpos2p 3*ones(1,S.CellInfo.NumCells)];
                elseif strcmp(expt(ianimal).area{iseries},'AL')
                    cellprop.Cellpos2p = [cellprop.Cellpos2p 4*ones(1,S.CellInfo.NumCells)];
                else
                    cellprop.Cellpos2p = [cellprop.Cellpos2p zeros(1,S.CellInfo.NumCells)];
                end
                if ~isempty(S.CellInfo.RFXpos)
                    cellprop.RFXpos = cat(1,cellprop.RFXpos,S.CellInfo.RFXpos);
                    cellprop.Xpos = [cellprop.Xpos S.CellInfo.Xpos];
                    cellprop.XposZmax = [cellprop.XposZmax S.CellInfo.XposZmax];
                    cellprop.globalXposrep = cat(1,cellprop.globalXposrep,S.CellInfo.globalXposrep);
                    [~,xposmax] = max(mean(S.CellInfo.RFXpos(S.CellInfo.XposZmax>2 & S.CellInfo.Probe == 2,:),1));
                    cellprop.XposPop = [cellprop.XposPop xposmax*ones(1,S.CellInfo.NumCells)];

                    cellprop.RFYpos = cat(1,cellprop.RFYpos,S.CellInfo.RFYpos);
                    cellprop.Ypos = [cellprop.Ypos S.CellInfo.Ypos];
                    cellprop.YposZmax = [cellprop.YposZmax S.CellInfo.YposZmax];
                    cellprop.globalYposrep = cat(1,cellprop.globalYposrep,S.CellInfo.globalYposrep);
                    [~,yposmax] = max(mean(S.CellInfo.RFYpos(S.CellInfo.XposZmax>2 & S.CellInfo.Probe == 2,:),1));
                    cellprop.YposPop = [cellprop.YposPop yposmax*ones(1,S.CellInfo.NumCells)];
                else
                    cellprop.RFXpos = cat(1,cellprop.RFXpos,NaN(S.CellInfo.NumCells,nbXpos));
                    cellprop.Xpos = [cellprop.Xpos NaN(1,S.CellInfo.NumCells)];
                    cellprop.XposZmax = [cellprop.XposZmax NaN(1,S.CellInfo.NumCells)];
                    cellprop.globalXposrep = cat(1,cellprop.globalXposrep,NaN(S.CellInfo.NumCells,nVStimepts));
                    cellprop.XposPop = [cellprop.XposPop NaN(1,S.CellInfo.NumCells)];

                    cellprop.RFYpos = cat(1,cellprop.RFYpos,NaN(S.CellInfo.NumCells,nbYpos));
                    cellprop.Ypos = [cellprop.Ypos NaN(1,S.CellInfo.NumCells)];
                    cellprop.YposZmax = [cellprop.YposZmax NaN(1,S.CellInfo.NumCells)];
                    cellprop.globalYposrep = cat(1,cellprop.globalYposrep,NaN(S.CellInfo.NumCells,nVStimepts));
                    cellprop.YposPop = [cellprop.YposPop NaN(1,S.CellInfo.NumCells)];
                end
                
                for iprobe = 1:2
                    if ~isempty(S.CellInfo.LFP2Spike_phscorrMUA{iprobe})
                        cellprop.phsfieldMUA{iprobe} = cat(1,cellprop.phsfieldMUA{iprobe},S.CellInfo.phsfieldMUA{iprobe}');
                        cellprop.phsfieldMUASE{iprobe} = cat(1,cellprop.phsfieldMUASE{iprobe},S.CellInfo.phsfieldMUASE{iprobe}');
                        cellprop.phsfieldZMUA{iprobe} = cat(1,cellprop.phsfieldZMUA{iprobe},S.CellInfo.phsfieldZMUA{iprobe});
                        cellprop.phsfieldPosMUA{iprobe} = cat(1,cellprop.phsfieldPosMUA{iprobe},S.CellInfo.phsfieldPosMUA{iprobe});
                        cellprop.LFP2Spike_phscorrMUA{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMUA{iprobe},S.CellInfo.LFP2Spike_phscorrMUA{iprobe});
                        
                        cellprop.phsfieldMUAnorm{iprobe} = cat(1,cellprop.phsfieldMUAnorm{iprobe},S.CellInfo.phsfieldMUAnorm{iprobe});
                        cellprop.phsfieldPosMUAnorm{iprobe} = cat(1,cellprop.phsfieldPosMUAnorm{iprobe},S.CellInfo.phsfieldPosMUAnorm{iprobe});
                        cellprop.LFP2Spike_phscorrMUAnorm{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMUAnorm{iprobe},S.CellInfo.LFP2Spike_phscorrMUAnorm{iprobe});
                        
                        cellprop.phsfieldMean{iprobe} = cat(1,cellprop.phsfieldMean{iprobe},S.CellInfo.phsfieldMean{iprobe});
                        cellprop.phsfieldPosMean{iprobe} = cat(1,cellprop.phsfieldPosMean{iprobe},S.CellInfo.phsfieldPosMean{iprobe});
                        cellprop.LFP2Spike_phscorrMean{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMean{iprobe},S.CellInfo.LFP2Spike_phscorrMean{iprobe});
                    else
                        cellprop.phsfieldMUA{iprobe} = cat(2,cellprop.phsfieldMUA{iprobe},NaN(size(S.CellInfo.phsfieldMUA{1}')));
                        cellprop.phsfieldMUASE{iprobe} = cat(2,cellprop.phsfieldMUASE{iprobe},NaN(size(S.CellInfo.phsfieldMUASE{1}')));
                        cellprop.phsfieldZMUA{iprobe} = cat(1,cellprop.phsfieldZMUA{iprobe},NaN(size(S.CellInfo.phsfieldZMUA{1})));
                        cellprop.phsfieldPosMUA{iprobe} = cat(1,cellprop.phsfieldPosMUA{iprobe},NaN(size(S.CellInfo.phsfieldPosMUA{1})));
                        cellprop.LFP2Spike_phscorrMUA{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMUA{iprobe},NaN(size(S.CellInfo.LFP2Spike_phscorrMUA{1})));
                        
                        cellprop.phsfieldMUAnorm{iprobe} = cat(2,cellprop.phsfieldMUAnorm{iprobe},NaN(size(S.CellInfo.phsfieldMUAnorm{1})));
                        cellprop.phsfieldPosMUAnorm{iprobe} = cat(1,cellprop.phsfieldPosMUAnorm{iprobe},NaN(size(S.CellInfo.phsfieldPosMUAnorm{1})));
                        cellprop.LFP2Spike_phscorrMUAnorm{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMUAnorm{iprobe},NaN(size(S.CellInfo.LFP2Spike_phscorrMUAnorm{1})));
                        
                        cellprop.phsfieldMean{iprobe} = cat(2,cellprop.phsfieldMean{iprobe},NaN(size(S.CellInfo.phsfieldMean{1})));
                        cellprop.phsfieldPosMean{iprobe} = cat(1,cellprop.phsfieldPosMean{iprobe},NaN(size(S.CellInfo.phsfieldPosMean{1})));
                        cellprop.LFP2Spike_phscorrMean{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMean{iprobe},NaN(size(S.CellInfo.LFP2Spike_phscorrMean{1})));
                    end
                end
                
                allcontidx = size(S.CellInfo.field,1);
                for g = 1:3
                    if g <= size(S.CellInfo.field,2)
                        cellprop.field{g} = cat(1,cellprop.field{g}, S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldSE{g} = cat(1,cellprop.fieldSE{g}, S.CellInfo.fieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldZ{g} = cat(2,cellprop.fieldZ{g}, S.CellInfo.fieldZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldXcorr{g} = cat(1,cellprop.fieldXcorr{g}, S.CellInfo.fieldXcorr{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldPos{g} = cat(2,cellprop.fieldPos{g}, S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldPosSE{g} = cat(2,cellprop.fieldPosSE{g}, S.CellInfo.fieldPosSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldCOM{g} = cat(2,cellprop.fieldCOM{g}, S.CellInfo.fieldCOM{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.rate{g} = cat(2,cellprop.rate{g}, S.CellInfo.rate{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldMax{g} = cat(2,cellprop.fieldMax{g}, S.CellInfo.fieldMax{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldMin{g} = cat(2,cellprop.fieldMin{g}, S.CellInfo.fieldMin{allcontidx,g,1,S.CellInfo.outcomeVal == 2});

                        cellprop.SSI{g} = cat(2,cellprop.SSI{g}, S.CellInfo.SSI{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.SpatialInfo{g} = cat(2,cellprop.SpatialInfo{g}, S.CellInfo.SpatialInfo{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.SpatialInfoPerSpike{g} = cat(2,cellprop.SpatialInfoPerSpike{g}, S.CellInfo.SpatialInfoPerSpike{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
%                         cellprop.thetaPhaseMax{g} = cat(2,cellprop.thetaPhaseMax{g},S.CellInfo.PhaseMax{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXPhstheta{g} = cat(1,cellprop.field2dXPhstheta{g},S.CellInfo.field2dXPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXPhsthetaSE{g} = cat(1,cellprop.field2dXPhsthetaSE{g},S.CellInfo.field2dXPhsthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXthetapos{g} = cat(1,cellprop.field2dXthetapos{g},S.CellInfo.field2dXthetapos{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXthetaposSE{g} = cat(1,cellprop.field2dXthetaposSE{g},S.CellInfo.field2dXthetaposSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXthetaposNorm{g} = cat(1,cellprop.field2dXthetaposNorm{g},S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXthetaposSENorm{g} = cat(1,cellprop.field2dXthetaposSENorm{g},S.CellInfo.field2dXthetaposSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        cellprop.field2dXcorrtheta{g} = cat(1,cellprop.field2dXcorrtheta{g},S.CellInfo.field2dXcorrtheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetaSE{g} = cat(1,cellprop.field2dXcorrthetaSE{g},S.CellInfo.field2dXcorrthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetamax{g} = cat(1,cellprop.field2dXcorrthetamax{g},S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetamaxSE{g} = cat(1,cellprop.field2dXcorrthetamaxSE{g},S.CellInfo.field2dXcorrthetamaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        thetafieldset1 = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                        thetafieldset2 = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                        cellprop.field2dXcorrthetamax_set1{g} = cat(1,cellprop.field2dXcorrthetamax_set1{g},thetafieldset1);
                        cellprop.field2dXcorrthetamax_set2{g} = cat(1,cellprop.field2dXcorrthetamax_set2{g},thetafieldset2);
                        thetareliabilityCorr = sum((thetafieldset1-repmat(nanmean(thetafieldset1,2),[1 size(thetafieldset1,2)])).*...
                                               (thetafieldset2-repmat(nanmean(thetafieldset2,2),[1 size(thetafieldset2,2)])),2)./...
                                               sqrt(sum((thetafieldset1-repmat(nanmean(thetafieldset1,2),[1 size(thetafieldset1,2)])).^2,2).*...
                                               sum((thetafieldset2-repmat(nanmean(thetafieldset2,2),[1 size(thetafieldset2,2)])).^2,2));
                        cellprop.thetareliabilityCorr{g} = cat(1,cellprop.thetareliabilityCorr{g},thetareliabilityCorr);                        
                        
                        cellprop.phsfieldPos{g} = cat(2,cellprop.phsfieldPos{g},S.CellInfo.phsfieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldPosSE{g} = cat(2,cellprop.phsfieldPosSE{g},S.CellInfo.phsfieldPosSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldCOM{g} = cat(2,cellprop.phsfieldCOM{g},S.CellInfo.phsfieldCOM{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldCOMSE{g} = cat(2,cellprop.phsfieldCOMSE{g},S.CellInfo.phsfieldCOMSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfield{g} = cat(1,cellprop.phsfield{g},S.CellInfo.phsfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldSE{g} = cat(1,cellprop.phsfieldSE{g},S.CellInfo.phsfieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldZ{g} = cat(2,cellprop.phsfieldZ{g},S.CellInfo.phsfieldZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsmodulation{g} = cat(2,cellprop.phsmodulation{g},S.CellInfo.phsmodulation{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        ZXposthetaMax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        ZXposthetaMin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        XposthetaMax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        XposthetaMin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        ZXPhsthetaMax = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        ZXPhsthetaMin = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        PhsthetaMod = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        Phsthetacross = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        phaseXmin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        phaseXmax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        phaseXahead = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        phaseXbehind = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        phaseXcross = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        nphsbins = size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},2);
                        for icell = 1:size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
                            %                         Xpostheta = S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            %                         XposthetaSE = S.CellInfo.field2dXthetaposSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            Xpostheta = S.CellInfo.field2dXcorrthetamaxNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            XposthetaSE = S.CellInfo.field2dXcorrthetamaxSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            [~, imax] = max(Xpostheta);
                            ZXposthetaMax(icell) = abs((Xpostheta(imax) - mean(Xpostheta))./XposthetaSE(imax));
                            [~, imin] = min(Xpostheta);
                            ZXposthetaMin(icell) = abs((Xpostheta(imin) - mean(Xpostheta))./XposthetaSE(imin));
                            
                            Xpostheta = S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            [XposthetaMax(icell), ~] = max(Xpostheta);
                            [XposthetaMin(icell), ~] = min(Xpostheta);
                            
                            Phstheta = S.CellInfo.field2dPhsthetaNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            PhsthetaSE = S.CellInfo.field2dPhsthetaSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            [~, imax] = max(Phstheta);
                            ZXPhsthetaMax(icell) = abs((Phstheta(imax) - mean(Phstheta))./PhsthetaSE(imax));
                            [~, imin] = min(Phstheta);
                            ZXPhsthetaMin(icell) = abs((Phstheta(imin) - mean(Phstheta))./PhsthetaSE(imin));
%                             Phsthetacross(icell) = (imax-1)*360/nphsbins;
                            
                            Phstheta = S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:)./mean(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:));
                            Phi = (0:(nphsbins-1))/nphsbins*2*pi;
                            Phi_num = sum(Phstheta.*sin(Phi));
                            Phi_den = sum(Phstheta.*cos(Phi));
                            Phsthetacross(icell) = mod(90 - 360*atan2(Phi_den,Phi_num)/(2*pi),360);
                            [PhsthetaMod(icell), ~] = max(Phstheta);
                            
                            Xpostheta = squeeze(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:));
                            imax = getCircularAverage(Xpostheta',0,0.01,0.05);%[~, imax] = max(Xpostheta);
                            phaseXmax(icell) = imax;
                            imin = getCircularAverage(-Xpostheta',0,0.01,0.05);%[~, imin] = min(Xpostheta);
                            phaseXmin(icell) = imin;
                            
                            Xpostheta_ahead = Xpostheta - mean(Xpostheta);
                            Xpostheta_ahead(Xpostheta_ahead>0) = 0;
                            phaseXahead(icell) = getCircularAverage(abs(Xpostheta_ahead)',0,1);
                            Xpostheta_behind = Xpostheta - mean(Xpostheta);
                            Xpostheta_behind(Xpostheta_behind<0) = 0;
                            phaseXbehind(icell) = getCircularAverage(abs(Xpostheta_behind)',0,1);
                            
                            Phi = (0:(nphsbins-1))/nphsbins*2*pi;
                            Phi_num = sum(Xpostheta.*sin(Phi));
                            Phi_den = sum(Xpostheta.*cos(Phi));
                            phaseXcross(icell) = 180 - 360*atan2(Phi_den,Phi_num)/(2*pi);
                        end
                        cellprop.ZXmaxthetapos{g} = cat(1,cellprop.ZXmaxthetapos{g},ZXposthetaMax);
                        cellprop.ZXminthetapos{g} = cat(1,cellprop.ZXminthetapos{g},ZXposthetaMin);
                        cellprop.Xmaxthetapos{g} = cat(1,cellprop.Xmaxthetapos{g},XposthetaMax);
                        cellprop.Xminthetapos{g} = cat(1,cellprop.Xminthetapos{g},XposthetaMin);
                        cellprop.ZPhsmaxtheta{g} = cat(1,cellprop.ZPhsmaxtheta{g},ZXPhsthetaMax);
                        cellprop.ZPhsmintheta{g} = cat(1,cellprop.ZPhsmintheta{g},ZXPhsthetaMin);
                        cellprop.thetaPhaseMax{g} = cat(1,cellprop.thetaPhaseMax{g},Phsthetacross);
                        cellprop.PhsthetaMod{g} = cat(1,cellprop.PhsthetaMod{g},PhsthetaMod);
                        cellprop.PhsXmaxthetapos{g} = cat(1,cellprop.PhsXmaxthetapos{g},phaseXmax);
                        cellprop.PhsXminthetapos{g} = cat(1,cellprop.PhsXminthetapos{g},phaseXmin);
                        cellprop.PhsXaheadthetapos{g} = cat(1,cellprop.PhsXaheadthetapos{g},phaseXahead);
                        cellprop.PhsXbehindthetapos{g} = cat(1,cellprop.PhsXbehindthetapos{g},phaseXbehind);
                        cellprop.PhsXcrossthetapos{g} = cat(1,cellprop.PhsXcrossthetapos{g},phaseXcross);
                        
%                         ZXposthetaMax = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                         ZXposthetaMin = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                         phaseXmin = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                         phaseXmax = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                         for icell = 1:size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                             Xpostheta = S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
%                             XposthetaSE = S.CellInfo.field2dXthetaposSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
%                             [~, imax] = max(Xpostheta);
%                             phaseXmax(icell) = imax;
%                             ZXposthetaMax(icell) = abs((Xpostheta(imax) - mean(Xpostheta))./XposthetaSE(imax));
%                             [~, imin] = min(Xpostheta);
%                             ZXposthetaMin(icell) = abs((Xpostheta(imin) - mean(Xpostheta))./XposthetaSE(imin));
%                             phaseXmin(icell) = imin;
%                         end
%                         cellprop.ZXmaxthetaposNorm{g} = cat(1,cellprop.ZXmaxthetaposNorm{g},ZXposthetaMax);
%                         cellprop.ZXminthetaposNorm{g} = cat(1,cellprop.ZXminthetaposNorm{g},ZXposthetaMin);
%                         cellprop.PhsXmaxthetaposNorm{g} = cat(1,cellprop.PhsXmaxthetaposNorm{g},phaseXmax);
%                         cellprop.PhsXminthetaposNorm{g} =cat(1,cellprop.PhsXminthetaposNorm{g},phaseXmin);
                        
                        cellprop.PhaseRho{g} = cat(2,cellprop.PhaseRho{g},S.CellInfo.PhaseRho{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.PhasePval{g} = cat(2,cellprop.PhasePval{g},S.CellInfo.PhasePval{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dPhstheta{g} = cat(1,cellprop.field2dPhstheta{g},S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.PhaseRayleighPval{g} = cat(2,cellprop.PhaseRayleighPval{g},S.CellInfo.PhaseRayleighPval{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.PhaseRayleighZ{g} = cat(2,cellprop.PhaseRayleighZ{g},S.CellInfo.PhaseRayleighZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        cellprop.field_half1{g} = cat(1,cellprop.field_half1{g}, S.CellInfo.field_half1{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field_half2{g} = cat(1,cellprop.field_half2{g}, S.CellInfo.field_half2{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        fieldset1 = squeeze(S.CellInfo.field_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                        fieldset2 = squeeze(S.CellInfo.field_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                        cellprop.field_set1{g} = cat(1,cellprop.field_set1{g}, fieldset1);
                        cellprop.field_set2{g} = cat(1,cellprop.field_set2{g}, fieldset2);
                        reliabilityCorr = sum((fieldset1-repmat(nanmean(fieldset1,2),[1 size(fieldset1,2)])).*...
                                               (fieldset2-repmat(nanmean(fieldset2,2),[1 size(fieldset2,2)])),2)./...
                                               sqrt(sum((fieldset1-repmat(nanmean(fieldset1,2),[1 size(fieldset1,2)])).^2,2).*...
                                               sum((fieldset2-repmat(nanmean(fieldset2,2),[1 size(fieldset2,2)])).^2,2));
                        cellprop.reliabilityCorr{g} = cat(1,cellprop.reliabilityCorr{g},reliabilityCorr);
                        
                        
%                         maxpos = [];
%                         for icell = 1:size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                             map = repmat(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:),[1 3]);
%                             xorig = 0.5:size(map,2)-0.5;
%                             xinterp = 0.1:0.1:size(map,2);
%                             mapinterp = interp1(xorig,map,xinterp,'spline');
%                             mapinterp = mapinterp((numel(mapinterp)/3+1):(2*numel(mapinterp)/3));
%                             [~,imax] = max(mapinterp);
%                             maxpos(icell) = xinterp(imax);
%                         end
%                         cellprop.fieldPos{g} = cat(2,cellprop.fieldPos{g},maxpos(:)');
                        
                        maxcorr{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                        maxcorrstd{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                        if ~isempty(Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model)
                            for icell = 1:numel(Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning)
                                map = Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).meanrespModel;
                                mapbase = Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).meanrespModel;
                                
                                map = map - nanmean(map);
                                map = map./sqrt(sum(map.^2));
                                mapbase = mapbase - nanmean(mapbase);
                                mapbase = mapbase./sqrt(sum(mapbase.^2));
                                fieldXcorr_all = zeros(1,numel(mapbase));
                                xshiftlim = floor(numel(mapbase)/2);%20;
                                ishift = 0;%floor(numel(mapbase)/2)-xshiftlim-1;
                                for xshift = -xshiftlim:xshiftlim
                                    ishift = ishift + 1;
                                    fieldXcorr_all(ishift) = map*circshift(mapbase,xshift)';
                                end
                                
                                %                             xorig = 0.5:numel(fieldXcorr_all)-0.5;
                                %                             xinterp = 0:0.1:numel(fieldXcorr_all);
                                %
                                %                             if sum(isnan(fieldXcorr_all)) == 0
                                %                                 fieldXcorr_all = interp1(xorig,fieldXcorr_all,xinterp,'spline');
                                %                             else
                                %                                 fieldXcorr_all = interp1(xorig,fieldXcorr_all,xinterp);
                                %                             end
                                maxcorr_all = getCircularAverage(fieldXcorr_all(:),0,1);
                                %                             [~ , maxcorr_all] = max(fieldXcorr_all);
                                %                             maxcorr_all = xinterp(maxcorr_all) - xinterp(floor(size(fieldXcorr_all,2)/2) + 1);
                                
                                stdmaxcorr = 0;
                                kfold = size(Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel,1);
                                %SEM computed using Jacknife method
                                for i = 1:kfold
                                    map_iter = Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel(i,:);%map;%
                                    mapbase_iter = mapbase;%Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel(i,:);%
                                    
                                    map_iter = map_iter - nanmean(map_iter);
                                    map_iter = map_iter./sqrt(sum(map_iter.^2));
                                    mapbase_iter = mapbase_iter - nanmean(mapbase_iter);
                                    mapbase_iter = mapbase_iter./sqrt(sum(mapbase_iter.^2));
                                    fieldXcorr_iter = zeros(1,numel(mapbase_iter));
                                    xshiftlim = floor(numel(mapbase_iter)/2);%20;
                                    ishift = 0;%floor(numel(mapbase_iter)/2)-xshiftlim-1;
                                    for xshift = -xshiftlim:xshiftlim
                                        ishift = ishift + 1;
                                        fieldXcorr_iter(ishift) = map_iter*circshift(mapbase_iter,xshift)';
                                    end
                                    
                                    %                                 if sum(isnan(fieldXcorr_iter)) == 0
                                    %                                     fieldXcorr_iter = interp1(xorig,fieldXcorr_iter,xinterp,'spline');
                                    %                                 else
                                    %                                     fieldXcorr_iter = interp1(xorig,fieldXcorr_iter,xinterp);
                                    %                                 end
                                    maxcorr_iter = getCircularAverage(fieldXcorr_iter(:),0,1);
                                    %                                 [~ , maxcorr_iter] = max(fieldXcorr_iter);
                                    %                                 maxcorr_iter = xinterp(maxcorr_iter) - xinterp(floor(size(fieldXcorr_iter,2)/2) + 1);
                                    stdmaxcorr = stdmaxcorr + (kfold - 1)/kfold*(maxcorr_iter - maxcorr_all).^2;
                                end
                                stdmaxcorr = (stdmaxcorr).^0.5;
                                
                                maxcorr{g}(icell) = maxcorr_all;
                                maxcorrstd{g}(icell) = stdmaxcorr;
                            end
                        end
                        cellprop.fieldXcorrMax{g} = cat(2, cellprop.fieldXcorrMax{g}, maxcorr{g});
                        cellprop.fieldXcorrMaxSE{g} = cat(2, cellprop.fieldXcorrMaxSE{g}, maxcorrstd{g});
                        
                        if ~isempty(Slfpspkcoherence.resCA1V1)
                            if isempty(Slfpspkcoherence.resCA1V1(2).spk_CohSpecChAll)
                                lfpspkcoherence = Slfpspkcoherence.resCA1V1(1).spk_CohSpecChAll';
                                lfpPhsspkcoherence = Slfpspkcoherence.resCA1V1(1).spk_PhsCohSpecChAll';
                                lfpspkthetaPhscoherence = circ_mean(lfpPhsspkcoherence(:,Slfpspkcoherence.resCA1V1(1).f>=6 & Slfpspkcoherence.resCA1V1(1).f<=9),[],2);
                            else
                                lfpspkcoherence = (cat(2,Slfpspkcoherence.resCA1V1(1).spk_CohSpecChAll,Slfpspkcoherence.resCA1V1(2).spk_CohSpecChAll))';
                                lfpPhsspkcoherence = (cat(2,Slfpspkcoherence.resCA1V1(1).spk_PhsCohSpecChAll,Slfpspkcoherence.resCA1V1(2).spk_PhsCohSpecChAll))';
                                lfpspkthetaPhscoherence = circ_mean(lfpPhsspkcoherence(:,Slfpspkcoherence.resCA1V1(1).f>=6 & Slfpspkcoherence.resCA1V1(1).f<=9),[],2);
                            end
                            cellprop.lfpSpkCoherence{g} = cat(1,cellprop.lfpSpkCoherence{g},lfpspkcoherence);
                            cellprop.lfpSpkThetaPhsCoherence{g} = cat(1,cellprop.lfpSpkThetaPhsCoherence{g},lfpspkthetaPhscoherence);
                        end
                    else
                        cellprop.field{g} = cat(1,cellprop.field{g}, NaN(size(S.CellInfo.field{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldSE{g} = cat(1,cellprop.fieldSE{g}, NaN(size(S.CellInfo.fieldSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldZ{g} = cat(2,cellprop.fieldZ{g}, NaN(size(S.CellInfo.fieldZ{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldXcorr{g} = cat(1,cellprop.fieldXcorr{g}, NaN(size(S.CellInfo.fieldXcorr{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldPos{g} = cat(2,cellprop.fieldPos{g}, NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldPosSE{g} = cat(2,cellprop.fieldPosSE{g}, NaN(size(S.CellInfo.fieldPosSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldCOM{g} = cat(2,cellprop.fieldCOM{g}, NaN(size(S.CellInfo.fieldCOM{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.rate{g} = cat(2,cellprop.rate{g}, NaN(size(S.CellInfo.rate{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldMax{g} = cat(2,cellprop.fieldMax{g}, NaN(size(S.CellInfo.fieldMax{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldMin{g} = cat(2,cellprop.fieldMin{g}, NaN(size(S.CellInfo.fieldMin{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));

                        cellprop.SSI{g} = cat(2,cellprop.SSI{g}, NaN(size(S.CellInfo.SSI{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.SpatialInfo{g} = cat(2,cellprop.SpatialInfo{g}, NaN(size(S.CellInfo.SpatialInfo{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.SpatialInfoPerSpike{g} = cat(2,cellprop.SpatialInfoPerSpike{g}, NaN(size(S.CellInfo.SpatialInfoPerSpike{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.thetaPhaseMax{g} = cat(1,cellprop.thetaPhaseMax{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.field2dXPhstheta{g} = cat(1,cellprop.field2dXPhstheta{g},NaN(size(S.CellInfo.field2dXPhstheta{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXPhsthetaSE{g} = cat(1,cellprop.field2dXPhsthetaSE{g},NaN(size(S.CellInfo.field2dXPhsthetaSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXthetapos{g} = cat(1,cellprop.field2dXthetapos{g},NaN(size(S.CellInfo.field2dXthetapos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXthetaposSE{g} = cat(1,cellprop.field2dXthetaposSE{g},NaN(size(S.CellInfo.field2dXthetaposSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXthetaposNorm{g} = cat(1,cellprop.field2dXthetaposNorm{g},NaN(size(S.CellInfo.field2dXthetaposNorm{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXthetaposSENorm{g} = cat(1,cellprop.field2dXthetaposSENorm{g},NaN(size(S.CellInfo.field2dXthetaposSENorm{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dXcorrtheta{g} = cat(1,cellprop.field2dXcorrtheta{g},NaN(size(S.CellInfo.field2dXcorrtheta{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaSE{g} = cat(1,cellprop.field2dXcorrthetaSE{g},NaN(size(S.CellInfo.field2dXcorrthetaSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetamax{g} = cat(1,cellprop.field2dXcorrthetamax{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetamaxSE{g} = cat(1,cellprop.field2dXcorrthetamaxSE{g},NaN(size(S.CellInfo.field2dXcorrthetamaxSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dXcorrthetamax_set1{g} = cat(1,cellprop.field2dXcorrthetamax_set1{g},NaN(size(squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,2,1,S.CellInfo.outcomeVal == 2}(:,1,:)))));
                        cellprop.field2dXcorrthetamax_set2{g} = cat(1,cellprop.field2dXcorrthetamax_set2{g},NaN(size(squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,2,1,S.CellInfo.outcomeVal == 2}(:,2,:)))));
                        
                        cellprop.phsfieldPos{g} = cat(2,cellprop.phsfieldPos{g},NaN(size(S.CellInfo.phsfieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldPosSE{g} = cat(2,cellprop.phsfieldPosSE{g},NaN(size(S.CellInfo.phsfieldPosSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldCOM{g} = cat(2,cellprop.phsfieldCOM{g},NaN(size(S.CellInfo.phsfieldCOM{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldCOMSE{g} = cat(2,cellprop.phsfieldCOMSE{g},NaN(size(S.CellInfo.phsfieldCOMSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfield{g} = cat(1,cellprop.phsfield{g},NaN(size(S.CellInfo.phsfield{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldSE{g} = cat(1,cellprop.phsfieldSE{g},NaN(size(S.CellInfo.phsfieldSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldZ{g} = cat(2,cellprop.phsfieldZ{g},NaN(size(S.CellInfo.phsfieldZ{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsmodulation{g} = cat(2,cellprop.phsmodulation{g},NaN(size(S.CellInfo.phsmodulation{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.ZXmaxthetapos{g} = cat(1,cellprop.ZXmaxthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.ZXminthetapos{g} = cat(1,cellprop.ZXminthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.Xmaxthetapos{g} = cat(1,cellprop.Xmaxthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.Xminthetapos{g} = cat(1,cellprop.Xminthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.ZPhsmaxtheta{g} = cat(1,cellprop.ZPhsmaxtheta{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.ZPhsmintheta{g} = cat(1,cellprop.ZPhsmintheta{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsthetaMod{g} = cat(1,cellprop.PhsthetaMod{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsXmaxthetapos{g} = cat(1,cellprop.PhsXmaxthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsXminthetapos{g} = cat(1,cellprop.PhsXminthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsXaheadthetapos{g} = cat(1,cellprop.PhsXaheadthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsXbehindthetapos{g} = cat(1,cellprop.PhsXbehindthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsXcrossthetapos{g} = cat(1,cellprop.PhsXcrossthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        
                        cellprop.PhaseRho{g} = cat(2,cellprop.PhaseRho{g},NaN(size(S.CellInfo.PhaseRho{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.PhasePval{g} = cat(2,cellprop.PhasePval{g},NaN(size(S.CellInfo.PhasePval{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dPhstheta{g} = cat(1,cellprop.field2dPhstheta{g},NaN(size(S.CellInfo.field2dPhstheta{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.PhaseRayleighPval{g} = cat(2,cellprop.PhaseRayleighPval{g},NaN(size(S.CellInfo.PhaseRayleighPval{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.PhaseRayleighZ{g} = cat(2,cellprop.PhaseRayleighZ{g},NaN(size(S.CellInfo.PhaseRayleighZ{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field_half1{g} = cat(1,cellprop.field_half1{g}, NaN(size(S.CellInfo.field_half1{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field_half2{g} = cat(1,cellprop.field_half2{g}, NaN(size(S.CellInfo.field_half2{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field_set1{g} = cat(1,cellprop.field_set1{g}, NaN(size(S.CellInfo.field_set1{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field_set2{g} = cat(1,cellprop.field_set2{g}, NaN(size(S.CellInfo.field_set1{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        reliabilityCorr = NaN(numel(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2}),1);
                        cellprop.reliabilityCorr{g} = cat(1,cellprop.reliabilityCorr{g},reliabilityCorr);
                        cellprop.thetareliabilityCorr{g} = cat(1,cellprop.thetareliabilityCorr{g},reliabilityCorr);
                        
                        maxcorr{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2}));
                        maxcorrstd{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2}));
                        cellprop.fieldXcorrMax{g} = cat(2, cellprop.fieldXcorrMax{g}, maxcorr{g});
                        cellprop.fieldXcorrMaxSE{g} = cat(2, cellprop.fieldXcorrMaxSE{g}, maxcorrstd{g});
                        
                        lfpspkcoherence = NaN(size(Slfpspkcoherence.resCA1V1(1).spk_CohSpecChAll'));
                        lfpspkthetaPhscoherence = NaN(size(lfpspkcoherence,1),1);
                        cellprop.lfpSpkCoherence{g} = cat(1,cellprop.lfpSpkCoherence{g},lfpspkcoherence);
                        cellprop.lfpSpkThetaPhsCoherence{g} = cat(1,cellprop.lfpSpkThetaPhsCoherence{g},lfpspkthetaPhscoherence);
                    end
                end
            end
        end
    end
end
end