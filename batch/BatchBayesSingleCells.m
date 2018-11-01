function EXP = BatchBayesSingleCells(EXP, DIRS)
strlistvarname = {'2p data','electrophys data'};
[varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'dataset', 'SelectionMode', 'single', 'InitialValue', 1);
if ok && varnamesel == 1
    batch2p = true;
    expt = getExperimentList2p;
    strlistvarname = {'V1','PPC', 'AL'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'area', 'SelectionMode', 'single', 'InitialValue', 1);
    area_str = strlistvarname{varnamesel};
elseif ok
    batch2p = false;
    expt = getExperimentList;
    area_str = 'CA1V1';
end
prompt = {'Speed Threshold';'nthetaphsbins';'nspdbins';'Window size (ms)';'Spatial smth (%)';'# of X bins';'FGoodcluster';'FMUAcluster';'FoverwriteCellinfo'};
dlg_title = 'Parameters';
num_lines = 1;
def = {'1';'0';'5';'150';'4';'100';'1';'0';'1'};
nruns = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(nruns)
    SpeedThreshold = str2num(nruns{1});
    nthetaphsbins = str2num(nruns{2});
    nspeedbins = str2num(nruns{3});
    Tsmthwin = str2num(nruns{4});
    Xsmthwin = str2num(nruns{5});
    nDecbins = str2num(nruns{6});
    FGoodcluster = logical(str2num(nruns{7}));
    FMUAcluster = logical(str2num(nruns{8}));
    FoverwriteCellinfo = logical(str2num(nruns{9}));
end
Fprocessagain = true;

if FGoodcluster && ~FMUAcluster
    cellstr = 'goodonly';
elseif FGoodcluster && FMUAcluster
    cellstr = 'All';
end

nanimals = numel(expt);%1;
Nperm_cellprop = 2;%100

FUnsortedcluster = 0;
maxRate = inf;
zth = -inf;
SSImin = -inf;


filesuffix_EXP = ['Twin' num2str(Tsmthwin) '_' 'Xwin' num2str(Xsmthwin) '_' 'spdth' num2str(SpeedThreshold) '_' num2str(nspeedbins) 'speedbins' '_' num2str(nthetaphsbins) 'thetabins' '_' cellstr];
filesuffix_cellprop = ['Twin' num2str(Tsmthwin) '_' 'Xwin' num2str(Xsmthwin) '_' 'spdth' num2str(SpeedThreshold) '_' num2str(nthetaphsbins) 'thetabins' '_' cellstr];
for ianimal = 1:nanimals
    for iseries = 1:numel(expt(ianimal).series)
        if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
            EXP = TVRData;
            EXP.animal = expt(ianimal).animal;
            EXP.series = num2str(expt(ianimal).series{iseries});
            EXP.iseries = expt(ianimal).series{iseries};
            EXP.exptList = expt(ianimal).exp{iseries};
            if ~batch2p
                dDIRname = ['D:\DATA\batch'  filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];%[DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];
            else
                dDIRname = ['D:\DATA\batch'  filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];%[DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];
            end
            if ~isdir(dDIRname)
                mkdir(dDIRname)
            end
            savedfile = [dDIRname filesep 'EXP_' filesuffix_EXP '.mat'];
            disp([EXP.animal ' series ' num2str(EXP.iseries)]);
            
            if ~exist(savedfile,'file') || Fprocessagain
                if strcmp(EXP.animal,'M160114C_BALL') && EXP.iseries == 323
                    shanknum = [0:15 0];
                    suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
                else
                    shanknum = [0:7 0];
                    suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
                end
                EXP.LoadVRData(shanknum, suffix, SpeedThreshold, nthetaphsbins);
                
                savedcellinfo = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_cellProperties.mat'];
                savedmaps1d = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d.mat'];
                savedVS = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_VS.mat'];
                
                if exist(savedcellinfo) && ~FoverwriteCellinfo
                    disp('loading cellInfo from saved file');
                    S = load(savedcellinfo);
                    EXP.CellInfo = S.CellInfo;
                    S = load(savedVS);
                    EXP.data.VStuning = S.VStuning;
                    EXP.data.VSstim = S.VSstim;
                    S = load(savedmaps1d);
                    EXP.maps1d = S.maps1d;
                else
                    EXP.CalculateStimTuning([], shanknum, suffix);
                    VStuning = EXP.data.VStuning;
                    VSstim = EXP.data.VSstim;
                    save(savedVS,'VStuning','VSstim');
                    delay = 0;
                    EXP.Calculate1Dmaps('trajPercent', Tsmthwin,Xsmthwin,delay);
                    maps1d = EXP.maps1d;
                    save(savedmaps1d,'maps1d');
                    EXP.defineCellInfo(EXP.data.es.spikeIDs, EXP.data.es.chanIDs,EXP.data.es.ProbeIDs);
                    EXP.defineCellProp(Nperm_cellprop);
                    CellInfo = EXP.CellInfo;
                    save(savedcellinfo,'CellInfo');
                end
                
                traincont = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
                traingain = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
                trainroomlength = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
                trainoutcome = find(EXP.SubsetVal.outcome == 2);
                type = 'mean';
                Flookuptable = false;
                nruns = 1;
                goodidx_th = 30;
                Tsmth_win = Tsmthwin;%20;%
                Xsmth_win = Xsmthwin;
                numbins = nDecbins;
                thetachannel = 34;
                nthetabins = 1;%nthetaphsbins;%1;%6;%
                nspdbins = nspeedbins;%1;%
                smth_spd = Tsmth_win;
                speed_th = SpeedThreshold;
                
                kfold = 20;
                FoptiSmooth = 0;
                latcorrection = 0;
                alpha = 0;
                delta = 0;
                nbProbe = numel(unique(EXP.CellInfo.Probe));
                Spiketrain_orig = EXP.data.es.spikeTrain;
                EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
                    'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'smth_spd', smth_spd, 'latcorrection', latcorrection,...
                    'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                    'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin);
                
                if ~exist(savedfile,'file') || Fprocessagain
                    save(savedfile,'EXP','-v7.3');
                end
            else
                S = load(savedfile);
                EXP.Copyobj(S.EXP);
            end
            
            resAll.DecCells = cell(1,nbProbe);
            resAll.Xpred_ave = cell(1,nbProbe);
            resAll.Xpred_max = cell(1,nbProbe);
            for iprobe = 1:nbProbe
                resAll.DecCells{iprobe} = EXP.Bayes.DecCellidx{iprobe};
                amp_th = 1;maxtol = 1;
                resAll.Xpred_ave{iprobe} = getCircularAverage(EXP.Bayes.Posterior0{iprobe}',amp_th,maxtol);
                amp_th = 1;maxtol = 0.1;
                resAll.Xpred_max{iprobe} = getCircularAverage(EXP.Bayes.Posterior0{iprobe}',amp_th,maxtol);
            end
            
            res.Xpred_ave = cell(1,EXP.CellInfo.NumCells);
            res.Xpred_max = cell(1,EXP.CellInfo.NumCells);
            for iprobe = 1:nbProbe
                try
                cellprobeidx = find(EXP.Bayes.DecCellidx{iprobe} & (EXP.CellInfo.Probe == iprobe) & EXP.CellInfo.Goodcluster);
                for icell = 1:numel(cellprobeidx)
                    disp(['Processing cell #' num2str(icell) 'out of ' num2str(numel(cellprobeidx))]);
                    EXP.data.es.spikeTrain = Spiketrain_orig;
                    EXP.data.es.spikeTrain(:,cellprobeidx(icell)) = 0;
                    EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
                        'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'smth_spd', smth_spd, 'latcorrection', latcorrection,...
                        'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                        'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin, 'ProbeID', iprobe);
                    
                    res.DecCells{iprobe,cellprobeidx(icell)} = EXP.Bayes.DecCellidx{iprobe} & (EXP.CellInfo.Probe == iprobe);
                    res.DecCells{iprobe,cellprobeidx(icell)}(cellprobeidx(icell)) = false;
                    
                    amp_th = 0;maxtol = 1;
                    res.Xpred_ave{cellprobeidx(icell)} = getCircularAverage(EXP.Bayes.Posterior0{iprobe}',amp_th,maxtol);
                    amp_th = 0;maxtol = 0.1;
                    res.Xpred_max{cellprobeidx(icell)} = getCircularAverage(EXP.Bayes.Posterior0{iprobe}',amp_th,maxtol);
                end
                catch
                    keyboard
                end
            end
            
            nbcont = numel(EXP.SubsetVal.contrast) + 1;
            nbgain = numel(EXP.SubsetVal.gain)+1;
            nbroomlength = numel(EXP.SubsetVal.roomlength);
            nboutcome = numel(EXP.SubsetVal.outcome);
            
            CellInfo.Probe = EXP.CellInfo.Probe;
            CellInfo.SSI_traj = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfo_traj = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfoPerSpike_traj = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.field_traj = cell(nbcont, nbgain, nbroomlength, nboutcome);
            
            CellInfo.SSI_avedec = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfo_avedec = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfoPerSpike_avedec = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.field_avedec = cell(nbcont, nbgain, nbroomlength, nboutcome);
            
            CellInfo.SSI_maxdec = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfo_maxdec = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfoPerSpike_maxdec = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.field_maxdec = cell(nbcont, nbgain, nbroomlength, nboutcome);
            
            CellInfo.SSI_avedecAll = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfo_avedecAll = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfoPerSpike_avedecAll = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.field_avedecAll = cell(nbcont, nbgain, nbroomlength, nboutcome);
            
            CellInfo.SSI_maxdecAll = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfo_maxdecAll = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfoPerSpike_maxdecAll = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.field_maxdecAll = cell(nbcont, nbgain, nbroomlength, nboutcome);
            
            CellInfo.SSI_avedecCross = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfo_avedecCross = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfoPerSpike_avedecCross = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.field_avedecCross= cell(nbcont, nbgain, nbroomlength, nboutcome);
            
            CellInfo.SSI_maxdecCross = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfo_maxdecCross = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.SpatialInfoPerSpike_maxdecCross = cell(nbcont, nbgain, nbroomlength, nboutcome);
            CellInfo.field_maxdecCross = cell(nbcont, nbgain, nbroomlength, nboutcome);
            
            nbXbinsmth = round(1/(Xsmthwin/100));
            for c = 1:nbcont
                for g = 1:nbgain
                    for r = 1:nbroomlength
                        for o = [1 2 3]
                            if c > numel(EXP.SubsetVal.contrast)
                                contidx = find(EXP.SubsetVal.contrast>0);
                            else
                                contidx = c;
                            end
                            if g > numel(EXP.SubsetVal.gain)
                                gainidx = find(EXP.SubsetVal.gain>0);
                            else
                                gainidx = g;
                            end
                            if o == 1
                                oidx = [1 2 4 5];%error trials
                            elseif o == 2
                                oidx = 1:5;%all trials
                            elseif o == 3
                                oidx = o;%successful trials
                            end
                            dx = 1;
                            samplerate = 1./EXP.data.es.sampleSize;
                            CellInfo.field_traj{c, g, r, o} = NaN(EXP.CellInfo.NumCells,nDecbins);
                            CellInfo.field_avedec{c, g, r, o} = NaN(EXP.CellInfo.NumCells,nDecbins);
                            CellInfo.field_maxdec{c, g, r, o} = NaN(EXP.CellInfo.NumCells,nDecbins);
                            CellInfo.field_avedecAll{c, g, r, o} = NaN(EXP.CellInfo.NumCells,nDecbins);
                            CellInfo.field_maxdecAll{c, g, r, o} = NaN(EXP.CellInfo.NumCells,nDecbins);
                            CellInfo.field_avedecCross{c, g, r, o} = NaN(EXP.CellInfo.NumCells,nDecbins);
                            CellInfo.field_maxdecCross{c, g, r, o} = NaN(EXP.CellInfo.NumCells,nDecbins);
                            CellInfo.SSI_traj{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SSI_avedec{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SSI_maxdec{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SSI_avedecAll{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SSI_maxdecAll{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SSI_avedecCross{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SSI_maxdecCross{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfo_traj{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfo_avedec{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfo_maxdec{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfo_avedecAll{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfo_maxdecAll{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfo_avedecCross{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfo_maxdecCross{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfoPerSpike_traj{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfoPerSpike_avedec{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfoPerSpike_maxdec{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfoPerSpike_avedecAll{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfoPerSpike_maxdecAll{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfoPerSpike_avedecCross{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            CellInfo.SpatialInfoPerSpike_maxdecCross{c, g, r, o} = NaN(1,EXP.CellInfo.NumCells);
                            
                            for icell = 1:EXP.CellInfo.NumCells
                                if ~isempty(res.Xpred_ave{icell})
                                    tidx = EXP.getSubsets(contidx, gainidx, r, oidx);% & abs(res.Xpred_ave{icell}-1 - EXP.data.es.trajPercent)<30;
                                    spktrain = smthInTime(EXP.data.es.spikeTrain(:,icell), mean(samplerate), Tsmthwin, 'same', [], 'boxcar_centered');%EXP.data.es.spikeTrain(:,icell);
                                    
                                    varX = EXP.data.es.trajPercent;
                                    Xtemp = unwrap(varX/max(varX)*2*pi)*max(varX)/(2*pi);
                                    Xtemp = smthInTime(Xtemp, mean(samplerate), Tsmthwin, 'same', [], 'boxcar_centered');
                                    Xtemp = mod(Xtemp,max(varX));
                                    varX = Xtemp;

                                    [SInfo, SInfoperSpk, map] = SpatialInformation(varX(tidx), spktrain(tidx), dx, samplerate(tidx),nbXbinsmth,EXP.data.es.CircularMaze);
                                    SSI = (max(map) - mean(map))/mean(map);
                                    map = map(:)';
                                    if ~isempty(map)
                                        CellInfo.SSI_traj{c, g, r, o}(icell) = SSI;
                                        CellInfo.SpatialInfo_traj{c, g, r, o}(icell) = SInfo;
                                        CellInfo.SpatialInfoPerSpike_traj{c, g, r, o}(icell) = SInfoperSpk;
                                        CellInfo.field_traj{c, g, r, o}(icell,1:min(numel(map),nDecbins)) = map(1:min(end,nDecbins));
                                    end
                                    
                                    varX = resAll.Xpred_ave{EXP.CellInfo.Probe(icell)}-1;
                                    if ~isempty(varX)
                                        [SInfo, SInfoperSpk, map] = SpatialInformation(varX(tidx & ~isnan(varX)), spktrain(tidx & ~isnan(varX)), dx, samplerate(tidx & ~isnan(varX)),nbXbinsmth,EXP.data.es.CircularMaze);
                                        SSI = (max(map) - mean(map))/mean(map);
                                        map = map(:)';
                                        if ~isempty(map)
                                            CellInfo.SSI_avedecAll{c, g, r, o}(icell) = SSI;
                                            CellInfo.SpatialInfo_avedecAll{c, g, r, o}(icell) = SInfo;
                                            CellInfo.SpatialInfoPerSpike_avedecAll{c, g, r, o}(icell) = SInfoperSpk;
                                            CellInfo.field_avedecAll{c, g, r, o}(icell,1:min(numel(map),nDecbins)) = map(1:min(end,nDecbins));
                                        end
                                    end
                                    
                                    varX = res.Xpred_ave{icell}-1;
                                    if ~isempty(varX)
                                        [SInfo, SInfoperSpk, map] = SpatialInformation(varX(tidx & ~isnan(varX)), spktrain(tidx & ~isnan(varX)), dx, samplerate(tidx & ~isnan(varX)),nbXbinsmth,EXP.data.es.CircularMaze);
                                        SSI = (max(map) - mean(map))/mean(map);
                                        map = map(:)';
                                        if ~isempty(map)
                                            CellInfo.SSI_avedec{c, g, r, o}(icell) = SSI;
                                            CellInfo.SpatialInfo_avedec{c, g, r, o}(icell) = SInfo;
                                            CellInfo.SpatialInfoPerSpike_avedec{c, g, r, o}(icell) = SInfoperSpk;
                                            CellInfo.field_avedec{c, g, r, o}(icell,1:min(numel(map),nDecbins)) = map(1:min(end,nDecbins));
                                        end
                                    end
                                    
                                    varX = res.Xpred_max{icell}-1;
                                    if ~isempty(varX)
                                        [SInfo, SInfoperSpk, map] = SpatialInformation(varX(tidx & ~isnan(varX)), spktrain(tidx & ~isnan(varX)), dx, samplerate(tidx & ~isnan(varX)),nbXbinsmth,EXP.data.es.CircularMaze);
                                        SSI = (max(map) - mean(map))/mean(map);
                                        map = map(:)';
                                        if ~isempty(map)
                                            CellInfo.SSI_maxdec{c, g, r, o}(icell) = SSI;
                                            CellInfo.SpatialInfo_maxdec{c, g, r, o}(icell) = SInfo;
                                            CellInfo.SpatialInfoPerSpike_maxdec{c, g, r, o}(icell) = SInfoperSpk;
                                            CellInfo.field_maxdec{c, g, r, o}(icell,1:min(numel(map),nDecbins)) = map(1:min(end,nDecbins));
                                        end
                                    end
                                    
                                    varX = resAll.Xpred_max{EXP.CellInfo.Probe(icell)}-1;
                                    if ~isempty(varX)
                                        [SInfo, SInfoperSpk, map] = SpatialInformation(varX(tidx & ~isnan(varX)), spktrain(tidx & ~isnan(varX)), dx, samplerate(tidx & ~isnan(varX)),nbXbinsmth,EXP.data.es.CircularMaze);
                                        SSI = (max(map) - mean(map))/mean(map);
                                        map = map(:)';
                                        if ~isempty(map)
                                            CellInfo.SSI_maxdecAll{c, g, r, o}(icell) = SSI;
                                            CellInfo.SpatialInfo_maxdecAll{c, g, r, o}(icell) = SInfo;
                                            CellInfo.SpatialInfoPerSpike_maxdecAll{c, g, r, o}(icell) = SInfoperSpk;
                                            CellInfo.field_maxdecAll{c, g, r, o}(icell,1:min(numel(map),nDecbins)) = map(1:min(end,nDecbins));
                                        end
                                    end
                                    
                                    if nbProbe > 1
                                        if EXP.CellInfo.Probe(icell) > 1
                                            iprobe = 1;
                                        else
                                            iprobe = 2;
                                        end
                                        varX = resAll.Xpred_ave{iprobe}-1;
                                        if ~isempty(varX)
                                            [SInfo, SInfoperSpk, map] = SpatialInformation(varX(tidx & ~isnan(varX)), spktrain(tidx & ~isnan(varX)), dx, samplerate(tidx & ~isnan(varX)),nbXbinsmth,EXP.data.es.CircularMaze);
                                            SSI = (max(map) - mean(map))/mean(map);
                                            map = map(:)';
                                            if ~isempty(map)
                                                CellInfo.SSI_avedecCross{c, g, r, o}(icell) = SSI;
                                                CellInfo.SpatialInfo_avedecCross{c, g, r, o}(icell) = SInfo;
                                                CellInfo.SpatialInfoPerSpike_avedecCross{c, g, r, o}(icell) = SInfoperSpk;
                                                CellInfo.field_avedecCross{c, g, r, o}(icell,1:min(numel(map),nDecbins)) = map(1:min(end,nDecbins));
                                            end
                                        end
                                        
                                        varX = resAll.Xpred_max{iprobe}-1;
                                        if ~isempty(varX)
                                            [SInfo, SInfoperSpk, map] = SpatialInformation(varX(tidx & ~isnan(varX)), spktrain(tidx & ~isnan(varX)), dx, samplerate(tidx & ~isnan(varX)),nbXbinsmth,EXP.data.es.CircularMaze);
                                            SSI = (max(map) - mean(map))/mean(map);
                                            map = map(:)';
                                            if ~isempty(map)
                                                CellInfo.SSI_maxdecCross{c, g, r, o}(icell) = SSI;
                                                CellInfo.SpatialInfo_maxdecCross{c, g, r, o}(icell) = SInfo;
                                                CellInfo.SpatialInfoPerSpike_maxdecCross{c, g, r, o}(icell) = SInfoperSpk;
                                                CellInfo.field_maxdecCross{c, g, r, o}(icell,1:min(numel(map),nDecbins)) = map(1:min(end,nDecbins));
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            save([dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_EXP '_DecodingPopCell.mat'], 'resAll', 'res', 'CellInfo','-v7.3');
        end
    end
end