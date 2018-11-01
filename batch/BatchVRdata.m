function EXP = BatchVRdata
%Process VR data in batch, as defined in the getExperimentList2p function
%and save the processed data in the data directory
%Then, you can make figures from the saved processed data using the
%following functions in the command window:
% - for decoding maps,
% res = BatchBayesGrandAverage;
% popres = PopBayesAnalysis(res,[]);
% PopPlot(popres,[],[],'DistriXdecAve');
% - for single cell response profiles,
% cellprop = BatchCellInfo;
% PopPlotcellprop(cellprop);
% - for latencies and lick distributions
% coming soon

SetDirs;

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
    area_str = 'CA1';
end
prompt = {'Speed Threshold';'Window size (ms)';'Spatial smth (%)';'# of X bins';'FoverwriteCellinfo'};
dlg_title = 'Parameters';
num_lines = 1;
def = {'1';'150';'4';'100';'1'};
nruns = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(nruns)
    SpeedThreshold = str2num(nruns{1});
    Tsmthwin = str2num(nruns{2});
    Xsmthwin = str2num(nruns{3});
    nDecbins = str2num(nruns{4});
    FoverwriteCellinfo = logical(str2num(nruns{5}));
end

Nperm_cellprop = 2;%100
Nthetaphsbins = 0;


%batch files
for ianimal = 1:numel(expt)
    for iseries = 1:numel(expt(ianimal).series)
        %create the main data object
        EXP = TVRData;
        EXP.animal = expt(ianimal).animal;
        EXP.series = num2str(expt(ianimal).series{iseries});
        EXP.iseries = expt(ianimal).series{iseries};
        EXP.exptList = expt(ianimal).exp{iseries};
        if strcmp(EXP.animal,'M160114C_BALL') && EXP.iseries == 323
            shanknum = [0:15 0];
            suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
        else
            shanknum = [0:7 0];
            suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
        end
        EXP.LoadVRData(shanknum, suffix, SpeedThreshold, Nthetaphsbins);
        if ~EXP.data.twophoton
            dDIRname = [DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries)];
        else
            dDIRname = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries)];
        end
        
        savedcellinfo = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_cellProperties.mat'];
        savedmaps1d = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_maps1d.mat'];
        savedVS = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_VS.mat'];
        
        if exist(savedcellinfo) && exist(savedmaps1d) && exist(savedVS) && ~FoverwriteCellinfo
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
            EXP.Calculate1Dmaps('trajPercent',Tsmthwin,Xsmthwin,delay);
            maps1d = EXP.maps1d;
            save(savedmaps1d,'maps1d');
            EXP.defineCellInfo(EXP.data.es.spikeIDs, EXP.data.es.chanIDs,EXP.data.es.ProbeIDs);
            EXP.defineCellProp(Nperm_cellprop);
            CellInfo = EXP.CellInfo;
            save(savedcellinfo,'CellInfo');
        end
        
        %             delay = 800;
        %             EXP.SimulPlaceFields(delay);
        %             EXP.Calculate1Dmaps('trajPercent',PlotVar1DMaps.SmthTimeWindow,delay);
        %             EXP.defineCellProp;
        
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
        nthetabins = 1;%1;%6;%
        nspdbins = 5;%1;%
        smth_spd = Tsmth_win;
        speed_th = SpeedThreshold;
        
        kfold = 20;
        FoptiSmooth = 0;
        FGoodcluster = 1;
        FUnsortedcluster = 0;
        FMUAcluster = 0;
        maxRate = inf;
        zth = -inf;
        SSImin = -inf;
        cellselstr = 'goodonly';%'goodandMUA';%'allbutINs';%'goodonly';
        
        %             [nspdbins, smth_spd] = getOptSpdParams(EXP);
        %             latcorrection = getOptLatParams(EXP);
        
        latcorrection = 0;
        alpha = 0;
        delta = 0;
        
        EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
            'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'smth_spd', smth_spd, 'latcorrection', latcorrection,...
            'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
            'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin);
        
        maxtol = 0.1;
        EXP = BayesDecoderAverages(EXP,maxtol);%ComputeBayesAverage(EXP,Nperm_bayes);
        save([dDIRname filesep 'EXP_win' num2str(Tsmth_win) '_' 'speedThresh' num2str(SpeedThreshold) '_' num2str(nspdbins) 'speedbins_' cellselstr '.mat'],'EXP','-v7.3');
    end
end
end