function [EXP, PlotVarDecoder, DecodingDialog] = BatchBayesDecoder(EXP, DIRS, PlotVarDecoder, DecodingDialog)
strlistvarname = {'2p data','electrophys data'};
[varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'dataset', 'SelectionMode', 'single', 'InitialValue', 1);
if ok && varnamesel == 1
    batch2p = true;
    expt = getExperimentList2p;
    strlistvarname = {'V1medial','V1lateral','PPC', 'AL', 'V1medialV1lateral', 'V1medialV1lateralPPCAL'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'area', 'SelectionMode', 'single', 'InitialValue', 1);
    area_str = strlistvarname{varnamesel};
elseif ok
    batch2p = false;
    expt = getExperimentList;
    area_str = 'CA1V1';
end
prompt = {'Speed Threshold';'nthetaphsbins';'nspdbins';'neyebins';'nphsbins';'Tsmthwin_field (ms)';'Spatial smth (%)';'Tsmthwin_dec (ms)';'# of X bins';'FGoodcluster';'FMUAcluster';'FoverwriteCellinfo'};
dlg_title = 'Parameters';
num_lines = 1;
def = {'5';'0';'1';'1';'1';'15';'4';'50';'100';'1';'0';'1'};
nruns = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(nruns)
    SpeedThreshold = str2num(nruns{1});
    nthetaphsbins = str2num(nruns{2});
    nspeedbins = str2num(nruns{3});
    neyeXbins = str2num(nruns{4});
    nphsbins = str2num(nruns{5});
    Tsmthwin = str2num(nruns{6});
    Xsmthwin = str2num(nruns{7});
    Tsmthwin_dec = str2num(nruns{8});
    nDecbins = str2num(nruns{9});
    FGoodcluster = logical(str2num(nruns{10}));
    FMUAcluster = logical(str2num(nruns{11}));
    FoverwriteCellinfo = logical(str2num(nruns{12}));
end

nanimals = numel(expt);
Nperm_cellprop = 0;%100

FUnsortedcluster = 0;
maxRate = inf;
zth = -inf;
SSImin = -inf;

if FGoodcluster && ~FMUAcluster
    cellstr = 'goodonly';%'goodonly_unwrapped';%'goodonly_spktimes';%
elseif FGoodcluster && FMUAcluster
    cellstr = 'All';
end

filesuffix_EXP = ['Twin' num2str(Tsmthwin) '_' 'Xwin' num2str(Xsmthwin) '_' 'spdth' num2str(SpeedThreshold) '_' 'Decwin' num2str(Tsmthwin_dec) '_' num2str(nspeedbins) 'speedbins' '_' num2str(neyeXbins) 'eyebins' '_' num2str(nphsbins) 'thetabins' '_' cellstr '_spdquantilesRand'];
filesuffix_cellprop = ['Twin' num2str(Tsmthwin) '_' 'Xwin' num2str(Xsmthwin) '_' 'spdth' num2str(SpeedThreshold) '_' num2str(nthetaphsbins) 'thetabins' '_' cellstr];

for ianimal = 1:numel(expt)
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
                dDIRname = ['D:\DATA\batch\2p'  filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];%[DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];
            end
            if ~isdir(dDIRname)
                mkdir(dDIRname)
            end
            
            savedfile = [dDIRname filesep 'EXP_' filesuffix_EXP '.mat'];
            disp([EXP.animal ' series ' num2str(EXP.iseries)]);
            
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
            savedmaps1d_2fold = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d_2fold.mat'];
            savedmaps2d_spd = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2d_spd.mat'];
            savedmaps1d_phs = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d_phs.mat'];
            savedmaps1d_phs_2fold = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d_phs_2fold.mat'];
            savedmaps2d_phs = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2d_phs.mat'];
            savedmaps2d_phs_2fold = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2d_phs_2fold.mat'];
            savedVS = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_VS.mat'];
            
            if exist(savedcellinfo) && ~FoverwriteCellinfo
                disp('loading cellInfo from saved file');
                S = load(savedcellinfo);
                EXP.CellInfo = S.CellInfo;
                S = load(savedVS);
                EXP.data.VStuning = S.VStuning;
                EXP.data.VSstim = S.VSstim;
                S = load(savedmaps1d);
                EXP.maps1d.trajPercent = S.maps1d;
                S = load(savedmaps1d_2fold);
                EXP.maps1d.trajPercent_2fold = S.maps1d;
%                 S = load(savedmaps2d_spd);
%                 EXP.maps2d.trajPercent_smthBallSpd = S.maps2d;
                S = load(savedmaps1d_phs);
                EXP.maps2d.LFPphase2 = S.maps1d;
                S = load(savedmaps1d_phs_2fold);
                EXP.maps2d.LFPphase2_2fold = S.maps1d;
                S = load(savedmaps2d_phs);
                EXP.maps2d.trajPercent_LFPphase2 = S.maps2d;
                S = load(savedmaps2d_phs_2fold);
                EXP.maps2d.trajPercent_LFPphase2_2fold = S.maps2d;
                
                EXP.defineCellInfo(EXP.data.es.spikeIDs, EXP.data.es.chanIDs,EXP.data.es.ProbeIDs);
                %we call definecellprop here to get the correction factor for the LFPphase
                if isfield(EXP.data.es,'LFPphase')
                    EXP.defineCellProp(Nperm_cellprop);
                    nProbe = numel(unique(EXP.CellInfo.Probe));
                    EXP.data.es.LFPphase2 = NaN(size(EXP.data.es.LFPphase,1),nProbe);
                    for iprobe = 1:nProbe
                        EXP.data.es.LFPphase2(:,iprobe) = mod(EXP.data.es.LFPphase - round(EXP.CellInfo.LFP2Spike_phscorrMUAnorm{1}),360);%- EXP.CellInfo.LFP2Spike_phscorrMUA{iprobe},360);
                    end
                end
            else
                EXP.maps1d = [];
                EXP.maps2d = [];
                EXP.defineCellInfo(EXP.data.es.spikeIDs, EXP.data.es.chanIDs,EXP.data.es.ProbeIDs);
                %we call definecellprop here to get the correction factor for the LFPphase
                if isfield(EXP.data.es,'LFPphase')
                    EXP.defineCellProp(Nperm_cellprop);
                    nProbe = numel(unique(EXP.CellInfo.Probe));
                    EXP.data.es.LFPphase2 = NaN(size(EXP.data.es.LFPphase,1),nProbe);
                    for iprobe = 1:nProbe
                        EXP.data.es.LFPphase2(:,iprobe) = mod(EXP.data.es.LFPphase - round(EXP.CellInfo.LFP2Spike_phscorrMUAnorm{1}),360);%- EXP.CellInfo.LFP2Spike_phscorrMUA{iprobe},360);
                    end
                end
                EXP.CalculateStimTuning([], shanknum, suffix);
                VStuning = EXP.data.VStuning;
                VSstim = EXP.data.VSstim;
                save(savedVS,'VStuning','VSstim');
                delay = 0;
                Xbinsize = 1;
                tic
                EXP.Calculate1Dmaps('trajPercent', Tsmthwin, Xbinsize, Xsmthwin,delay, EXP.data.es.CircularMaze);
                toc
%                 Spdbinsize = 0.1;
%                 SmthSpdWindow = 0.1;
%                 EXP.Calculate2Dmaps('trajPercent', 'smthBallSpd', Tsmthwin, Xbinsize, Spdbinsize, Xsmthwin, SmthSpdWindow, delay, EXP.data.es.CircularMaze, false);
                if isfield(EXP.data.es,'LFPphase2')
                    Phsbinsize = 20;
                    SmthPhsWindow = 40;
                    tic
                    EXP.Calculate1Dmaps('LFPphase2', Tsmthwin, Phsbinsize, SmthPhsWindow, delay, true);
                    toc
                    tic
                    EXP.Calculate2Dmaps('trajPercent', 'LFPphase2', Tsmthwin, Xbinsize, Phsbinsize, Xsmthwin, SmthPhsWindow, delay, EXP.data.es.CircularMaze, true);
                    toc
                end
                maps1d = EXP.maps1d.trajPercent;
                save(savedmaps1d,'maps1d');
                maps1d = EXP.maps1d.trajPercent_2fold;
                save(savedmaps1d_2fold,'maps1d');
%                 maps2d = EXP.maps2d.trajPercent_smthBallSpd;
%                 save(savedmaps2d_spd,'maps2d');
                maps1d = EXP.maps1d.LFPphase2;
                save(savedmaps1d_phs,'maps1d');
                maps1d = EXP.maps1d.LFPphase2_2fold;
                save(savedmaps1d_phs_2fold,'maps1d');
                maps2d = EXP.maps2d.trajPercent_LFPphase2;
                save(savedmaps2d_phs,'maps2d');
                maps2d = EXP.maps2d.trajPercent_LFPphase2_2fold;
                save(savedmaps2d_phs_2fold,'maps2d');
                EXP.defineCellProp(Nperm_cellprop);
                CellInfo = EXP.CellInfo;
                save(savedcellinfo,'CellInfo');
            end
            
%             EXP.defineCellProp(Nperm_cellprop);
%             CellInfo = EXP.CellInfo;
%             save(savedcellinfo,'CellInfo');
            
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
            Tsmth_win = Tsmthwin_dec;%20;%
            Xsmth_win = Xsmthwin;
            numbins = nDecbins;
            thetachannel = 34;
            nthetabins = nphsbins;%1;%6;%
            nspdbins = nspeedbins;%1;%
            neyebins = neyeXbins;%1;%
            Tsmth_field = Tsmthwin;
            speed_th = SpeedThreshold;
            
            kfold = 20;
            FoptiSmooth = 0;            
            %             [nspdbins, smth_spd] = getOptSpdParams(EXP);
            %             latcorrection = getOptLatParams(EXP);
            
            latcorrection = 0;
            alpha = 0;
            delta = 0;
            
            EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
                'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'neyebins', neyebins, 'Tsmth_field', Tsmth_field, 'latcorrection', latcorrection,...
                'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin);
            
%             numbinsX = nDecbins;
%             EXP.RunBayesDecoder2D('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
%                 'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBinsX', numbinsX, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'neyebins', neyebins, 'smth_spd', smth_spd, 'latcorrection', latcorrection,...
%                 'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
%                 'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin);
            
%             maxtol = 0.1;
%             EXP = BayesDecoderAverages(EXP,maxtol);%ComputeBayesAverage(EXP,Nperm_bayes);
            
            DecodingDialog.UpdateUIcontrol('Plots','String', PlotVarDecoder.PlotObjList, 'max', numel(PlotVarDecoder.PlotObjList), 'min', 0, 'Val', 1);
            PlotVarDecoder.ChosenObj = PlotVarDecoder.PlotObjList(1);
            DecodingDialog.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', 1:numel(EXP.SubsetVal.contrast));
            contvalidx = find(EXP.SubsetVal.contrast > 0);
            PlotVarDecoder.ChosenContrast = contvalidx(:);
            DecodingDialog.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', 1:numel(EXP.SubsetVal.gain));
            gainvalidx = find(EXP.SubsetVal.gain > 0 & EXP.SubsetVal.gain < 1);
            PlotVarDecoder.ChosenGain = gainvalidx(:)';
            DecodingDialog.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', 1);
            PlotVarDecoder.ChosenRoomlength = (1:EXP.SubsetVal.roomlength);
            DecodingDialog.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', find(EXP.SubsetVal.outcome == 2));
            PlotVarDecoder.ChosenOutcome = find(EXP.SubsetVal.outcome == 2);
            
            save(savedfile,'EXP','-v7.3');
        end
    end
end
end