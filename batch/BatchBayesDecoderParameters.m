function EXP = BatchBayesDecoderParameters(EXP, DIRS)
model = 'Data';%'Vision';%'Data';%'Local Distance';%'Global Distance';%
if strcmp(model,'Data')
    latencylist = [];
else
    latencylist = 300;%[0 200 400 600 800 1000 1200];
end
Nperm_cellprop = 2;
ParamsMeanErr = cell(1,max(1,numel(latencylist)));
ParamsCorr = cell(1,max(1,numel(latencylist)));

batch2p = false;
if batch2p
    expt = getExperimentList2p;
else
    expt = getExperimentList;
end

Tsmthwin = 250;
Xsmthwin = 4;
SpeedThreshold = 5;
nspeedbins = 3;
neyeXbins = 3;
nthetaphsbins = 0;
cellstr = 'goodonly';
filesuffix_EXP = ['Twin' num2str(Tsmthwin) '_' 'Xwin' num2str(Xsmthwin) '_' 'spdth' num2str(SpeedThreshold) '_' num2str(nspeedbins) 'speedbins' '_' num2str(neyeXbins) 'eyebins' '_' num2str(nthetaphsbins) 'thetabins' '_' cellstr];
            
nanimals = numel(expt);
Nperm_cellprop = 2;%100
filesuffix = 'win150_5speedbins_goodonly';
for ianimal = 1:numel(expt)
    EXP.animal = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
%         try
            EXP.series = num2str(expt(ianimal).series{iseries});
            EXP.iseries = expt(ianimal).series{iseries};
            EXP.exptList = expt(ianimal).exp{iseries};
            
            dDIRname = ['D:\DATA\batch'  filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];%[DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries)];
            savedfile = [dDIRname filesep 'EXP_' filesuffix_EXP '.mat'];
            
%             savedfile = [DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'EXP_' filesuffix '.mat'];
%             savedfil2p = [DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep 'EXP_' filesuffix '.mat'];
            disp([EXP.animal ' series ' num2str(EXP.iseries)]);
            if ~exist(savedfile,'file')
                if strcmp(EXP.animal,'M160114C_BALL') && EXP.iseries == 323
                    shanknum = [0:15 0];
                    suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
                else
                    shanknum = [0:7 0];
                    suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
                end
                EXP.LoadVRData(shanknum, suffix);
                
                
                savedcellinfo = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_cellProperties.mat'];
                savedmaps1d = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d.mat'];
                savedVS = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_VS.mat'];
                
                FoverwriteCellinfo = false;
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
                    try
                        EXP.CalculateStimTuning([], shanknum, suffix);
                    catch
                        warning('something''s wrong with visual stimuli')
                        EXP.data.VStuning = [];
                        EXP.data.VSstim = [];
                    end
                    delay = 0;
                    Tsmthwin = 150;
                    Xsmthwin = 4;
                    EXP.Calculate1Dmaps('trajPercent',Tsmthwin,Xsmthwin,delay);
                    EXP.defineCellInfo(EXP.data.es.spikeIDs, EXP.data.es.chanIDs,EXP.data.es.ProbeIDs);
                    EXP.defineCellProp(Nperm_cellprop);
                    CellInfo = EXP.CellInfo;
                    maps1d = EXP.maps1d;
                    VStuning = EXP.data.VStuning;
                    VSstim = EXP.data.VSstim;
                    save(savedcellinfo,'CellInfo','VStuning','VSstim','maps1d');
                end
            else
                S = load(savedfile);
                EXP.Copyobj(S.EXP);
            end
        
        for iSim = 1:max(1,numel(latencylist))
            if ~strcmp(model,'Data')
                latency = latencylist(iSim);
                EXP.SimulPlaceFields(latency,model,3);
            end
            
            latparams.nspdbinslist = [3];%1;%[1 3 5 10 15];[nspdbins, smth_spd] = getOptSpdParams(EXP);            
            latparams.smth_spdlist = [250];%[500];
            latparams.latcorrectionlist = [0:50:1500];%[0:50:1500];%getOptLatParams(EXP);latcorrectionlist = latcorrectionlist(:);%[0 100 200 300 400 500 600 700 800 900 1000];%
            latparams.alphalist = [0];
            latparams.deltalist = [0];
            
            [xParams, latParamsMeanErr0, latParamsMeanErr1, latParamsMeanErr2,...
             latParamsCorr, latParamsErr1, latParamsErr2, latParamsLickdistri,...
             latParamsLickdistriXpred, latParamsPostbeforelick,...
             latParamsMeanPostbeforelick] = TestDecodingParams(EXP, latparams.nspdbinslist, latparams.smth_spdlist, latparams.latcorrectionlist, latparams.alphalist, latparams.deltalist);            
            
            latparams.xParams = xParams;
            latEstim.ParamsMeanErr0 = latParamsMeanErr0;
            latEstim.ParamsMeanErr1 = latParamsMeanErr1;
            latEstim.ParamsMeanErr2 = latParamsMeanErr2;
            latEstim.ParamsErr1 = latParamsErr1;
            latEstim.ParamsErr2 = latParamsErr2;
            latEstim.ParamsLickdistri = latParamsLickdistri;
            latEstim.ParamsLickdistriXpred = latParamsLickdistriXpred;
            latEstim.ParamsPostbeforelick = latParamsPostbeforelick;
            latEstim.ParamsMeanPostbeforelick = latParamsMeanPostbeforelick;
            
            for iprobe = 1:size(latEstim.ParamsMeanErr2,1)
                for g = 1:numel(EXP.SubsetVal.gain)
                    [~,I] = max(squeeze(latEstim.ParamsMeanErr0(iprobe,:,:,:,:,:,g,:)),[],2);
                    latEstim.meanPostXVec{iprobe,g} = squeeze(I);
                end
                latEstim.OptLatencyVec{iprobe} = max(squeeze(latEstim.ParamsMeanErr2(iprobe,:,:,:,1,:,EXP.SubsetVal.gain == mode(EXP.data.es.gain),:)),[],2);%smooth(max(squeeze(latEstim.ParamsMeanErr2(iprobe,:,:,:,1,:,EXP.SubsetVal.gain == mode(EXP.data.es.gain),:)),[],2),3);
                
                
                latEstim.MaxLickdistri{iprobe} = zeros(numel(EXP.SubsetVal.gain),1);
                latEstim.MaxLickdistriXpred{iprobe} = zeros(numel(EXP.SubsetVal.gain),1);
                latEstim.MaxPostbeforeLickdistri{iprobe} = zeros(numel(EXP.SubsetVal.gain),size(latEstim.ParamsPostbeforelick,8));
                latEstim.MaxMeanPostbeforeLickdistri{iprobe} = zeros(numel(EXP.SubsetVal.gain),size(latEstim.ParamsPostbeforelick,8));
                for g = 1:numel(EXP.SubsetVal.gain)
                    [~, idxmax] = max(squeeze(latEstim.ParamsLickdistri(iprobe,latparams.latcorrectionlist==0,:,:,1,:,g,:)));
                    if ~isempty(idxmax)
                        latEstim.MaxLickdistri{iprobe}(g) = xParams(idxmax);
                    else
                        latEstim.MaxLickdistri{iprobe}(g) = NaN;
                    end
                    [~, idxmax] = max(squeeze(latEstim.ParamsLickdistri(iprobe,latparams.latcorrectionlist==0,:,:,1,:,g,:)));
                    if ~isempty(idxmax)
                        latEstim.MaxLickdistriXpred{iprobe}(g) = xParams(idxmax);
                    else
                        latEstim.MaxLickdistriXpred{iprobe}(g) = NaN;
                    end
                    for tt = 1:size(latEstim.ParamsPostbeforelick,8)
                        [~, idxmax] = max(squeeze(latEstim.ParamsPostbeforelick(iprobe,latparams.latcorrectionlist==0,:,:,1,:,g,tt,:)));
                        if ~isempty(idxmax)
                            latEstim.MaxPostbeforeLickdistri{iprobe}(g,tt) = xParams(idxmax);
                        else
                            latEstim.MaxPostbeforeLickdistri{iprobe}(g,tt) = NaN;
                        end
                        [~, idxmax] = max(squeeze(latEstim.ParamsMeanPostbeforelick(iprobe,latparams.latcorrectionlist==0,:,:,1,:,g,tt,:)));
                        if ~isempty(idxmax)
                            latEstim.MaxMeanPostbeforeLickdistri{iprobe}(g,tt) = xParams(idxmax);
                        else
                            latEstim.MaxMeanPostbeforeLickdistri{iprobe}(g,tt) = NaN;
                        end
                    end
                end
                
                [~, idxopt] = max(latEstim.OptLatencyVec{iprobe});
                latEstim.OptLatency(iprobe) = latparams.latcorrectionlist(idxopt);
                latEstim.OptLatency(iprobe) = max(0,latEstim.OptLatency(iprobe));
            end
            
%             kalparams.nspdbinslist = [5];%1;%[1 3 5 10 15];[nspdbins, smth_spd] = getOptSpdParams(EXP);            
%             kalparams.smth_spdlist = [0];%[500];
%             OptLatency = 0;
% %             if sum(OptLatency==0)==0
% %                 OptLatency = [OptLatency 0];
% %             end
%             kalparams.latcorrectionlist = [OptLatency];%getOptLatParams(EXP);latcorrectionlist = latcorrectionlist(:);%[0 100 200 300 400 500 600 700 800 900 1000];%
%             %latcorrectionlist = [0 100 200 300 400 500 600 700 800 900 1000];%getOptLatParams(EXP);
%             kalparams.alphalist = [0:0.1:2];
%             kalparams.deltalist = [0];
%             
%             [xParams, kalParamsMeanErr0, kalParamsMeanErr1, kalParamsMeanErr2, kalParamsCorr, kalParamsErr1, kalParamsErr2, kalParamsLickdistri, kalParamsLickdistriXpred, kalParamsPostbeforelick, kalParamsMeanPostbeforelick] = TestDecodingParams(EXP, kalparams.nspdbinslist, kalparams.smth_spdlist, kalparams.latcorrectionlist, kalparams.alphalist, kalparams.deltalist);
%             
%             kalparams.xParams = xParams;
%             kalEstim.ParamsMeanErr0 = kalParamsMeanErr0;
%             kalEstim.ParamsMeanErr1 = kalParamsMeanErr1;
%             kalEstim.ParamsMeanErr2 = kalParamsMeanErr2;
%             kalEstim.ParamsErr1 = kalParamsErr1;
%             kalEstim.ParamsErr2 = kalParamsErr2;
%             kalEstim.ParamsLickdistri = kalParamsLickdistri;
%             kalEstim.ParamsLickdistriXpred = kalParamsLickdistriXpred;
%             kalEstim.ParamsPostbeforelick = kalParamsPostbeforelick;
%             kalEstim.ParamsMeanPostbeforelick = kalParamsMeanPostbeforelick;
%             
%             kalEstim.OptVisstdVec = cell(size(kalEstim.ParamsMeanErr2,1),numel(EXP.SubsetVal.gain));
%             kalEstim.MaxLickdistri = zeros(size(kalEstim.ParamsMeanErr2,1),numel(EXP.SubsetVal.gain));
%             kalEstim.MaxLickdistriXpred = zeros(size(kalEstim.ParamsMeanErr2,1),numel(EXP.SubsetVal.gain));
%             kalEstim.MaxPostbeforeLickdistri = zeros(size(kalEstim.ParamsMeanErr2,1),numel(EXP.SubsetVal.gain),size(kalEstim.ParamsPostbeforelick,8));
%             kalEstim.MaxMeanPostbeforeLickdistri = zeros(size(kalEstim.ParamsMeanErr2,1),numel(EXP.SubsetVal.gain),size(kalEstim.ParamsPostbeforelick,8));
%             for iprobe = 1:size(kalEstim.ParamsMeanErr2,1)
%                 for g = 1:numel(EXP.SubsetVal.gain)
%                     [~,I] = max(squeeze(kalParamsMeanErr0(iprobe,1,:,:,:,:,g,:)),[],2);
%                     kalEstim.meanPostXVec{iprobe,g} = squeeze(I);
%                 end
%                 
%                 for g = 1:numel(EXP.SubsetVal.gain)
%                     [~, optKidx] = min(abs(kalEstim.meanPostXVec{iprobe,g}-kalEstim.meanPostXVec{iprobe,EXP.SubsetVal.gain == mode(EXP.data.es.gain)}));
%                     kalEstim.OptVisstd(iprobe,g) = kalparams.alphalist(optKidx);
%                     
%                     [~, idxmax] = max(squeeze(kalEstim.ParamsLickdistri(iprobe,1,:,:,optKidx,:,g,:)));
%                     kalEstim.MaxLickdistri(iprobe,g) = xParams(idxmax);
%                     [~, idxmax] = max(squeeze(kalEstim.ParamsLickdistriXpred(iprobe,1,:,:,optKidx,:,g,:)));
%                     kalEstim.MaxLickdistriXpred(iprobe,g) = xParams(idxmax);
%                     for tt = 1:size(kalEstim.ParamsPostbeforelick,8)
%                         [~, idxmax] = max(squeeze(kalEstim.ParamsPostbeforelick(iprobe,1,:,:,optKidx,:,g,tt,:)));
%                         kalEstim.MaxPostbeforeLickdistri(iprobe,g,tt) = xParams(idxmax);
%                         [~, idxmax] = max(squeeze(kalEstim.ParamsMeanPostbeforelick(iprobe,1,:,:,optKidx,:,g,tt,:)));
%                         kalEstim.MaxMeanPostbeforeLickdistri(iprobe,g,tt) = xParams(idxmax);
%                     end
%                 end
%             end
           
%             [MeanErr] = TestDecodingParams2(EXP, nspdbinslist, smth_spdlist, latcorrectionlist);
            
%             ParamsMeanErr{iSim} = MeanErr;
%             ParamsCorr{iSim} = Corr;
        end
        if strcmp(model,'Data')
            paramsfilename = 'decoderParams';%'decoderEyeParams';%'decoderLatParams';%'decoderVisSpdParams';
        else
            paramsfilename = ['decoderLatParams_' model];%['decoderSpdParams_' model];
        end
        
%         assignin('base','latParamsMeanErr0',latParamsMeanErr0)
%         assignin('base','latParamsMeanErr1',latParamsMeanErr1)
%         assignin('base','latParamsMeanErr2',latParamsMeanErr2)
%         assignin('base','latParamsErr1',latParamsErr1)
%         assignin('base','latParamsErr2',latParamsErr2)
        
%         assignin('base','kalParamsMeanErr0',kalParamsMeanErr0)
%         assignin('base','kalParamsMeanErr1',kalParamsMeanErr1)
%         assignin('base','kalParamsMeanErr2',kalParamsMeanErr2)
%         assignin('base','kalParamsErr1',kalParamsErr1)
%         assignin('base','kalParamsErr2',kalParamsErr2)
%         save([DIRS.multichanspikes filesep EXP.animal filesep EXP.series filesep paramsfilename '_kalman.mat'], 'kalparams', 'kalEstim','-v7.3');
        save([dDIRname filesep  EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_EXP '_' paramsfilename '_latency.mat'], 'latparams', 'latEstim','-v7.3');
    end
end