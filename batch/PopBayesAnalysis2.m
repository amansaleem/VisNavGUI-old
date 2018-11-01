function  popres = PopBayesAnalysis2(res,batch2p,average_type,FshowsessionFig)
if nargin < 2
    batch2p = [];
end
if nargin < 3
    average_type = 'sessions';%'trials'
end
if nargin < 4
    FshowsessionFig = false;
end
cl{1} = 'c';
cl{2} = 'k';
cl{3} = 'm';

popres.average_type = average_type;

resPosName = 'PostXSum';
resErrName = 'ErrXSum';
resBiasName = 'BiasXSum';

popres.DecAveType = res.DecAveType;
popres.s_FsmoothOutput = false;
popres.FsmoothOutput = true;

% resPosName = 'MaxDecXSum';
% resErrName = 'ErrDecXSum';

% resPosName = 'MeanDecXSum';
% resErrName = 'BiasDecXSum';

if FshowsessionFig
    figdirname = 'C:\Users\experiment\Desktop\figure_temp';
    prompt = {'Save single session fig?';'figure directory'};
    dlg_title = 'Parameters';
    num_lines = 1;
    def = {'0';figdirname};
    nruns = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(nruns)
        FshowsessionFig = logical(str2num(nruns{1}));
        figdirname = nruns{2};
        if ~isdir(figdirname)
            mkdir(figdirname)
        end
    end
end

popres.lambdaSmooth = 1;
if isempty(batch2p)
    strlistvarname = {'2p data','electrophys data'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'dataset', 'SelectionMode', 'single', 'InitialValue', 1);
    if ok && varnamesel == 1
        batch2p = true;
    elseif ok
        batch2p = false;
    end
end
if batch2p
    expt = getExperimentList2p;
    strlistvarname = {'V1medial','V1lateral','PPC', 'AL'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'area', 'InitialValue', 1);
    area_str = cell2mat(strlistvarname(varnamesel));
    nProbe = 1;
    popres.probestr{1} = area_str;
    popres.probestr{2} = 'none';
elseif ~batch2p
    expt = getExperimentList;
    area_str = 'CA1V1';
    nProbe = 2;
    popres.probestr{1} = 'CA1';
    popres.probestr{2} = 'V1';
end

if batch2p
    popres.maxtol = 1;%0.1;%
    popres.amp_th = 0;%1;%
    popres.interpol = 1;%0.05;%
    popres.maxtolpred = 1;%0.1;%
    popres.amp_thpred = 0;%1;%
    popres.interpolpred = 1;%0.05;%
    popres.maxtolcorr = 0.1;%1;
    popres.amp_thcorr = 0;%0;
    popres.interpolcorr = 0.05;%1;
else
    popres.maxtol = 1;%0.1;%
    popres.amp_th = 0;%1;%
    popres.interpol = 1;%0.05;%
    popres.maxtolpred = 1;%0.1;%
    popres.amp_thpred = 0;%1;%
    popres.interpolpred = 1;%0.05;%
    popres.maxtolcorr = 0.1;%1;
    popres.amp_thcorr = 1;%0;
    popres.interpolcorr = 0.05;%1;
end

popres.Tsmthwin = res.Tsmthwin;
popres.Xsmthwin = 2;%res.Xsmthwin;
popres.SpeedThreshold = res.SpeedThreshold;
popres.nSpdbins = res.nSpdbins;
popres.nVisbins = res.nVisbins;
popres.nEyebins = res.nEyebins;
popres.nPhsbins = res.nPhsbins;

popres.dx = 1;%0.1;

nanimal = numel(expt);
ngain = 3;
popres.PostXAll = cell(nProbe,ngain);
popres.ErrXAll = cell(nProbe,ngain);
popres.BiasXAll = cell(nProbe,ngain);
popres.Xsum = cell(nProbe,ngain);
popres.X = cell(nProbe,ngain);
popres.runspeed = cell(nProbe,ngain);
popres.visspeed = cell(nProbe,ngain);
popres.eyeXpos = cell(nProbe,ngain);
popres.licks = cell(nProbe,ngain);

popres.s_PostXAll = cell(nProbe,ngain);
popres.s_ErrXAll = cell(nProbe,ngain);
popres.s_BiasXAll = cell(nProbe,ngain);
popres.s_Xsum = cell(nProbe,ngain);
popres.s_X = cell(nProbe,ngain);
popres.s_trialID = cell(nProbe,ngain);
popres.s_runspeed = cell(nProbe,ngain);
popres.s_visspeed = cell(nProbe,ngain);
popres.s_licks = cell(nProbe,ngain);
popres.s_meanlickX = cell(nProbe,ngain);
popres.s_correcttrials = cell(nProbe,ngain);
popres.s_incorrecttrials = cell(nProbe,ngain);
popres.performance = cell(nProbe,ngain);
for iprobe = 1:nProbe
    for g = 1:ngain
        popres.PostXAll{iprobe,g} = 0;
        popres.ErrXAll{iprobe,g} = 0;
        popres.BiasXAll{iprobe,g} = 0;
        popres.Xsum{iprobe,g} = 0;
        popres.PostXAllSpd{iprobe,g} = 0;
        popres.ErrXAllSpd{iprobe,g} = 0;
        popres.BiasXAllSpd{iprobe,g} = 0;
        popres.XsumSpd{iprobe,g} = 0;
        popres.PostXAllVis{iprobe,g} = 0;
        popres.ErrXAllVis{iprobe,g} = 0;
        popres.BiasXAllVis{iprobe,g} = 0;
        popres.XsumVis{iprobe,g} = 0;
        popres.PostXAllEye{iprobe,g} = 0;
        popres.ErrXAllEye{iprobe,g} = 0;
        popres.BiasXAllEye{iprobe,g} = 0;
        popres.XsumEye{iprobe,g} = 0;
        if isfield(res,'XSumPhsAll')
            popres.ErrXPhsAll{iprobe,g} = 0;
            popres.BiasXPhsAll{iprobe,g} = 0;
            popres.XsumPhsAll{iprobe,g} = 0;
            for ispd = 1:res.nSpdPhsbins
                popres.ErrXAllPhs{iprobe,ispd,g} = 0;
                popres.BiasXAllPhs{iprobe,ispd,g} = 0;
                popres.XsumPhs{iprobe,ispd,g} = 0;
            end
        end
        popres.X{iprobe,g} = [];
        popres.runspeed{iprobe,g} = [];
        popres.visspeed{iprobe,g} = [];
        popres.thetaperiod{iprobe,g} = [];
        popres.eyeXpos{iprobe,g} = [];
        popres.licks{iprobe,g} = [];
        for kiter = 1:size(res.XSumCVO,5)
            popres.PostXAllCV{iprobe,g,kiter} = 0;
            popres.ErrXAllCV{iprobe,g,kiter} = 0;
            popres.BiasXAllCV{iprobe,g,kiter} = 0;
            popres.XsumCV{iprobe,g,kiter} = 0;
        end
        if isfield(res,'XSumPhsAll')
            for kiter = 1:size(res.XSumPhsAllCVO,5)
                popres.ErrXPhsAllCV{iprobe,g,kiter} = 0;
                popres.BiasXPhsAllCV{iprobe,g,kiter} = 0;
                popres.XsumPhsAllCV{iprobe,g,kiter} = 0;
            end
        end
        popres.ErrXpred_std{iprobe,g} = 0;
        for ianimal = 1:nanimal
            for iseries = 1:numel(expt(ianimal).series)
                if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                    if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
                        popres.s_PostXAll{ianimal,iseries,iprobe,g} = 0;
                        popres.s_ErrXAll{ianimal,iseries,iprobe,g} = 0;
                        popres.s_BiasXAll{ianimal,iseries,iprobe,g} = 0;
                        popres.s_Xsum{ianimal,iseries,iprobe,g} = 0;
                        popres.s_PostXAllSpd{ianimal,iseries,iprobe,g} = 0;
                        popres.s_ErrXAllSpd{ianimal,iseries,iprobe,g} = 0;
                        popres.s_BiasXAllSpd{ianimal,iseries,iprobe,g} = 0;
                        popres.s_XsumSpd{ianimal,iseries,iprobe,g} = 0;
                        popres.s_PostXAllVis{ianimal,iseries,iprobe,g} = 0;
                        popres.s_ErrXAllVis{ianimal,iseries,iprobe,g} = 0;
                        popres.s_BiasXAllVis{ianimal,iseries,iprobe,g} = 0;
                        popres.s_XsumVis{ianimal,iseries,iprobe,g} = 0;
                        popres.s_PostXAllEye{ianimal,iseries,iprobe,g} = 0;
                        popres.s_ErrXAllEye{ianimal,iseries,iprobe,g} = 0;
                        popres.s_BiasXAllEye{ianimal,iseries,iprobe,g} = 0;
                        popres.s_XsumEye{ianimal,iseries,iprobe,g} = 0;
                        if isfield(res,'XSumPhsAll')
                            popres.s_ErrXPhsAll{ianimal,iseries,iprobe,g} = 0;
                            popres.s_BiasXPhsAll{ianimal,iseries,iprobe,g} = 0;
                            popres.s_XsumPhsAll{ianimal,iseries,iprobe,g} = 0;
                            for ispd = 1:res.nSpdPhsbins
                                popres.s_ErrXAllPhs{ianimal,iseries,iprobe,ispd,g} = 0;
                                popres.s_BiasXAllPhs{ianimal,iseries,iprobe,ispd,g} = 0;
                                popres.s_XsumPhs{ianimal,iseries,iprobe,ispd,g} = 0;
                            end
                        end
                        popres.s_X{ianimal,iseries,iprobe,g} = [];
                        popres.s_trialID{ianimal,iseries,iprobe,g} = [];
                        popres.s_runspeed{ianimal,iseries,iprobe,g} = [];
                        popres.s_visspeed{ianimal,iseries,iprobe,g} = [];
                        popres.s_thetaperiod{ianimal,iseries,iprobe,g} = [];
                        popres.s_licks{ianimal,iseries,iprobe,g} = [];
                        popres.s_meanlickX{ianimal,iseries,iprobe,g} = [];
                        popres.s_correcttrials{ianimal,iseries,iprobe,g} = [];
                        popres.s_incorrecttrials{ianimal,iseries,iprobe,g} = [];
                        
                        for kiter = 1:size(res.XSumCVO,5)
                            popres.s_PostXAllCV{ianimal,iseries,iprobe,g,kiter} = 0;
                            popres.s_ErrXAllCV{ianimal,iseries,iprobe,g,kiter} = 0;
                            popres.s_BiasXAllCV{ianimal,iseries,iprobe,g,kiter} = 0;
                            popres.s_XsumCV{ianimal,iseries,iprobe,g,kiter} = 0;
                        end
                        if isfield(res,'XSumPhsAll')
                            for kiter = 1:size(res.XSumPhsAllCVO,5)
                                popres.s_ErrXPhsAllCV{ianimal,iseries,iprobe,g,kiter} = 0;
                                popres.s_BiasXPhsAllCV{ianimal,iseries,iprobe,g,kiter} = 0;
                                popres.s_XsumPhsAllCV{ianimal,iseries,iprobe,g,kiter} = 0;
                            end
                        end
                        
                        popres.s_performance{ianimal,iseries,iprobe,g} = [];
                    end
                end
            end
        end
    end
end

if FshowsessionFig
    f_post = figure('Name','Post');
end

for iprobe = 1:nProbe
    for ianimal = 1:nanimal
        disp(ianimal)
        nseries = 0;
        for iseries = 1:numel(expt(ianimal).series)
            if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
                    nseries = nseries + 1;
                end
            end
        end
        iseries_valid = 0;
        for iseries = 1:numel(expt(ianimal).series)
            if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
                    iseries_valid = iseries_valid + 1;
                    Xrange = max(res.X{ianimal,iseries}(~isnan(res.X{ianimal,iseries})));
                    for g = [2 1 3]
                        if ((iprobe == 1 && expt(ianimal).goodCA1dec{g}{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1dec{g}{iseries} == 1))
                            if ~isempty(res.XSum{ianimal,iseries,iprobe,g})
                                res.Phsbin{ianimal,iseries} = ones(size(res.traj{ianimal,iseries}));
                                tidx = res.tidx{ianimal,iseries,iprobe,g};% & ismember(res.trialgainchange{ianimal,iseries}, [1:5]);
                                
                                popres.runspeed{iprobe,g} = [popres.runspeed{iprobe,g} ; res.runSpeed{ianimal,iseries}(tidx)];
                                popres.visspeed{iprobe,g} = [popres.visspeed{iprobe,g} ; res.visSpeed{ianimal,iseries}(tidx)];
                                popres.eyeXpos{iprobe,g} = [popres.eyeXpos{iprobe,g} ; res.eyeXpos{ianimal,iseries}(tidx)];
                                popres.licks{iprobe,g} = [popres.licks{iprobe,g} ; res.firstgoodlick{ianimal,iseries}(tidx)];
                                
                                popres.X{iprobe,g} = [popres.X{iprobe,g} ; res.X{ianimal,iseries}(tidx)];
                                popres.PostXAll{iprobe,g} = popres.PostXAll{iprobe,g} + res.(resPosName){ianimal,iseries,iprobe,g};
                                popres.ErrXAll{iprobe,g} = popres.ErrXAll{iprobe,g} + res.(resErrName){ianimal,iseries,iprobe,g};
%                                 popres.BiasXAll{iprobe,g} = popres.BiasXAll{iprobe,g} + res.(resBiasName){ianimal,iseries,iprobe,g};
                                popres.Xsum{iprobe,g} = popres.Xsum{iprobe,g} + res.XSum{ianimal,iseries,iprobe,g};
                                
                                nfold = size(res.([resPosName 'CVO']),5);
                                for kiter = 1:nfold
                                    popres.PostXAllCV{iprobe,g,kiter} = popres.PostXAllCV{iprobe,g,kiter} + res.([resPosName 'CVO']){ianimal,iseries,iprobe,g,kiter};
                                    popres.ErrXAllCV{iprobe,g,kiter} = popres.ErrXAllCV{iprobe,g,kiter} + res.([resErrName 'CVO']){ianimal,iseries,iprobe,g,kiter};
%                                     popres.BiasXAllCV{iprobe,g,kiter} = popres.BiasXAllCV{iprobe,g,kiter} + res.([resBiasName 'CVO']){ianimal,iseries,iprobe,g,kiter};
                                    popres.XsumCV{iprobe,g,kiter} = popres.XsumCV{iprobe,g,kiter} + res.XSumCVO{ianimal,iseries,iprobe,g,kiter};
                                end
                                
%                                 res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,:,1) = nansum(res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,:,1:2),3);
%                                 res.XSumSpd{ianimal,iseries,iprobe,g}(:,:,1) = nansum(res.XSumSpd{ianimal,iseries,iprobe,g}(:,:,1:2),3);
%                                 res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,:,2) = nansum(res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,:,3:4),3);
%                                 res.XSumSpd{ianimal,iseries,iprobe,g}(:,:,2) = nansum(res.XSumSpd{ianimal,iseries,iprobe,g}(:,:,3:4),3);
%                                 res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,:,3) = nansum(res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,:,5:6),3);
%                                 res.XSumSpd{ianimal,iseries,iprobe,g}(:,:,3) = nansum(res.XSumSpd{ianimal,iseries,iprobe,g}(:,:,5:6),3);
%                                 res.PostXSumSpd{ianimal,iseries,iprobe,g} = res.PostXSumSpd{ianimal,iseries,iprobe,g}(:,:,1:3);
%                                 res.XSumSpd{ianimal,iseries,iprobe,g} = res.XSumSpd{ianimal,iseries,iprobe,g}(:,:,1:3);
                                popres.PostXAllSpd{iprobe,g} = popres.PostXAllSpd{iprobe,g} + res.([resPosName 'Spd']){ianimal,iseries,iprobe,g};
                                popres.ErrXAllSpd{iprobe,g} = popres.ErrXAllSpd{iprobe,g} + res.([resErrName 'Spd']){ianimal,iseries,iprobe,g};
%                                 popres.BiasXAllSpd{iprobe,g} = popres.BiasXAllSpd{iprobe,g} + res.([resBiasName 'Spd']){ianimal,iseries,iprobe,g};
                                popres.XsumSpd{iprobe,g} = popres.XsumSpd{iprobe,g} + res.XSumSpd{ianimal,iseries,iprobe,g};
                                
%                                 res.PostXSumEye{ianimal,iseries,iprobe,g}(:,:,1) = nansum(res.PostXSumEye{ianimal,iseries,iprobe,g}(:,:,1:2),3);
%                                 res.XSumEye{ianimal,iseries,iprobe,g}(:,:,1) = nansum(res.XSumEye{ianimal,iseries,iprobe,g}(:,:,1:2),3);
%                                 res.PostXSumEye{ianimal,iseries,iprobe,g}(:,:,2) = nansum(res.PostXSumEye{ianimal,iseries,iprobe,g}(:,:,3:4),3);
%                                 res.XSumEye{ianimal,iseries,iprobe,g}(:,:,2) = nansum(res.XSumEye{ianimal,iseries,iprobe,g}(:,:,3:4),3);
%                                 res.PostXSumEye{ianimal,iseries,iprobe,g}(:,:,3) = nansum(res.PostXSumEye{ianimal,iseries,iprobe,g}(:,:,5:6),3);
%                                 res.XSumEye{ianimal,iseries,iprobe,g}(:,:,3) = nansum(res.XSumEye{ianimal,iseries,iprobe,g}(:,:,5:6),3);
%                                 res.PostXSumEye{ianimal,iseries,iprobe,g} = res.PostXSumEye{ianimal,iseries,iprobe,g}(:,:,1:3);
%                                 res.XSumEye{ianimal,iseries,iprobe,g} = res.XSumEye{ianimal,iseries,iprobe,g}(:,:,1:3);
                                popres.PostXAllEye{iprobe,g} = popres.PostXAllEye{iprobe,g} + res.([resPosName 'Eye']){ianimal,iseries,iprobe,g};
                                popres.ErrXAllEye{iprobe,g} = popres.ErrXAllEye{iprobe,g} + res.([resErrName 'Eye']){ianimal,iseries,iprobe,g};
%                                 popres.BiasXAllEye{iprobe,g} = popres.BiasXAllEye{iprobe,g} + res.([resBiasName 'Eye']){ianimal,iseries,iprobe,g};
                                popres.XsumEye{iprobe,g} = popres.XsumEye{iprobe,g} + res.XSumEye{ianimal,iseries,iprobe,g};
                                
%                                 res.PostXSumVis{ianimal,iseries,iprobe,g}(:,:,1) = nansum(res.PostXSumVis{ianimal,iseries,iprobe,g}(:,:,1:2),3);
%                                 res.XSumVis{ianimal,iseries,iprobe,g}(:,:,1) = nansum(res.XSumVis{ianimal,iseries,iprobe,g}(:,:,1:2),3);
%                                 res.PostXSumVis{ianimal,iseries,iprobe,g}(:,:,2) = nansum(res.PostXSumVis{ianimal,iseries,iprobe,g}(:,:,3:4),3);
%                                 res.XSumVis{ianimal,iseries,iprobe,g}(:,:,2) = nansum(res.XSumVis{ianimal,iseries,iprobe,g}(:,:,3:4),3);
%                                 res.PostXSumVis{ianimal,iseries,iprobe,g}(:,:,3) = nansum(res.PostXSumVis{ianimal,iseries,iprobe,g}(:,:,5:6),3);
%                                 res.XSumVis{ianimal,iseries,iprobe,g}(:,:,3) = nansum(res.XSumVis{ianimal,iseries,iprobe,g}(:,:,5:6),3);
%                                 res.PostXSumVis{ianimal,iseries,iprobe,g} = res.PostXSumVis{ianimal,iseries,iprobe,g}(:,:,1:3);
%                                 res.XSumVis{ianimal,iseries,iprobe,g} = res.XSumVis{ianimal,iseries,iprobe,g}(:,:,1:3);
                                popres.PostXAllVis{iprobe,g} = popres.PostXAllVis{iprobe,g} + res.([resPosName 'Vis']){ianimal,iseries,iprobe,g};
                                popres.ErrXAllVis{iprobe,g} = popres.ErrXAllVis{iprobe,g} + res.([resErrName 'Vis']){ianimal,iseries,iprobe,g};
%                                 popres.BiasXAllVis{iprobe,g} = popres.BiasXAllVis{iprobe,g} + res.([resBiasName 'Vis']){ianimal,iseries,iprobe,g};
                                popres.XsumVis{iprobe,g} = popres.XsumVis{iprobe,g} + res.XSumVis{ianimal,iseries,iprobe,g};
                                
                                if isfield(res,'XSumPhsAll')
%                                     popres.thetaperiod{iprobe,g} = [popres.thetaperiod{iprobe,g} ; res.thetaperiod{ianimal,iseries,iprobe,g}(tidx)];

                                    Prange = size(res.([resErrName 'PhsAll']){ianimal,iseries,iprobe,g},1);
                                    Phsrange = size(res.([resErrName 'PhsAll']){ianimal,iseries,iprobe,g},2);
                                    if strcmp(popres.DecAveType,'DecPost') || strcmp(popres.DecAveType,'DecError')
                                        popres.ErrXPhsAll{iprobe,g} = popres.ErrXPhsAll{iprobe,g} + squeeze(nanmean(res.([resErrName 'PhsAll']){ianimal,iseries,iprobe,g}...
                                                                                                                                                        ./res.XSumPhsAll{ianimal,iseries,iprobe,g},3));
%                                         popres.BiasXPhsAll{iprobe,g} = popres.BiasXPhsAll{iprobe,g} + squeeze(nanmean(res.([resBiasName 'PhsAll']){ianimal,iseries,iprobe,g}...
%                                                                                                                                                           ./res.XSumPhsAll{ianimal,iseries,iprobe,g},3));
%                                     elseif strcmp(popres.DecAveType,'DecError')
%                                         popres.ErrXPhsAll{iprobe,g} = popres.ErrXPhsAll{iprobe,g} + squeeze(nanmean(special_smooth2D(res.([resErrName 'PhsAll']){ianimal,iseries,iprobe,g},[2/Prange 2/Phsrange],[true true])...
%                                                                                                                                                        ./special_smooth2D(res.XSumPhsAll{ianimal,iseries,iprobe,g},[2/Prange 2/Phsrange],[true true]),3));
%                                         popres.BiasXPhsAll{iprobe,g} = popres.BiasXPhsAll{iprobe,g} + squeeze(nanmean(special_smooth2D(res.([resBiasName 'PhsAll']){ianimal,iseries,iprobe,g},[2/Prange 2/Phsrange],[true true])...
%                                                                                                                                                           ./special_smooth2D(res.XSumPhsAll{ianimal,iseries,iprobe,g},[2/Prange 2/Phsrange],[true true]),3));
                                    end
                                    popres.XsumPhsAll{iprobe,g} = popres.XsumPhsAll{iprobe,g} + ones(size(res.XSumPhsAll{ianimal,iseries,iprobe,g},1),size(res.XSumPhsAll{ianimal,iseries,iprobe,g},2));
                                    
                                    for ispd = 1:res.nSpdPhsbins
                                        popres.ErrXAllPhs{iprobe,ispd,g} = popres.ErrXAllPhs{iprobe,ispd,g} + res.([resErrName 'Phs']){ianimal,iseries,iprobe,ispd,g};
%                                         popres.BiasXAllPhs{iprobe,ispd,g} = popres.BiasXAllPhs{iprobe,ispd,g} + res.([resBiasName 'Phs']){ianimal,iseries,iprobe,ispd,g};
                                        popres.XsumPhs{iprobe,ispd,g} = popres.XsumPhs{iprobe,ispd,g} + res.XSumPhs{ianimal,iseries,iprobe,ispd,g};
                                    end
                                    
                                    nfold = size(res.XSumCVO,5);
                                    for kiter = 1:nfold
%                                         popres.s_PostXPhsAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_PostXPhsAllCV{ianimal,iseries,iprobe,g,kiter} + res.ErrXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter};
%                                         popres.s_XsumPhsAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_XsumPhsAllCV{ianimal,iseries,iprobe,g,kiter} + res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter};
                                        if strcmp(popres.DecAveType,'DecPost') || strcmp(popres.DecAveType,'DecError')
                                            popres.ErrXPhsAllCV{iprobe,g,kiter} = popres.ErrXPhsAllCV{iprobe,g,kiter} + squeeze(nanmean(res.([resErrName 'PhsAllCVO']){ianimal,iseries,iprobe,g,kiter}...
                                                                                                                                                                            ./res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},3));
%                                             popres.BiasXPhsAllCV{iprobe,g,kiter} = popres.BiasXPhsAllCV{iprobe,g,kiter} + squeeze(nanmean(res.([resBiasName 'PhsAllCVO']){ianimal,iseries,iprobe,g,kiter}...
%                                                                                                                                                                               ./res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},3));
%                                         elseif strcmp(popres.DecAveType,'DecError')
%                                             popres.ErrXPhsAllCV{iprobe,g,kiter} = popres.ErrXPhsAllCV{iprobe,g,kiter} + squeeze(nanmean(special_smooth2D(res.([resErrName 'PhsAllCVO']){ianimal,iseries,iprobe,g,kiter},[2/Prange 2/Phsrange],[true true])...
%                                                                                                                                                                                              ./special_smooth2D(res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},[2/Prange 2/Phsrange],[true true]),3));
%                                             popres.BiasXPhsAllCV{iprobe,g,kiter} = popres.BiasXPhsAllCV{iprobe,g,kiter} + squeeze(nanmean(special_smooth2D(res.([resBiasName 'PhsAllCVO']){ianimal,iseries,iprobe,g,kiter},[2/Prange 2/Phsrange],[true true])...
%                                                                                                                                                                                                ./special_smooth2D(res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},[2/Prange 2/Phsrange],[true true]),3));
                                        end
                                        popres.XsumPhsAllCV{iprobe,g,kiter} = popres.XsumPhsAllCV{iprobe,g,kiter} + ones(size(res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},1),size(res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},2));
                                    end
                                end
                                
                                popres.s_runspeed{ianimal,iseries,iprobe,g} = [popres.s_runspeed{ianimal,iseries,iprobe,g} ; res.runSpeed{ianimal,iseries}(tidx)];
                                popres.s_visspeed{ianimal,iseries,iprobe,g} = [popres.s_visspeed{ianimal,iseries,iprobe,g} ; res.visSpeed{ianimal,iseries}(tidx)];
                                popres.s_licks{ianimal,iseries,iprobe,g} = [popres.s_licks{ianimal,iseries,iprobe,g} ; res.firstgoodlick{ianimal,iseries}(tidx)];
                                
                                popres.s_X{ianimal,iseries,iprobe,g} = [popres.s_X{ianimal,iseries,iprobe,g} ; res.X{ianimal,iseries}(tidx)];
                                popres.s_PostXAll{ianimal,iseries,iprobe,g} = popres.s_PostXAll{ianimal,iseries,iprobe,g} + res.(resPosName){ianimal,iseries,iprobe,g};
                                popres.s_ErrXAll{ianimal,iseries,iprobe,g} = popres.s_ErrXAll{ianimal,iseries,iprobe,g} + res.(resErrName){ianimal,iseries,iprobe,g};
%                                 popres.s_BiasXAll{ianimal,iseries,iprobe,g} = popres.s_BiasXAll{ianimal,iseries,iprobe,g} + res.(resBiasName){ianimal,iseries,iprobe,g};
                                popres.s_Xsum{ianimal,iseries,iprobe,g} = popres.s_Xsum{ianimal,iseries,iprobe,g} + res.XSum{ianimal,iseries,iprobe,g};
                                
                                nfold = size(res.XSumCVO,5);
                                for kiter = 1:nfold
                                    popres.s_PostXAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_PostXAllCV{ianimal,iseries,iprobe,g,kiter} + res.([resPosName 'CVO']){ianimal,iseries,iprobe,g,kiter};
                                    popres.s_ErrXAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_ErrXAllCV{ianimal,iseries,iprobe,g,kiter} + res.([resErrName 'CVO']){ianimal,iseries,iprobe,g,kiter};
%                                     popres.s_BiasXAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_BiasXAllCV{ianimal,iseries,iprobe,g,kiter} + res.([resBiasName 'CVO']){ianimal,iseries,iprobe,g,kiter};
                                    popres.s_XsumCV{ianimal,iseries,iprobe,g,kiter} = popres.s_XsumCV{ianimal,iseries,iprobe,g,kiter} + res.XSumCVO{ianimal,iseries,iprobe,g,kiter};
                                end
                                
                                popres.s_PostXAllSpd{ianimal,iseries,iprobe,g} = popres.s_PostXAllSpd{ianimal,iseries,iprobe,g} + res.([resPosName 'Spd']){ianimal,iseries,iprobe,g};
                                popres.s_ErrXAllSpd{ianimal,iseries,iprobe,g} = popres.s_ErrXAllSpd{ianimal,iseries,iprobe,g} + res.([resErrName 'Spd']){ianimal,iseries,iprobe,g};
%                                 popres.s_BiasXAllSpd{ianimal,iseries,iprobe,g} = popres.s_BiasXAllSpd{ianimal,iseries,iprobe,g} + res.([resBiasName 'Spd']){ianimal,iseries,iprobe,g};
                                popres.s_XsumSpd{ianimal,iseries,iprobe,g} = popres.s_XsumSpd{ianimal,iseries,iprobe,g} + res.XSumSpd{ianimal,iseries,iprobe,g};
                                
                                popres.s_PostXAllEye{ianimal,iseries,iprobe,g} = popres.s_PostXAllEye{ianimal,iseries,iprobe,g} + res.([resPosName 'Eye']){ianimal,iseries,iprobe,g};
                                popres.s_ErrXAllEye{ianimal,iseries,iprobe,g} = popres.s_ErrXAllEye{ianimal,iseries,iprobe,g} + res.([resErrName 'Eye']){ianimal,iseries,iprobe,g};
%                                 popres.s_BiasXAllEye{ianimal,iseries,iprobe,g} = popres.s_BiasXAllEye{ianimal,iseries,iprobe,g} + res.([resBiasName 'Eye']){ianimal,iseries,iprobe,g};
                                popres.s_XsumEye{ianimal,iseries,iprobe,g} = popres.s_XsumEye{ianimal,iseries,iprobe,g} + res.XSumEye{ianimal,iseries,iprobe,g};
                                
                                popres.s_PostXAllVis{ianimal,iseries,iprobe,g} = popres.s_PostXAllVis{ianimal,iseries,iprobe,g} + res.([resPosName 'Vis']){ianimal,iseries,iprobe,g};
                                popres.s_ErrXAllVis{ianimal,iseries,iprobe,g} = popres.s_ErrXAllVis{ianimal,iseries,iprobe,g} + res.([resErrName 'Vis']){ianimal,iseries,iprobe,g};
%                                 popres.s_BiasXAllVis{ianimal,iseries,iprobe,g} = popres.s_BiasXAllVis{ianimal,iseries,iprobe,g} + res.([resBiasName 'Vis']){ianimal,iseries,iprobe,g};
                                popres.s_XsumVis{ianimal,iseries,iprobe,g} = popres.s_XsumVis{ianimal,iseries,iprobe,g} + res.XSumVis{ianimal,iseries,iprobe,g};
                                
                                if isfield(res,'XSumPhsAll')
%                                     popres.s_thetaperiod{ianimal,iseries,iprobe,g} = [popres.s_thetaperiod{ianimal,iseries,iprobe,g} ; res.thetaperiod{ianimal,iseries,iprobe,g}(tidx)];
                                    
%                                     popres.s_PostXPhsAll{ianimal,iseries,iprobe,g} = popres.s_PostXPhsAll{ianimal,iseries,iprobe,g} + res.ErrXSumPhsAll{ianimal,iseries,iprobe,g};
%                                     popres.s_XsumPhsAll{ianimal,iseries,iprobe,g} = popres.s_XsumPhsAll{ianimal,iseries,iprobe,g} + res.XSumPhsAll{ianimal,iseries,iprobe,g};
                                    Prange = size(res.([resErrName 'PhsAll']){ianimal,iseries,iprobe,g},1);
                                    Phsrange = size(res.([resErrName 'PhsAll']){ianimal,iseries,iprobe,g},2);
                                    if strcmp(popres.DecAveType,'DecPost') || strcmp(popres.DecAveType,'DecError')
                                        popres.s_ErrXPhsAll{ianimal,iseries,iprobe,g} = popres.s_ErrXPhsAll{ianimal,iseries,iprobe,g} + squeeze(nanmean(res.([resErrName 'PhsAll']){ianimal,iseries,iprobe,g}...
                                                                                                                                                        ./res.XSumPhsAll{ianimal,iseries,iprobe,g},3));
%                                         popres.s_BiasXPhsAll{ianimal,iseries,iprobe,g} = popres.s_BiasXPhsAll{ianimal,iseries,iprobe,g} + squeeze(nanmean(res.([resBiasName 'PhsAll']){ianimal,iseries,iprobe,g}...
%                                                                                                                                                           ./res.XSumPhsAll{ianimal,iseries,iprobe,g},3));
%                                     elseif strcmp(popres.DecAveType,'DecError')
%                                         popres.s_ErrXPhsAll{ianimal,iseries,iprobe,g} = popres.s_ErrXPhsAll{ianimal,iseries,iprobe,g} + squeeze(nanmean(special_smooth2D(res.([resErrName 'PhsAll']){ianimal,iseries,iprobe,g},[2/Prange 2/Phsrange],[true true])...
%                                                                                                                                                        ./special_smooth2D(res.XSumPhsAll{ianimal,iseries,iprobe,g},[2/Prange 2/Phsrange],[true true]),3));
%                                         popres.s_BiasXPhsAll{ianimal,iseries,iprobe,g} = popres.s_BiasXPhsAll{ianimal,iseries,iprobe,g} + squeeze(nanmean(special_smooth2D(res.([resBiasName 'PhsAll']){ianimal,iseries,iprobe,g},[2/Prange 2/Phsrange],[true true])...
%                                                                                                                                                           ./special_smooth2D(res.XSumPhsAll{ianimal,iseries,iprobe,g},[2/Prange 2/Phsrange],[true true]),3));
                                    end
                                    popres.s_XsumPhsAll{ianimal,iseries,iprobe,g} = popres.s_XsumPhsAll{ianimal,iseries,iprobe,g} + ones(size(res.XSumPhsAll{ianimal,iseries,iprobe,g},1),size(res.XSumPhsAll{ianimal,iseries,iprobe,g},2));
                                    
                                    for ispd = 1:res.nSpdPhsbins
                                        popres.s_ErrXAllPhs{ianimal,iseries,iprobe,ispd,g} = popres.s_ErrXAllPhs{ianimal,iseries,iprobe,ispd,g} + res.([resErrName 'Phs']){ianimal,iseries,iprobe,ispd,g};
%                                         popres.s_BiasXAllPhs{ianimal,iseries,iprobe,ispd,g} = popres.s_BiasXAllPhs{ianimal,iseries,iprobe,ispd,g} + res.([resBiasName 'Phs']){ianimal,iseries,iprobe,ispd,g};
                                        popres.s_XsumPhs{ianimal,iseries,iprobe,ispd,g} = popres.s_XsumPhs{ianimal,iseries,iprobe,ispd,g} + res.XSumPhs{ianimal,iseries,iprobe,ispd,g};
                                    end
                                    
                                    nfold = size(res.XSumCVO,5);
                                    for kiter = 1:nfold
%                                         popres.s_PostXPhsAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_PostXPhsAllCV{ianimal,iseries,iprobe,g,kiter} + res.ErrXSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter};
%                                         popres.s_XsumPhsAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_XsumPhsAllCV{ianimal,iseries,iprobe,g,kiter} + res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter};
                                        if strcmp(popres.DecAveType,'DecPost') || strcmp(popres.DecAveType,'DecError')
                                            popres.s_ErrXPhsAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_ErrXPhsAllCV{ianimal,iseries,iprobe,g,kiter} + squeeze(nanmean(res.([resErrName 'PhsAllCVO']){ianimal,iseries,iprobe,g,kiter}...
                                                                                                                                                                            ./res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},3));
%                                             popres.s_BiasXPhsAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_BiasXPhsAllCV{ianimal,iseries,iprobe,g,kiter} + squeeze(nanmean(res.([resBiasName 'PhsAllCVO']){ianimal,iseries,iprobe,g,kiter}...
%                                                                                                                                                                               ./res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},3));
%                                         elseif strcmp(popres.DecAveType,'DecError')
%                                             popres.s_ErrXPhsAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_ErrXPhsAllCV{ianimal,iseries,iprobe,g,kiter} + squeeze(nanmean(special_smooth2D(res.([resErrName 'PhsAllCVO']){ianimal,iseries,iprobe,g,kiter},[2/Prange 2/Phsrange],[true true])...
%                                                                                                                                                                                              ./special_smooth2D(res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},[2/Prange 2/Phsrange],[true true]),3));
%                                             popres.s_BiasXPhsAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_BiasXPhsAllCV{ianimal,iseries,iprobe,g,kiter} + squeeze(nanmean(special_smooth2D(res.([resBiasName 'PhsAllCVO']){ianimal,iseries,iprobe,g,kiter},[2/Prange 2/Phsrange],[true true])...
%                                                                                                                                                                                                ./special_smooth2D(res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},[2/Prange 2/Phsrange],[true true]),3));
                                        end
                                        popres.s_XsumPhsAllCV{ianimal,iseries,iprobe,g,kiter} = popres.s_XsumPhsAllCV{ianimal,iseries,iprobe,g,kiter} + ones(size(res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},1),size(res.XSumPhsAllCVO{ianimal,iseries,iprobe,g,kiter},2));
                                    end
                                end
                                popres.s_trialID{ianimal,iseries,iprobe,g} = [popres.s_trialID{ianimal,iseries,iprobe,g} ; res.trialID{ianimal,iseries}(tidx)];
                                
                                %transfer the performance in Batch
                                %BayesGrandAverage
                                if sum(ismember(res.outcomeVal{ianimal,iseries}, [3]))>0
                                    earlytrials = ismember(res.outcome{ianimal,iseries}(res.gain{ianimal,iseries} == res.gainVal{ianimal,iseries}(g)), 3);
                                else
                                    earlytrials = false(size(res.tidx{ianimal,iseries,iprobe,g}));
                                end
                                if sum(ismember(res.outcomeVal{ianimal,iseries}, [4]))>0
                                    latetrials = ismember(res.outcome{ianimal,iseries}(res.gain{ianimal,iseries} == res.gainVal{ianimal,iseries}(g)), 4);
                                else
                                    latetrials = false(sum(res.gain{ianimal,iseries} == res.gainVal{ianimal,iseries}(g)),1);
                                end
                                tidx_incorrect = (earlytrials | latetrials);
                                popres.s_correcttrials{ianimal,iseries,iprobe,g} = [popres.s_correcttrials{ianimal,iseries,iprobe,g} ; unique(res.trialID{ianimal,iseries}(res.tidx{ianimal,iseries,iprobe,g}))];
                                popres.s_incorrecttrials{ianimal,iseries,iprobe,g} = [popres.s_incorrecttrials{ianimal,iseries,iprobe,g} ; unique(res.trialID{ianimal,iseries}(tidx_incorrect))];
                            end
                            
                        end
                        %transfer the performance in Batch
                                %BayesGrandAverage
                        if numel(unique([popres.s_correcttrials{ianimal,iseries,iprobe,g};popres.s_incorrecttrials{ianimal,iseries,iprobe,g}])) > 0
                            perf = numel(popres.s_correcttrials{ianimal,iseries,iprobe,g})...
                                /numel(([popres.s_correcttrials{ianimal,iseries,iprobe,g};popres.s_incorrecttrials{ianimal,iseries,iprobe,g}]));
                            popres.performance{ianimal,iseries,iprobe,g} = perf;
                        end
                        
                        dxlick = 1;
                        Xsmthlick = 100/popres.Xsmthwin;
                        maxtol_licks = 0.1;
                        amp_th_licks = 0;
                        [lickmap,x,~] = fast1Dmap(popres.s_X{ianimal,iseries,iprobe,g},popres.s_licks{ianimal,iseries,iprobe,g},dxlick,1,Xsmthlick,1);
                        [lickmapref,x,~] = fast1Dmap(popres.s_X{ianimal,iseries,iprobe,2},popres.s_licks{ianimal,iseries,iprobe,2},dxlick,1,Xsmthlick,1);
                        if ~isempty(lickmap) && ~isempty(lickmapref)
                            lickcorr = NaN(size(lickmap));
                            ishift = 0;
                            for xshift = (-floor(numel(lickmap)/2)+1):floor(numel(lickmap)/2)
                                ishift = ishift + 1;
                                lickcorr(ishift) = sum(lickmap.*circshift(lickmapref,xshift));
                            end
%                             lickmap = circshift(lickmap/sum(lickmap)*numel(x)/dxlick,round(numel(x)/dxlick/2));
%                             if ianimal == 5 && iseries == 7 && g == 3
%                                 lickmap(1:30) = 0;
%                             end
%                             [~,imax] = max(lickmap);
                            lickpos = getCircularAverage(lickcorr,popres.amp_thcorr,popres.maxtolcorr,popres.interpolcorr);%getCircularAverage(lickcorr,amp_th_licks,maxtol_licks,0.05);
%                             if ~isnan(lickpos)
%                                 lickpos = x(floor(lickpos)) + rem(lickpos,1)*(x(2) - x(1));
%                             end
                            popres.s_meanlickX{ianimal,iseries,iprobe,g} = lickpos;
                        else
                            popres.s_meanlickX{ianimal,iseries,iprobe,g} = NaN;
                        end
                        
                        
                        Prange = size(popres.s_PostXAll{ianimal,iseries,iprobe,g},1);
                        if numel(popres.s_Xsum{ianimal,iseries,iprobe,g}(:)) > 1
                            if g~=2
                                popref.ErrXAll =  popres.s_ErrXAll{ianimal,iseries,iprobe,2};
                                popref.ErrXpred =  popres.s_ErrXpred{ianimal,iseries,iprobe,2};
                                popref.meanErrX_Post = popres.s_meanErrX_Post{ianimal,iseries,iprobe,2};
                            else
                                popref = [];
                            end
                            [popres.s_PostXAll{ianimal,iseries,iprobe,g},~,popres.s_PostXpred{ianimal,iseries,iprobe,g},~,~,~,~,~]...
                                = ComputePostAve(popres.s_PostXAll{ianimal,iseries,iprobe,g}, popres.s_Xsum{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput);
                            
                            [popres.s_ErrXAll{ianimal,iseries,iprobe,g},...
                                popres.s_meanErrX_Post{ianimal,iseries,iprobe,g},...
                                popres.s_ErrXpred{ianimal,iseries,iprobe,g},...
                                popres.s_CorrXAll{ianimal,iseries,iprobe,g},...
                                popres.s_CorrXpred{ianimal,iseries,iprobe,g},...
                                popres.s_meanCorrX_Post{ianimal,iseries,iprobe,g},...
                                popres.s_meanErrX_Postdiff{ianimal,iseries,iprobe,g},...
                                popres.s_ErrXpreddiff{ianimal,iseries,iprobe,g}]...
                                = ComputePostAve(popres.s_ErrXAll{ianimal,iseries,iprobe,g}, popres.s_Xsum{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, popref);
                            
%                             [popres.s_BiasXAll{ianimal,iseries,iprobe,g},...
%                                 popres.s_meanBiasX_Post{ianimal,iseries,iprobe,g},...
%                                 popres.s_PostBiaspred{ianimal,iseries,iprobe,g},~,~,~,...
%                                  popres.s_meanBiasX_Postdiff{ianimal,iseries,iprobe,g},...
%                                 popres.s_PostBiaspreddiff{ianimal,iseries,iprobe,g}]...
%                                 = ComputePostAve(popres.s_BiasXAll{ianimal,iseries,iprobe,g}, popres.s_Xsum{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, popref);
                            
                            nfold = size(popres.s_PostXAllCV,5);
                            if nfold > 1
                                s_PostXAll_std = 0;
                                s_ErrXAll_std = 0;
                                s_meanErrX_Post_std = 0;
                                s_meanErrX_Postdiff_std = 0;
                                s_PostXpred_std = 0;
                                s_ErrXpred_std = 0;
                                s_ErrXpreddiff_std = 0;
                                s_CorrXAll_std = 0;
                                s_CorrXpred_std = 0;
                                s_meanCorr_std = 0;
                                s_BiasXAll_std = 0;
                                s_meanBiasX_Post_std = 0;
                                s_meanBiasX_Postdiff_std = 0;
                                s_PostBiaspred_std = 0;
                                s_PostBiaspreddiff_std = 0;
                                for kiter = 1:nfold
                                    [s_PostXAll_iter,~,s_PostXpred_iter,~,~,~,~,~]...
                                        = ComputePostAve(popres.s_PostXAllCV{ianimal,iseries,iprobe,g,kiter}, popres.s_XsumCV{ianimal,iseries,iprobe,g,kiter}, popres, popres.s_FsmoothOutput);
                                    
                                    [s_ErrXAll_iter,...
                                        s_meanErrX_Post_iter,...
                                        s_ErrXpred_iter,...
                                        s_CorrXAll_iter,...
                                        s_CorrXpred_iter,...
                                        s_meanCorr_iter,...
                                        s_meanErrX_Postdiff_iter,...
                                        s_ErrXpreddiff_iter]...
                                        = ComputePostAve(popres.s_ErrXAllCV{ianimal,iseries,iprobe,g,kiter}, popres.s_XsumCV{ianimal,iseries,iprobe,g,kiter}, popres, popres.s_FsmoothOutput, popref);
                                    
%                                     [s_BiasXAll_iter,...
%                                         s_meanBiasX_Post_iter,...
%                                         s_PostBiaspred_iter,~,~,~,...
%                                         s_meanBiasX_Postdiff_iter,...
%                                         s_PostBiaspreddiff_iter]...
%                                         = ComputePostAve(popres.s_BiasXAllCV{ianimal,iseries,iprobe,g}, popres.s_XsumCV{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, popref);
                                    
                                    s_PostXAll_std = s_PostXAll_std + (nfold - 1)/nfold*(s_PostXAll_iter - popres.s_PostXAll{ianimal,iseries,iprobe,g}).^2;
                                    s_ErrXAll_std = s_ErrXAll_std + (nfold - 1)/nfold*(s_ErrXAll_iter - popres.s_ErrXAll{ianimal,iseries,iprobe,g}).^2;
                                    s_meanErrX_Post_std = s_meanErrX_Post_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_meanErrX_Post_iter,2*pi/Prange*popres.s_meanErrX_Post{ianimal,iseries,iprobe,g})).^2;
                                    s_meanErrX_Postdiff_std = s_meanErrX_Postdiff_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_meanErrX_Postdiff_iter,2*pi/Prange*popres.s_meanErrX_Postdiff{ianimal,iseries,iprobe,g})).^2;
                                    s_PostXpred_std = s_PostXpred_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_PostXpred_iter,2*pi/Prange*popres.s_PostXpred{ianimal,iseries,iprobe,g})).^2;
                                    s_ErrXpred_std = s_ErrXpred_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_ErrXpred_iter,2*pi/Prange*popres.s_ErrXpred{ianimal,iseries,iprobe,g})).^2;
                                    s_ErrXpreddiff_std = s_ErrXpreddiff_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_ErrXpreddiff_iter,2*pi/Prange*popres.s_ErrXpreddiff{ianimal,iseries,iprobe,g})).^2;
                                    s_CorrXAll_std = s_CorrXAll_std + (nfold - 1)/nfold*(s_CorrXAll_iter - popres.s_CorrXAll{ianimal,iseries,iprobe,g}).^2;
                                    s_CorrXpred_std = s_CorrXpred_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_CorrXpred_iter,2*pi/Prange*popres.s_CorrXpred{ianimal,iseries,iprobe,g})).^2;
                                    s_meanCorr_std = s_meanCorr_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_meanCorr_iter,2*pi/Prange*popres.s_meanCorrX_Post{ianimal,iseries,iprobe,g})).^2;
%                                     s_BiasXAll_std = s_BiasXAll_std + (nfold - 1)/nfold*(s_BiasXAll_iter - popres.s_BiasXAll{ianimal,iseries,iprobe,g}).^2;
%                                     s_meanBiasX_Post_std = s_meanBiasX_Post_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_meanBiasX_Post_iter,2*pi/Prange*popres.s_meanBiasX_Post{ianimal,iseries,iprobe,g})).^2;
%                                     s_meanBiasX_Postdiff_std = s_meanBiasX_Postdiff_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_meanBiasX_Postdiff_iter,2*pi/Prange*popres.s_meanBiasX_Postdiff{ianimal,iseries,iprobe,g})).^2;
%                                     s_PostBiaspred_std = s_PostBiaspred_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_PostBiaspred_iter,2*pi/Prange*popres.s_PostBiaspred{ianimal,iseries,iprobe,g})).^2;
%                                     s_PostBiaspreddiff_std = s_PostBiaspreddiff_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_PostBiaspreddiff_iter,2*pi/Prange*popres.s_PostBiaspreddiff{ianimal,iseries,iprobe,g})).^2;
                                end
                                popres.s_PostXAll_SE{ianimal,iseries,iprobe,g} = sqrt(s_PostXAll_std);
                                popres.s_ErrXAll_SE{ianimal,iseries,iprobe,g} = sqrt(s_ErrXAll_std);
                                popres.s_meanErrX_Post_SE{ianimal,iseries,iprobe,g} = sqrt(s_meanErrX_Post_std);
                                popres.s_meanErrX_Postdiff_SE{ianimal,iseries,iprobe,g} = sqrt(s_meanErrX_Postdiff_std);
                                popres.s_PostXpred_SE{ianimal,iseries,iprobe,g} = sqrt(s_PostXpred_std);
                                popres.s_ErrXpred_SE{ianimal,iseries,iprobe,g} = sqrt(s_ErrXpred_std);
                                popres.s_ErrXpreddiff_SE{ianimal,iseries,iprobe,g} = sqrt(s_ErrXpreddiff_std);
                                popres.s_CorrXAll_SE{ianimal,iseries,iprobe,g} = sqrt(s_CorrXAll_std);
                                popres.s_CorrXpred_SE{ianimal,iseries,iprobe,g} = sqrt(s_CorrXpred_std);
                                popres.s_meanCorrX_Post_SE{ianimal,iseries,iprobe,g} = sqrt(s_meanCorr_std);
%                                 popres.s_BiasXAll_SE{ianimal,iseries,iprobe,g} = sqrt(s_BiasXAll_std);
%                                 popres.s_meanBiasX_Post_SE{ianimal,iseries,iprobe,g} = sqrt(s_meanBiasX_Post_std);
%                                 popres.s_meanBiasX_Postdiff_SE{ianimal,iseries,iprobe,g} = sqrt(s_meanBiasX_Postdiff_std);
%                                 popres.s_PostBiaspred_SE{ianimal,iseries,iprobe,g} = sqrt(s_PostBiaspred_std);
%                                 popres.s_PostBiaspreddiff_SE{ianimal,iseries,iprobe,g} = sqrt(s_PostBiaspreddiff_std);
                            end
                            
                            popres.meanErrXAllZscore{ianimal,iseries,iprobe,g} = (popres.s_meanErrX_Post{ianimal,iseries,iprobe,g}-popres.s_meanErrX_Post{ianimal,iseries,iprobe,2})/popres.s_meanErrX_Post_SE{ianimal,iseries,iprobe,g};
                            
                            if g~=2
                                popref.ErrXAll =  popres.s_ErrXAllSpd{ianimal,iseries,iprobe,2};
                                popref.ErrXpred =  popres.s_ErrXpredSpd{ianimal,iseries,iprobe,2};
                                popref.meanErrX_Post = popres.s_meanErrX_PostSpd{ianimal,iseries,iprobe,2};
                            else
                                popref = [];
                            end
                            
                            [popres.s_PostXAllSpd{ianimal,iseries,iprobe,g},~,popres.s_PostXpredSpd{ianimal,iseries,iprobe,g},~,~,~,~,~]...
                                = ComputePostAve(popres.s_PostXAllSpd{ianimal,iseries,iprobe,g}, popres.s_XsumSpd{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput);
                            try
                            [popres.s_ErrXAllSpd{ianimal,iseries,iprobe,g},...
                                popres.s_meanErrX_PostSpd{ianimal,iseries,iprobe,g},...
                                popres.s_ErrXpredSpd{ianimal,iseries,iprobe,g},...
                                popres.s_CorrXAllSpd{ianimal,iseries,iprobe,g},...
                                popres.s_CorrXpredSpd{ianimal,iseries,iprobe,g},...
                                popres.s_meanCorrX_PostSpd{ianimal,iseries,iprobe,g},...
                                popres.s_meanErrX_PostSpddiff{ianimal,iseries,iprobe,g},...
                                popres.s_ErrXpredSpddiff{ianimal,iseries,iprobe,g}]...
                                = ComputePostAve(popres.s_ErrXAllSpd{ianimal,iseries,iprobe,g}, popres.s_XsumSpd{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, popref);
                            catch
                                keyboard
                            end
%                             [popres.s_BiasXAllSpd{ianimal,iseries,iprobe,g},...
%                                 popres.s_meanBiasX_PostSpd{ianimal,iseries,iprobe,g},...
%                                 popres.s_PostBiaspredSpd{ianimal,iseries,iprobe,g},~,~,~,...
%                                 popres.s_meanBiasX_PostSpddiff{ianimal,iseries,iprobe,g},...
%                                 popres.s_PostBiaspredSpddiff{ianimal,iseries,iprobe,g}]...
%                                 = ComputePostAve(popres.s_BiasXAllSpd{ianimal,iseries,iprobe,g}, popres.s_XsumSpd{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, popref);
                            
                            if g~=2
                                popref.ErrXAll =  popres.s_ErrXAllEye{ianimal,iseries,iprobe,2};
                                popref.ErrXpred =  popres.s_ErrXpredEye{ianimal,iseries,iprobe,2};
                                popref.meanErrX_Post = popres.s_meanErrX_PostEye{ianimal,iseries,iprobe,2};
                            else
                                popref = [];
                            end
                            
                            [popres.s_PostXAllEye{ianimal,iseries,iprobe,g},~,popres.s_PostXpredEye{ianimal,iseries,iprobe,g},~,~,~,~,~]...
                                = ComputePostAve(popres.s_PostXAllEye{ianimal,iseries,iprobe,g}, popres.s_XsumEye{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput);
                            
                            [popres.s_ErrXAllEye{ianimal,iseries,iprobe,g},...
                                popres.s_meanErrX_PostEye{ianimal,iseries,iprobe,g},...
                                popres.s_ErrXpredEye{ianimal,iseries,iprobe,g},...
                                popres.s_CorrXAllEye{ianimal,iseries,iprobe,g},...
                                popres.s_CorrXpredEye{ianimal,iseries,iprobe,g},...
                                popres.s_meanCorrX_PostEye{ianimal,iseries,iprobe,g},...
                                popres.s_meanErrX_PostEyediff{ianimal,iseries,iprobe,g},...
                                popres.s_ErrXpredEyediff{ianimal,iseries,iprobe,g}]...
                                = ComputePostAve(popres.s_ErrXAllEye{ianimal,iseries,iprobe,g}, popres.s_XsumEye{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, popref);
                            
%                             [popres.s_BiasXAllEye{ianimal,iseries,iprobe,g},...
%                                 popres.s_meanBiasX_PostEye{ianimal,iseries,iprobe,g},...
%                                 popres.s_PostBiaspredEye{ianimal,iseries,iprobe,g},~,~,~,...
%                                 popres.s_meanBiasX_PostEyediff{ianimal,iseries,iprobe,g},...
%                                 popres.s_PostBiaspredEyediff{ianimal,iseries,iprobe,g}]...
%                                 = ComputePostAve(popres.s_BiasXAllEye{ianimal,iseries,iprobe,g}, popres.s_XsumEye{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, popref);
                            
                            if g~=2
                                popref.ErrXAll =  popres.s_ErrXAllVis{ianimal,iseries,iprobe,2};
                                popref.ErrXpred =  popres.s_ErrXpredVis{ianimal,iseries,iprobe,2};
                                popref.meanErrX_Post = popres.s_meanErrX_PostVis{ianimal,iseries,iprobe,2};
                            else
                                popref = [];
                            end
                            
                            [popres.s_PostXAllVis{ianimal,iseries,iprobe,g},~,popres.s_PostXpredVis{ianimal,iseries,iprobe,g},~,~,~,~,~]...
                                = ComputePostAve(popres.s_PostXAllVis{ianimal,iseries,iprobe,g}, popres.s_XsumVis{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput);
                            
                            [popres.s_ErrXAllVis{ianimal,iseries,iprobe,g},...
                                popres.s_meanErrX_PostVis{ianimal,iseries,iprobe,g},...
                                popres.s_ErrXpredVis{ianimal,iseries,iprobe,g},...
                                popres.s_CorrXAllVis{ianimal,iseries,iprobe,g},...
                                popres.s_CorrXpredVis{ianimal,iseries,iprobe,g},...
                                popres.s_meanCorrX_PostVis{ianimal,iseries,iprobe,g},...
                                popres.s_meanErrX_PostVisdiff{ianimal,iseries,iprobe,g},...
                                popres.s_ErrXpredVisdiff{ianimal,iseries,iprobe,g}]...
                                = ComputePostAve(popres.s_ErrXAllVis{ianimal,iseries,iprobe,g}, popres.s_XsumVis{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, popref);
                            
%                             [popres.s_BiasXAllVis{ianimal,iseries,iprobe,g},...
%                                 popres.s_meanBiasX_PostVis{ianimal,iseries,iprobe,g},...
%                                 popres.s_PostBiaspredVis{ianimal,iseries,iprobe,g},~,~,~,...
%                                 popres.s_meanBiasX_PostVisdiff{ianimal,iseries,iprobe,g},...
%                                 popres.s_PostBiaspredVisdiff{ianimal,iseries,iprobe,g}]...
%                                 = ComputePostAve(popres.s_BiasXAllVis{ianimal,iseries,iprobe,g}, popres.s_XsumVis{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, popref);
                            
                            if isfield(res,'XSumPhsAll')
                                if strcmp(popres.DecAveType,'DecPost') || strcmp(popres.DecAveType,'DecError')
                                    phsSmthWin = 2;
%                                 elseif strcmp(popres.DecAveType,'DecError')
%                                     phsSmthWin = NaN;
                                end
                                for ispd = 1:res.nSpdPhsbins
                                    [popres.s_ErrXAllPhs{ianimal,iseries,iprobe,ispd,g},~,...
                                        popres.s_ErrXpredAllPhs{ianimal,iseries,iprobe,ispd,g},...
                                        popres.s_CorrXAllPhs{ianimal,iseries,iprobe,ispd,g},...
                                        popres.s_CorrXpredAllPhs{ianimal,iseries,iprobe,ispd,g},~,~,~]...
                                        = ComputePostAve(popres.s_ErrXAllPhs{ianimal,iseries,iprobe,ispd,g}, popres.s_XsumPhs{ianimal,iseries,iprobe,ispd,g}, popres, popres.s_FsmoothOutput, [], 2, 2, true, true);
                                    
%                                     [popres.s_BiasXAllPhs{ianimal,iseries,iprobe,ispd,g},~,...
%                                         popres.s_PostBiaspredAllPhs{ianimal,iseries,iprobe,ispd,g},~,~,~,~,~]...
%                                         = ComputePostAve(popres.s_BiasXAllPhs{ianimal,iseries,iprobe,ispd,g}, popres.s_XsumPhs{ianimal,iseries,iprobe,ispd,g}, popres, popres.s_FsmoothOutput, [], 2, 2, true, true);
                                end
                                
                                [popres.s_ErrXPhsAll{ianimal,iseries,iprobe,g},~,...
                                        popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,g},...
                                        popres.s_CorrXPhsAll{ianimal,iseries,iprobe,g},...
                                        popres.s_CorrXpredPhsAll{ianimal,iseries,iprobe,g},~,~,~]...
                                        = ComputePostAve(popres.s_ErrXPhsAll{ianimal,iseries,iprobe,g}, popres.s_XsumPhsAll{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, [], phsSmthWin, phsSmthWin, true, true);% NaN, NaN, true, true);
                                    
%                                 [popres.s_BiasXPhsAll{ianimal,iseries,iprobe,g},~,...
%                                         popres.s_PostBiaspredPhsAll{ianimal,iseries,iprobe,g},~,~,~,~,~]...
%                                         = ComputePostAve(popres.s_BiasXPhsAll{ianimal,iseries,iprobe,g}, popres.s_XsumPhsAll{ianimal,iseries,iprobe,g}, popres, popres.s_FsmoothOutput, [], phsSmthWin, phsSmthWin, true, true);% NaN, NaN, true, true);
                                
                                popres.s_ErrXpredNormPhsAll{ianimal,iseries,iprobe,g} = popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,g} - mean(popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,g});
                                popres.s_ErrXpredNormPhsAll{ianimal,iseries,iprobe,g} = popres.s_ErrXpredNormPhsAll{ianimal,iseries,iprobe,g}./std(popres.s_ErrXpredNormPhsAll{ianimal,iseries,iprobe,g});
                                
                                if iprobe == 1 && g == 2
                                    nphsbins = numel(popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,2});
                                    ref = -sin((0:(nphsbins-1))/nphsbins*2*pi)';%popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,2};
                                else
                                    if expt(ianimal).goodCA1dec{2}{iseries} == 1
                                        ref = popres.s_ErrXpredPhsShiftAll{ianimal,iseries,1,2};%popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,2};
                                    else
                                        ref = 0;
                                    end
                                end
                                
                                [popres.s_ErrXPhsAllaligned{ianimal,iseries,iprobe,g},...
                                    popres.s_ErrXpredPhsAllaligned{ianimal,iseries,iprobe,g},...
                                    popres.s_CorrXPhsAllaligned{ianimal,iseries,iprobe,g},...
                                    popres.s_CorrXpredPhsAllaligned{ianimal,iseries,iprobe,g},popres.s_ErrXpredPhsShiftAll{ianimal,iseries,iprobe,g}]...
                                    = PhaseAlignement(ref,...
                                                      popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,g},...
                                                      popres.s_ErrXPhsAll{ianimal,iseries,iprobe,g},...
                                                      popres.s_CorrXPhsAll{ianimal,iseries,iprobe,g},...
                                                      popres.s_CorrXpredPhsAll{ianimal,iseries,iprobe,g});
                                
                                popres.s_ErrXpredNormPhsAllaligned{ianimal,iseries,iprobe,g} = popres.s_ErrXpredPhsAllaligned{ianimal,iseries,iprobe,g} - mean(popres.s_ErrXpredPhsAllaligned{ianimal,iseries,iprobe,g});
                                popres.s_ErrXpredNormPhsAllaligned{ianimal,iseries,iprobe,g} = popres.s_ErrXpredNormPhsAllaligned{ianimal,iseries,iprobe,g}./std(popres.s_ErrXpredNormPhsAllaligned{ianimal,iseries,iprobe,g});
                                
                                nfold = size(popres.s_ErrXPhsAllCV,5);
                                if nfold > 1
                                    s_ErrXPhsAll_std = 0;
                                    s_ErrXpredPhsAll_std = 0;
                                    s_ErrXpredNormPhsAll_std = 0;
                                    s_ErrXPhsAllaligned_std = 0;
                                    s_ErrXpredPhsAllaligned_std = 0;
                                    s_ErrXpredNormPhsAllaligned_std = 0;
                                    s_ErrXpredPhsShiftAll_std = 0;
                                    s_BiasXPhsAll_std = 0;
                                    s_PostBiaspredPhsAll_std = 0;
                                    for kiter = 1:nfold
                                        [s_ErrXPhsAll_iter,~,...
                                            s_ErrXpredPhsAll_iter,~,~,~]...
                                            = ComputePostAve(popres.s_ErrXPhsAllCV{ianimal,iseries,iprobe,g,kiter}, popres.s_XsumPhsAllCV{ianimal,iseries,iprobe,g,kiter}, popres, popres.s_FsmoothOutput, [], phsSmthWin, phsSmthWin, true, true);%NaN, NaN, true, true);%
                                        s_ErrXpredPhsAll_iter = s_ErrXpredPhsAll_iter - mean(s_ErrXpredPhsAll_iter);
                                        s_ErrXpredNormPhsAll_iter = s_ErrXpredPhsAll_iter./std(s_ErrXpredPhsAll_iter);
                                        
                                        s_ErrXPhsAll_std = s_ErrXPhsAll_std + (nfold - 1)/nfold*(s_ErrXPhsAll_iter - popres.s_ErrXPhsAll{ianimal,iseries,iprobe,g}).^2;
                                        s_ErrXpredPhsAll_std = s_ErrXpredPhsAll_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_ErrXpredPhsAll_iter,2*pi/Prange*(popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,g}-mean(popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,g})))).^2;
                                        s_ErrXpredNormPhsAll_std = s_ErrXpredNormPhsAll_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_ErrXpredNormPhsAll_iter,2*pi/Prange*popres.s_ErrXpredNormPhsAll{ianimal,iseries,iprobe,g})).^2;
                                        
                                        [s_ErrXPhsAll_aligned_iter,...
                                            s_ErrXpredPhsAll_aligned_iter,...
                                            ~,~,s_ErrXpredPhsShiftAll_iter]...
                                            = PhaseAlignement(ref,...
                                            s_ErrXpredPhsAll_iter,...
                                            s_ErrXPhsAll_iter,[],[]);
                                        s_ErrXpredPhsAll_aligned_iter = s_ErrXpredPhsAll_aligned_iter - mean(s_ErrXpredPhsAll_aligned_iter);
                                        s_ErrXpredNormPhsAll_aligned_iter = s_ErrXpredPhsAll_aligned_iter./std(s_ErrXpredPhsAll_aligned_iter);
                                        
                                        s_ErrXPhsAllaligned_std = s_ErrXPhsAllaligned_std + (nfold - 1)/nfold*(s_ErrXPhsAll_aligned_iter - popres.s_ErrXPhsAllaligned{ianimal,iseries,iprobe,g}).^2;
                                        s_ErrXpredPhsAllaligned_std = s_ErrXpredPhsAllaligned_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_ErrXpredPhsAll_aligned_iter,2*pi/Prange*(popres.s_ErrXpredPhsAllaligned{ianimal,iseries,iprobe,g}-mean(popres.s_ErrXpredPhsAllaligned{ianimal,iseries,iprobe,g})))).^2;
                                        s_ErrXpredNormPhsAllaligned_std = s_ErrXpredNormPhsAllaligned_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_ErrXpredNormPhsAll_aligned_iter,2*pi/Prange*popres.s_ErrXpredNormPhsAllaligned{ianimal,iseries,iprobe,g})).^2;
                                        s_ErrXpredPhsShiftAll_std = s_ErrXpredPhsShiftAll_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_ErrXpredPhsShiftAll_iter,2*pi/Prange*popres.s_ErrXpredPhsShiftAll{ianimal,iseries,iprobe,g})).^2;
                                        
%                                         [s_BiasXPhsAll_iter,~,...
%                                             s_PostBiaspredPhsAll_iter,~,~,~]...
%                                             = ComputePostAve(popres.s_BiasXPhsAllCV{ianimal,iseries,iprobe,g,kiter}, popres.s_XsumPhsAllCV{ianimal,iseries,iprobe,g,kiter}, popres, popres.s_FsmoothOutput, [], phsSmthWin, phsSmthWin, true, true);%NaN, NaN, true, true);%
%                                         s_BiasXPhsAll_std = s_BiasXPhsAll_std + (nfold - 1)/nfold*(s_BiasXPhsAll_iter - popres.s_BiasXPhsAll{ianimal,iseries,iprobe,g}).^2;
%                                         s_PostBiaspredPhsAll_std = s_PostBiaspredPhsAll_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*s_PostBiaspredPhsAll_iter,2*pi/Prange*(popres.s_PostBiaspredPhsAll{ianimal,iseries,iprobe,g}-mean(popres.s_PostBiaspredPhsAll{ianimal,iseries,iprobe,g})))).^2;
                                        
                                    end
                                    popres.s_ErrXPhsAll_SE{ianimal,iseries,iprobe,g} = sqrt(s_ErrXPhsAll_std);
                                    popres.s_ErrXpredPhsAll_SE{ianimal,iseries,iprobe,g} = sqrt(s_ErrXpredPhsAll_std);
                                    popres.s_ErrXpredNormPhsAll_SE{ianimal,iseries,iprobe,g} = sqrt(s_ErrXpredNormPhsAll_std);
                                    
                                    popres.s_ErrXPhsAllaligned_SE{ianimal,iseries,iprobe,g} = sqrt(s_ErrXPhsAllaligned_std);
                                    popres.s_ErrXpredPhsAllaligned_SE{ianimal,iseries,iprobe,g} = sqrt(s_ErrXpredPhsAllaligned_std);
                                    popres.s_ErrXpredNormPhsAllaligned_SE{ianimal,iseries,iprobe,g} = sqrt(s_ErrXpredNormPhsAllaligned_std);
                                    popres.s_ErrXpredPhsShiftAll_SE{ianimal,iseries,iprobe,g} = sqrt(s_ErrXpredPhsShiftAll_std);
                                    
%                                     popres.s_BiasXPhsAll_SE{ianimal,iseries,iprobe,g} = sqrt(s_BiasXPhsAll_std);
%                                     popres.s_PostBiaspredPhsAll_SE{ianimal,iseries,iprobe,g} = sqrt(s_PostBiaspredPhsAll_std);
                                end
                                
                                [~, imaxtheta] = max(popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,g}-mean(popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,g}));
                                popres.ErrXpredPhsAllZscore{ianimal,iseries,iprobe,g} = (popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,g}(imaxtheta)-mean(popres.s_ErrXpredPhsAll{ianimal,iseries,iprobe,g}))/popres.s_ErrXpredPhsAll_SE{ianimal,iseries,iprobe,g}(imaxtheta);
                                [~, imaxtheta] = max(popres.s_ErrXpredNormPhsAll{ianimal,iseries,iprobe,g}-mean(popres.s_ErrXpredNormPhsAll{ianimal,iseries,iprobe,g}));
                                popres.ErrXpredNormPhsAllZscore{ianimal,iseries,iprobe,g} = (popres.s_ErrXpredNormPhsAll{ianimal,iseries,iprobe,g}(imaxtheta)-mean(popres.s_ErrXpredNormPhsAll{ianimal,iseries,iprobe,g}))/popres.s_ErrXpredNormPhsAll_SE{ianimal,iseries,iprobe,g}(imaxtheta);
                            end
                            
                            if FshowsessionFig
                                figure(f_post)
                                subplot(4,nseries,iseries_valid)
                                hold on;plot([popres.s_meanErrX_Post{ianimal,iseries,iprobe,g}  popres.s_meanErrX_Post{ianimal,iseries,iprobe,g}],[0 3],cl{g})
                                set(gca,'PlotBoxAspectRatio', [1 1 1]);
                                title(num2str(expt(ianimal).series{iseries}));
                                subplot(4,nseries,g*nseries+iseries_valid)
                                imagesc(popres.s_PostXAll{ianimal,iseries,iprobe,g})
                                hold on;plot([0 Xrange],[0 Xrange],'k');
                                colormap(parula);%colormap(RedWhiteBlue);%
                                set(gca,'Clim',[0 2],'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
                            end
                        end
                    end
                    %                 pause
                    %                 delete(f)
                end
            end
        end
        if FshowsessionFig
            figure(f_post)
            savefig2pdf([figdirname filesep popres.probestr{iprobe} '-' expt(ianimal).animal 'post' '.pdf']);
            savefig2pdf([figdirname filesep popres.probestr{iprobe} '-' expt(ianimal).animal 'post']);
            clf(f_post);
        end
    end
    
    
    switch popres.average_type
        case 'sessions'
            for g = [2 1 3]
                popref = [];
                if g ~= 2
                    popref.ErrXAll = popres.ErrXAll{iprobe,2};
                    popref.PostXpred = popres.PostXpred{iprobe,2};
                    popref.ErrXpred = popres.ErrXpred{iprobe,2};
                    popref.meanErrX_Post = popres.meanErrX_Post{iprobe,2};
                    popref.BiasXpred = popres.BiasXpred{iprobe,2};
                    popref.meanBiasX_Post = popres.meanBiasX_Post{iprobe,2};
                end
                [popres.PostXAll{iprobe,g}, popres.ErrXAll{iprobe,g}, popres.BiasXAll{iprobe,g},...
                    popres.PostXpred{iprobe,g}, popres.ErrXpred{iprobe,g}, popres.BiasXpred{iprobe,g},...
                    popres.meanErrX_Post{iprobe,g}, popres.meanBiasX_Post{iprobe,g},...
                    popres.PostXpred_SE{iprobe,g}, popres.ErrXpred_SE{iprobe,g}, popres.BiasXpred_SE{iprobe,g},...
                    popres.meanErrX_Post_SE{iprobe,g}, popres.meanBiasX_Post_SE{iprobe,g},...
                    popres.PostXpreddiff{iprobe,g}, popres.ErrXpreddiff{iprobe,g},popres.meanErrX_Postdiff{iprobe,g},...
                    popres.PostXpreddiff_SE{iprobe,g}, popres.ErrXpreddiff_SE{iprobe,g}, popres.meanErrX_Postdiff_SE{iprobe,g},...
                    popres.CorrXAll{iprobe,g}, popres.CorrXpred{iprobe,g}, popres.meanCorrX_Post{iprobe,g},...
                    popres.CorrXpred_SE{iprobe,g}, popres.meanCorrX_Post_SE{iprobe,g}]...
                    = ComputePostGrandAve(popres.s_PostXAll(:,:,iprobe,g),popres.s_ErrXAll(:,:,iprobe,g),popres.s_BiasXAll(:,:,iprobe,g),popres.s_CorrXAll(:,:,iprobe,g),popref,popres);
                
                popref = [];
                if g ~= 2
                    popref.ErrXAll = popres.ErrXAllSpd{iprobe,2};
                    popref.PostXpred = popres.PostXpredSpd{iprobe,2};
                    popref.ErrXpred = popres.ErrXpredSpd{iprobe,2};
                    popref.meanErrX_Post = popres.meanErrX_PostSpd{iprobe,2};
                    popref.BiasXpred = popres.BiasXpredSpd{iprobe,2};
                    popref.meanBiasX_Post = popres.meanBiasX_PostSpd{iprobe,2};
                end
                [popres.PostXAllSpd{iprobe,g}, popres.ErrXAllSpd{iprobe,g}, popres.BiasXAllSpd{iprobe,g},...
                    popres.PostXpredSpd{iprobe,g}, popres.ErrXpredSpd{iprobe,g}, popres.BiasXpredSpd{iprobe,g},...
                    popres.meanErrX_PostSpd{iprobe,g}, popres.meanBiasX_PostSpd{iprobe,g},...
                    popres.PostXpredSpd_SE{iprobe,g}, popres.ErrXpredSpd_SE{iprobe,g}, popres.ErrXpredSpd_SE{iprobe,g},...
                    popres.meanErrX_PostSpd_SE{iprobe,g}, popres.meanBiasX_PostSpd_SE{iprobe,g},...
                    popres.PostXpredSpddiff{iprobe,g}, popres.ErrXpredSpddiff{iprobe,g},popres.meanErrX_PostSpddiff{iprobe,g},...
                    popres.PostXpredSpddiff_SE{iprobe,g}, popres.ErrXpredSpddiff_SE{iprobe,g}, popres.meanErrX_PostSpddiff_SE{iprobe,g},...
                    popres.CorrXAllSpd{iprobe,g}, popres.CorrXpredSpd{iprobe,g}, popres.meanCorrX_PostSpd{iprobe,g},...
                    popres.CorrXpredSpd_SE{iprobe,g}, popres.meanCorrX_PostSpd_SE{iprobe,g}]...
                    = ComputePostGrandAve(popres.s_PostXAllSpd(:,:,iprobe,g),popres.s_ErrXAllSpd(:,:,iprobe,g),popres.s_BiasXAllSpd(:,:,iprobe,g),popres.s_CorrXAllSpd(:,:,iprobe,g),popref,popres);
                
                popref = [];
                if g ~= 2
                    popref.ErrXAll = popres.ErrXAllVis{iprobe,2};
                    popref.PostXpred = popres.PostXpredVis{iprobe,2};
                    popref.ErrXpred = popres.ErrXpredVis{iprobe,2};
                    popref.meanErrX_Post = popres.meanErrX_PostVis{iprobe,2};
                    popref.BiasXpred = popres.BiasXpredVis{iprobe,2};
                    popref.meanBiasX_Post = popres.meanBiasX_PostVis{iprobe,2};
                end
                [popres.PostXAllVis{iprobe,g}, popres.ErrXAllVis{iprobe,g}, popres.BiasXAllVis{iprobe,g},...
                    popres.PostXpredVis{iprobe,g}, popres.ErrXpredVis{iprobe,g}, popres.BiasXpredVis{iprobe,g},...
                    popres.meanErrX_PostVis{iprobe,g}, popres.meanBiasX_PostVis{iprobe,g},...
                    popres.PostXpredVis_SE{iprobe,g}, popres.ErrXpredVis_SE{iprobe,g}, popres.BiasXpredVis_SE{iprobe,g},...
                    popres.meanErrX_PostVis_SE{iprobe,g}, popres.meanBiasX_PostVis_SE{iprobe,g},...
                    popres.PostXpredVisdiff{iprobe,g}, popres.ErrXpredVisdiff{iprobe,g},popres.meanErrX_PostVisdiff{iprobe,g},...
                    popres.PostXpredVisdiff_SE{iprobe,g}, popres.ErrXpredVisdiff_SE{iprobe,g}, popres.meanErrX_PostVisdiff_SE{iprobe,g},...
                    popres.CorrXAllVis{iprobe,g}, popres.CorrXpredVis{iprobe,g}, popres.meanCorrX_PostVis{iprobe,g},...
                    popres.CorrXpredVis_SE{iprobe,g}, popres.meanCorrX_PostVis_SE{iprobe,g}]...
                    = ComputePostGrandAve(popres.s_PostXAllVis(:,:,iprobe,g),popres.s_ErrXAllVis(:,:,iprobe,g),popres.s_BiasXAllVis(:,:,iprobe,g),popres.s_CorrXAllVis(:,:,iprobe,g),popref,popres);
                
                popref = [];
                if g ~= 2
                    popref.ErrXAll = popres.ErrXAllEye{iprobe,2};
                    popref.PostXpred = popres.PostXpredEye{iprobe,2};
                    popref.ErrXpred = popres.ErrXpredEye{iprobe,2};
                    popref.meanErrX_Post = popres.meanErrX_PostEye{iprobe,2};
                    popref.BiasXpred = popres.BiasXpredEye{iprobe,2};
                    popref.meanBiasX_Post = popres.meanBiasX_PostEye{iprobe,2};
                end
                [popres.PostXAllEye{iprobe,g}, popres.ErrXAllEye{iprobe,g}, popres.BiasXAllEye{iprobe,g},...
                    popres.PostXpredEye{iprobe,g}, popres.ErrXpredEye{iprobe,g}, popres.BiasXpredEye{iprobe,g},...
                    popres.meanErrX_PostEye{iprobe,g}, popres.meanBiasX_PostEye{iprobe,g},...
                    popres.PostXpredEye_SE{iprobe,g}, popres.ErrXpredEye_SE{iprobe,g}, popres.BiasXpredEye_SE{iprobe,g},...
                    popres.meanErrX_PostEye_SE{iprobe,g}, popres.meanBiasX_PostEye_SE{iprobe,g},...
                    popres.PostXpredEyediff{iprobe,g}, popres.ErrXpredEyediff{iprobe,g},popres.meanErrX_PostEyediff{iprobe,g},...
                    popres.PostXpredEyediff_SE{iprobe,g}, popres.ErrXpredEyediff_SE{iprobe,g}, popres.meanErrX_PostEyediff_SE{iprobe,g},...
                    popres.CorrXAllEye{iprobe,g}, popres.CorrXpredEye{iprobe,g}, popres.meanCorrX_PostEye{iprobe,g},...
                    popres.CorrXpredEye_SE{iprobe,g}, popres.meanCorrX_PostEye_SE{iprobe,g}]...
                    = ComputePostGrandAve(popres.s_PostXAllEye(:,:,iprobe,g),popres.s_ErrXAllEye(:,:,iprobe,g),popres.s_BiasXAllEye(:,:,iprobe,g),popres.s_CorrXAllEye(:,:,iprobe,g),popref,popres);
                
                if isfield(res,'XSumPhsAll')
                    for ispd = 1:res.nSpdPhsbins
                        popref = [];
                        if g ~= 2
                            popref.ErrXAll = popres.ErrXAllPhs{iprobe,ispd,2};
                            popref.PostXpred = popres.ErrXpredAllPhs{iprobe,ispd,2};
                            popref.ErrXpred = popres.ErrXpredAllPhs{iprobe,ispd,2};
                            popref.meanErrX_Post = popres.meanErrX_PostAllPhs{iprobe,ispd,2};
                        end
                        
                        [popres.ErrXAllPhs{iprobe,ispd,g}, ~,popres.BiasXAllPhs{iprobe,ispd,g},...
                            popres.ErrXpredAllPhs{iprobe,ispd,g}, ~,popres.BiasXpredAllPhs{iprobe,ispd,g},...
                            popres.meanErrX_PostAllPhs{iprobe,ispd,g}, popres.meanBiasX_PostAllPhs{iprobe,ispd,g},...
                            popres.ErrXpredAllPhs_SE{iprobe,ispd,g}, ~, popres.BiasXpredAllPhs_SE{iprobe,ispd,g},...
                            popres.meanErrX_PostAllPhs_SE{iprobe,ispd,g}, popres.meanBiasX_PostAllPhs_SE{iprobe,ispd,g},...
                            popres.ErrXpredAllPhsdiff{iprobe,ispd,g}, ~,popres.meanErrX_PostAllPhsdiff{iprobe,ispd,g},...
                            popres.ErrXpredAllPhsdiff_SE{iprobe,ispd,g}, ~, popres.meanErrX_PostAllPhsdiff_SE{iprobe,ispd,g},...
                            popres.CorrXAllPhs{iprobe,ispd,g}, popres.CorrXpredPhs{iprobe,ispd,g}, popres.meanCorrX_PostPhs{iprobe,ispd,g},...
                            popres.CorrXpredPhs_SE{iprobe,ispd,g}, popres.meanCorrX_PostPhs_SE{iprobe,ispd,g}]...
                            = ComputePostGrandAve(popres.s_ErrXAllPhs(:,:,iprobe,ispd,g),popres.s_ErrXAllPhs(:,:,iprobe,ispd,g),popres.s_BiasXAllPhs(:,:,iprobe,ispd,g),popres.s_CorrXAllPhs(:,:,iprobe,ispd,g),popref,popres, true, false);
                    end
                    
                    popref = [];
                    if g ~= 2
                        popref.ErrXAll = popres.ErrXPhsAll{iprobe,2};
                        popref.PostXpred = popres.ErrXpredPhsAll{iprobe,2};
                        popref.ErrXpred = popres.ErrXpredPhsAll{iprobe,2};
                        popref.meanErrX_Post = popres.meanErrX_PostPhsAll{iprobe,2};
                    end
                    [popres.ErrXPhsAll{iprobe,g}, ~, popres.BiasXPhsAll{iprobe,g},...
                        popres.ErrXpredPhsAll{iprobe,g}, ~, popres.BiasXpredPhsAll{iprobe,g},...
                        popres.meanErrX_PostPhsAll{iprobe,g}, popres.meanBiasX_PostPhsAll{iprobe,g},...
                        popres.ErrXpredPhsAll_SE{iprobe,g}, ~, popres.BiasXpredPhsAll_SE{iprobe,g},...
                        popres.meanErrX_PostPhsAll_SE{iprobe,g}, popres.meanBiasX_PostPhsAll_SE{iprobe,g},...
                        popres.ErrXpredPhsAlldiff{iprobe,g}, ~,popres.meanErrX_PostPhsAlldiff{iprobe,g},...
                        popres.ErrXpredPhsAlldiff_SE{iprobe,g}, ~, popres.meanErrX_PostPhsAlldiff_SE{iprobe,g},...
                        popres.CorrXAllPhsAll{iprobe,g}, popres.CorrXpredPhsAll{iprobe,g}, popres.meanCorrX_PostPhsAll{iprobe,g},...
                        popres.CorrXpredPhsAll_SE{iprobe,g}, popres.meanCorrX_PostPhsAll_SE{iprobe,g}]...
                        = ComputePostGrandAve(popres.s_ErrXPhsAll(:,:,iprobe,g),popres.s_ErrXPhsAll(:,:,iprobe,g),popres.s_BiasXPhsAll(:,:,iprobe,g),popres.s_CorrXPhsAll(:,:,iprobe,g),popref,popres, true, false);
                    
                    popref = [];
                    if g ~= 2
                        popref.ErrXAll = popres.ErrXPhsAll{iprobe,2};
                        popref.PostXpred = popres.ErrXpredNormPhsAll{iprobe,2};
                        popref.ErrXpred = popres.ErrXpredNormPhsAll{iprobe,2};
                        popref.meanErrX_Post = 0;
                    end
                    [~, ~, ~,...
                        popres.ErrXpredNormPhsAll{iprobe,g}, ~,~,...
                        ~,~,...
                        popres.ErrXpredNormPhsAll_SE{iprobe,g}, ~, ~,...
                        ~,~,...
                        popres.ErrXpredNormPhsAlldiff{iprobe,g}, ~,~,...
                        popres.ErrXpredNormPhsAlldiff_SE{iprobe,g}, ~, ~,...
                        ~, popres.CorrXpredNormPhsAll{iprobe,g}, ~,...
                        popres.CorrXpredNormPhsAll_SE{iprobe,g}, ~]...
                        = ComputePostGrandAve(popres.s_ErrXPhsAll(:,:,iprobe,g),popres.s_ErrXPhsAll(:,:,iprobe,g),popres.s_BiasXPhsAll(:,:,iprobe,g),popres.s_CorrXPhsAll(:,:,iprobe,g),popref,popres, true, true);
                    
                    popref = [];
                    if g ~= 2
                        popref.ErrXAll = popres.ErrXPhsAllaligned{iprobe,2};
                        popref.PostXpred = popres.ErrXpredPhsAllaligned{iprobe,2};
                        popref.ErrXpred = popres.ErrXpredPhsAllaligned{iprobe,2};
                        popref.meanErrX_Post = popres.meanErrX_PostPhsAllaligned{iprobe,2};
                    end
                    [popres.ErrXPhsAllaligned{iprobe,g}, ~,~,...
                        popres.ErrXpredPhsAllaligned{iprobe,g}, ~,~,...
                        popres.meanErrX_PostPhsAllaligned{iprobe,g},~,...
                        popres.ErrXpredPhsAllaligned_SE{iprobe,g}, ~,~,...
                        popres.meanErrX_PostPhsAllaligned_SE{iprobe,g},~,...
                        popres.ErrXpredPhsAllaligneddiff{iprobe,g}, ~,popres.meanErrX_PostPhsAllaligneddiff{iprobe,g},...
                        popres.ErrXpredPhsAllaligneddiff_SE{iprobe,g}, ~, popres.meanErrX_PostPhsAllaligneddiff_SE{iprobe,g},...
                        popres.CorrXAllPhsAllaligned{iprobe,g}, popres.CorrXpredPhsAllaligned{iprobe,g}, popres.meanCorrX_PostPhsAllaligned{iprobe,g},...
                        popres.CorrXpredPhsAllaligned_SE{iprobe,g}, popres.meanCorrX_PostPhsAllaligned_SE{iprobe,g}]...
                        = ComputePostGrandAve(popres.s_ErrXPhsAllaligned(:,:,iprobe,g),popres.s_ErrXPhsAllaligned(:,:,iprobe,g),popres.s_ErrXPhsAllaligned(:,:,iprobe,g),popres.s_CorrXPhsAllaligned(:,:,iprobe,g), popref, popres, true, false);
                    popref = [];
                    if g ~= 2
                        popref.ErrXAll = popres.ErrXPhsAllaligned{iprobe,2};
                        popref.PostXpred = popres.ErrXpredNormPhsAllaligned{iprobe,2};
                        popref.ErrXpred = popres.ErrXpredNormPhsAllaligned{iprobe,2};
                        popref.meanErrX_Post = 0;
                    end
                    [~, ~, ~,...
                        popres.ErrXpredNormPhsAllaligned{iprobe,g}, ~, ~,...
                        ~,~,...
                        popres.ErrXpredNormPhsAllaligned_SE{iprobe,g}, ~,~,...
                        ~,~,...
                        popres.ErrXpredNormPhsAllaligneddiff{iprobe,g}, ~,~,...
                        popres.ErrXpredNormPhsAllaligneddiff_SE{iprobe,g}, ~, ~,...
                        ~, popres.CorrXpredNormPhsAllaligned{iprobe,g}, ~,...
                        popres.CorrXpredNormPhsAllaligned_SE{iprobe,g}, ~]...
                        = ComputePostGrandAve(popres.s_ErrXPhsAllaligned(:,:,iprobe,g),popres.s_ErrXPhsAllaligned(:,:,iprobe,g),popres.s_ErrXPhsAllaligned(:,:,iprobe,g),popres.s_CorrXPhsAllaligned(:,:,iprobe,g),popref,popres, true, true);
                end
%                 [popres.ErrXpred_reg{iprobe,g},popres.xreg_interp] = Xpredcorrection(popres.PostXpred{iprobe,g},popres.PostXpred{iprobe,2},size(popres.PostXpred{iprobe,2},1),popres.dx);
            end
        case 'trials'
            for g = [2 1 3]
                if numel(popres.Xsum{iprobe,g}(:)) > 1
                    if g~=2
                        popref.ErrXAll =  popres.ErrXAll{iprobe,2};
                        popref.ErrXpred =  popres.ErrXpred{iprobe,2};
                        popref.meanErrX_Post = popres.meanErrX_Post{iprobe,2};
                    else
                        popref = [];
                    end
                    [popres.PostXAll{iprobe,g},~,popres.PostXpred{iprobe,g},~,~,~,~,~]...
                        = ComputePostAve(popres.PostXAll{iprobe,g}, popres.Xsum{iprobe,g}, popres, popres.FsmoothOutput);
                    
                    [popres.ErrXAll{iprobe,g},...
                        popres.meanErrX_Post{iprobe,g},...
                        popres.ErrXpred{iprobe,g},...
                        popres.CorrXAll{iprobe,g},...
                        popres.CorrXpred{iprobe,g},...
                        popres.meanCorrX_Post{iprobe,g},...
                        popres.meanErrX_Postdiff{iprobe,g},...
                        popres.ErrXpreddiff{iprobe,g}]...
                        = ComputePostAve(popres.ErrXAll{iprobe,g}, popres.Xsum{iprobe,g}, popres, popres.FsmoothOutput, popref);
                    
%                     [popres.BiasXAll{iprobe,g},...
%                         popres.meanBiasX_Post{iprobe,g},...
%                         popres.BiasXpred{iprobe,g},~,~,~,...
%                         popres.meanBiasX_Postdiff{iprobe,g},...
%                         popres.BiasXpreddiff{iprobe,g}]...
%                         = ComputePostAve(popres.BiasXAll{iprobe,g}, popres.Xsum{iprobe,g}, popres, popres.FsmoothOutput, popref);
                    
                    nfold = size(popres.PostXAllCV,3);
                    if nfold > 1
                        PostXAll_std = 0;
                        PostErrAll_std = 0;
                        meanErrX_Post_std = 0;
                        meanErrX_Postdiff_std = 0;
                        PostXpred_std = 0;
                        PostErrpred_std = 0;
                        PostErrpreddiff_std = 0;
                        PostCorrAll_std = 0;
                        PostCorrpred_std = 0;
                        meanCorr_std = 0;
                        PostBiasAll_std = 0;
                        meanBiasX_Post_std = 0;
                        meanBiasX_Postdiff_std = 0;
                        PostBiaspred_std = 0;
                        PostBiaspreddiff_std = 0;
                        for kiter = 1:nfold
                            [PostXAll_iter,~,PostXpred_iter,~,~,~,~,~]...
                                = ComputePostAve(popres.PostXAllCV{iprobe,g,kiter}, popres.XsumCV{iprobe,g,kiter}, popres, popres.FsmoothOutput);
                            
                            [PostErrAll_iter,...
                                meanErrX_Post_iter,...
                                PostErrpred_iter,...
                                PostCorrAll_iter,...
                                PostCorrpred_iter,...
                                meanCorr_iter,...
                                meanErrX_Postdiff_iter,...
                                PostErrpreddiff_iter]...
                                = ComputePostAve(popres.ErrXAllCV{iprobe,g,kiter}, popres.XsumCV{iprobe,g,kiter}, popres, popres.FsmoothOutput, popref);
                            
%                             [PostBiasAll_iter,...
%                                 meanBiasX_Post_iter,...
%                                 PostBiaspred_iter,~,~,~,...
%                                 meanBiasX_Postdiff_iter,...
%                                 PostBiaspreddiff_iter]...
%                                 = ComputePostAve(popres.BiasXAllCV{iprobe,g}, popres.XsumCV{iprobe,g}, popres, popres.FsmoothOutput, popref);
                            
                            PostXAll_std = PostXAll_std + (nfold - 1)/nfold*(PostXAll_iter - popres.PostXAll{iprobe,g}).^2;
                            PostErrAll_std = PostErrAll_std + (nfold - 1)/nfold*(PostErrAll_iter - popres.ErrXAll{iprobe,g}).^2;
                            meanErrX_Post_std = meanErrX_Post_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*meanErrX_Post_iter,2*pi/Prange*popres.meanErrX_Post{iprobe,g})).^2;
                            meanErrX_Postdiff_std = meanErrX_Postdiff_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*meanErrX_Postdiff_iter,2*pi/Prange*popres.meanErrX_Postdiff{iprobe,g})).^2;
                            PostXpred_std = PostXpred_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostXpred_iter,2*pi/Prange*popres.PostXpred{iprobe,g})).^2;
                            PostErrpred_std = PostErrpred_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpred_iter,2*pi/Prange*popres.ErrXpred{iprobe,g})).^2;
                            PostErrpreddiff_std = PostErrpreddiff_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpreddiff_iter,2*pi/Prange*popres.ErrXpreddiff{iprobe,g})).^2;
                            PostCorrAll_std = PostCorrAll_std + (nfold - 1)/nfold*(PostCorrAll_iter - popres.CorrXAll{iprobe,g}).^2;
                            PostCorrpred_std = PostCorrpred_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostCorrpred_iter,2*pi/Prange*popres.CorrXpred{iprobe,g})).^2;
                            meanCorr_std = meanCorr_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*meanCorr_iter,2*pi/Prange*popres.meanCorrX_Post{iprobe,g})).^2;
%                             PostBiasAll_std = PostBiasAll_std + (nfold - 1)/nfold*(PostBiasAll_iter - popres.BiasXAll{iprobe,g}).^2;
%                             meanBiasX_Post_std = meanBiasX_Post_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*meanBiasX_Post_iter,2*pi/Prange*popres.meanBiasX_Post{iprobe,g})).^2;
%                             meanBiasX_Postdiff_std = meanBiasX_Postdiff_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*meanBiasX_Postdiff_iter,2*pi/Prange*popres.meanBiasX_Postdiff{iprobe,g})).^2;
%                             PostBiaspred_std = PostBiaspred_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostBiaspred_iter,2*pi/Prange*popres.BiasXpred{iprobe,g})).^2;
%                             PostBiaspreddiff_std = PostBiaspreddiff_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostBiaspreddiff_iter,2*pi/Prange*popres.BiasXpreddiff{iprobe,g})).^2;
                        end
                        popres.PostXAll_SE{iprobe,g} = sqrt(PostXAll_std);
                        popres.ErrXAll_SE{iprobe,g} = sqrt(PostErrAll_std);
                        popres.meanErrX_Post_SE{iprobe,g} = sqrt(meanErrX_Post_std);
                        popres.meanErrX_Postdiff_SE{iprobe,g} = sqrt(meanErrX_Postdiff_std);
                        popres.PostXpred_SE{iprobe,g} = sqrt(PostXpred_std);
                        popres.ErrXpred_SE{iprobe,g} = sqrt(PostErrpred_std);
                        popres.ErrXpreddiff_SE{iprobe,g} = sqrt(PostErrpreddiff_std);
                        popres.CorrXAll_SE{iprobe,g} = sqrt(PostCorrAll_std);
                        popres.CorrXpred_SE{iprobe,g} = sqrt(PostCorrpred_std);
                        popres.meanCorrX_Post_SE{iprobe,g} = sqrt(meanCorr_std);
%                         popres.BiasXAll_SE{iprobe,g} = sqrt(PostBiasAll_std);
%                         popres.meanBiasX_Post_SE{iprobe,g} = sqrt(meanBiasX_Post_std);
%                         popres.meanBiasX_Postdiff_SE{iprobe,g} = sqrt(meanBiasX_Postdiff_std);
%                         popres.BiasXpred_SE{iprobe,g} = sqrt(PostBiaspred_std);
%                         popres.BiasXpreddiff_SE{iprobe,g} = sqrt(PostBiaspreddiff_std);
                    end
                    
                    if g~=2
                        popref.ErrXAll =  popres.ErrXAllSpd{iprobe,2};
                        popref.ErrXpred =  popres.ErrXpredSpd{iprobe,2};
                        popref.meanErrX_Post = popres.meanErrX_PostSpd{iprobe,2};
                    else
                        popref = [];
                    end
                    
                    [popres.PostXAllSpd{iprobe,g},~,popres.PostXpredSpd{iprobe,g},~,~,~,~,~]...
                        = ComputePostAve(popres.PostXAllSpd{iprobe,g}, popres.XsumSpd{iprobe,g}, popres, popres.FsmoothOutput);
                    
                    [popres.ErrXAllSpd{iprobe,g},...
                        popres.meanErrX_PostSpd{iprobe,g},...
                        popres.ErrXpredSpd{iprobe,g},...
                        popres.CorrXAllSpd{iprobe,g},...
                        popres.CorrXpredSpd{iprobe,g},...
                        popres.meanCorrX_PostSpd{iprobe,g},...
                        popres.meanErrX_PostSpddiff{iprobe,g},...
                        popres.ErrXpredSpddiff{iprobe,g}]...
                        = ComputePostAve(popres.ErrXAllSpd{iprobe,g}, popres.XsumSpd{iprobe,g}, popres, popres.FsmoothOutput, popref);
                    
%                     [popres.BiasXAllSpd{iprobe,g},...
%                         popres.meanBiasX_PostSpd{iprobe,g},...
%                         popres.BiasXpredSpd{iprobe,g},~,~,~,...
%                         popres.meanBiasX_PostSpddiff{iprobe,g},...
%                         popres.BiasXpredSpddiff{iprobe,g}]...
%                         = ComputePostAve(popres.BiasXAllSpd{iprobe,g}, popres.XsumSpd{iprobe,g}, popres, popres.FsmoothOutput, popref);
                    
                    if g~=2
                        popref.ErrXAll =  popres.ErrXAllEye{iprobe,2};
                        popref.ErrXpred =  popres.ErrXpredEye{iprobe,2};
                        popref.meanErrX_Post = popres.meanErrX_PostEye{iprobe,2};
                    else
                        popref = [];
                    end
                    
                    [popres.PostXAllEye{iprobe,g},~,popres.PostXpredEye{iprobe,g},~,~,~,~,~]...
                        = ComputePostAve(popres.PostXAllEye{iprobe,g}, popres.XsumEye{iprobe,g}, popres, popres.FsmoothOutput);
                    
                    [popres.ErrXAllEye{iprobe,g},...
                        popres.meanErrX_PostEye{iprobe,g},...
                        popres.ErrXpredEye{iprobe,g},...
                        popres.CorrXAllEye{iprobe,g},...
                        popres.CorrXpredEye{iprobe,g},...
                        popres.meanCorrX_PostEye{iprobe,g},...
                        popres.meanErrX_PostEyediff{iprobe,g},...
                        popres.ErrXpredEyediff{iprobe,g}]...
                        = ComputePostAve(popres.ErrXAllEye{iprobe,g}, popres.XsumEye{iprobe,g}, popres, popres.FsmoothOutput, popref);
                    
%                     [popres.BiasXAllEye{iprobe,g},...
%                         popres.meanBiasX_PostEye{iprobe,g},...
%                         popres.BiasXpredEye{iprobe,g},~,~,~,...
%                         popres.meanBiasX_PostEyediff{iprobe,g},...
%                         popres.BiasXpredEyediff{iprobe,g}]...
%                         = ComputePostAve(popres.BiasXAllEye{iprobe,g}, popres.XsumEye{iprobe,g}, popres, popres.FsmoothOutput, popref);
                    
                    
                    if g~=2
                        popref.ErrXAll =  popres.ErrXAllVis{iprobe,2};
                        popref.ErrXpred =  popres.ErrXpredVis{iprobe,2};
                        popref.meanErrX_Post = popres.meanErrX_PostVis{iprobe,2};
                    else
                        popref = [];
                    end
                    
                    [popres.PostXAllVis{iprobe,g},~,popres.PostXpredVis{iprobe,g},~,~,~,~,~]...
                        = ComputePostAve(popres.PostXAllVis{iprobe,g}, popres.XsumVis{iprobe,g}, popres, popres.FsmoothOutput);
                    
                    [popres.ErrXAllVis{iprobe,g},...
                        popres.meanErrX_PostVis{iprobe,g},...
                        popres.ErrXpredVis{iprobe,g},...
                        popres.CorrXAllVis{iprobe,g},...
                        popres.CorrXpredVis{iprobe,g},...
                        popres.meanCorrX_PostVis{iprobe,g},...
                        popres.meanErrX_PostVisdiff{iprobe,g},...
                        popres.ErrXpredVisdiff{iprobe,g}]...
                        = ComputePostAve(popres.ErrXAllVis{iprobe,g}, popres.XsumVis{iprobe,g}, popres, popres.FsmoothOutput, popref);
                    
%                     [popres.BiasXAllVis{iprobe,g},...
%                         popres.meanBiasX_PostVis{iprobe,g},...
%                         popres.BiasXpredVis{iprobe,g},~,~,~,...
%                         popres.meanBiasX_PostVisdiff{iprobe,g},...
%                         popres.BiasXpredVisdiff{iprobe,g}]...
%                         = ComputePostAve(popres.BiasXAllVis{iprobe,g}, popres.XsumVis{iprobe,g}, popres, popres.FsmoothOutput, popref);
                   
                    if isfield(res,'XSumPhsAll')
                        if strcmp(popres.DecAveType,'DecPost') || strcmp(popres.DecAveType,'DecError')
                            phsSmthWin = 2;
%                         elseif strcmp(popres.DecAveType,'DecError')
%                             phsSmthWin = NaN;
                        end
                        for ispd = 1:res.nSpdPhsbins
                            [popres.ErrXAllPhs{iprobe,ispd,g},~,...
                                popres.ErrXpredAllPhs{iprobe,ispd,g},...
                                popres.CorrXAllPhs{iprobe,ispd,g},...
                                popres.CorrXpredAllPhs{iprobe,ispd,g},~,~,~]...
                                = ComputePostAve(popres.ErrXAllPhs{iprobe,ispd,g}, popres.XsumPhs{iprobe,ispd,g}, popres, popres.FsmoothOutput, [], 2, 2, true, true);
                            
%                             [popres.BiasXAllPhs{iprobe,ispd,g},~,...
%                                 popres.BiasXpredAllPhs{iprobe,ispd,g},~,~,~,~,~]...
%                                 = ComputePostAve(popres.BiasXAllPhs{iprobe,ispd,g}, popres.XsumPhs{iprobe,ispd,g}, popres, popres.FsmoothOutput, [], 2, 2, true, true);
                        end
                        
                        [popres.ErrXPhsAll{iprobe,g},~,...
                            popres.ErrXpredPhsAll{iprobe,g},...
                            popres.CorrXPhsAll{iprobe,g},...
                            popres.CorrXpredPhsAll{iprobe,g},~,~,~]...
                            = ComputePostAve(popres.ErrXPhsAll{iprobe,g}, popres.XsumPhsAll{iprobe,g}, popres, popres.FsmoothOutput, [], phsSmthWin, phsSmthWin, true, true);% NaN, NaN, true, true);
                        
%                         [popres.BiasXPhsAll{iprobe,g},~,...
%                             popres.BiasXpredPhsAll{iprobe,g},~,~,~,~,~]...
%                             = ComputePostAve(popres.BiasXPhsAll{iprobe,g}, popres.XsumPhsAll{iprobe,g}, popres, popres.FsmoothOutput, [], phsSmthWin, phsSmthWin, true, true);% NaN, NaN, true, true);
                        
                        popres.ErrXpredNormPhsAll{iprobe,g} = popres.ErrXpredPhsAll{iprobe,g} - mean(popres.ErrXpredPhsAll{iprobe,g});
                        popres.ErrXpredNormPhsAll{iprobe,g} = popres.ErrXpredNormPhsAll{iprobe,g}./std(popres.ErrXpredNormPhsAll{iprobe,g});
                        
                        if iprobe == 1 && g == 2
                            nphsbins = numel(popres.ErrXpredPhsAll{iprobe,2});
                            ref = -sin((0:(nphsbins-1))/nphsbins*2*pi)';%popres.ErrXpredPhsAll{iprobe,2};
                        else
                            if expt(ianimal).goodCA1dec{2}{iseries} == 1
                                ref = popres.ErrXpredPhsShiftAll{1,2};%popres.ErrXpredPhsAll{iprobe,2};
                            else
                                ref = 0;
                            end
                        end
                        try
                        [popres.ErrXPhsAllaligned{iprobe,g},...
                            popres.ErrXpredPhsAllaligned{iprobe,g},...
                            popres.CorrXPhsAllaligned{iprobe,g},...
                            popres.CorrXpredPhsAllaligned{iprobe,g},popres.ErrXpredPhsShiftAll{iprobe,g}]...
                            = PhaseAlignement(ref,...
                            popres.ErrXpredPhsAll{iprobe,g},...
                            popres.ErrXPhsAll{iprobe,g},...
                            popres.CorrXPhsAll{iprobe,g},...
                            popres.CorrXpredPhsAll{iprobe,g});
                        catch
                            keyboard
                        end
                        popres.ErrXpredNormPhsAllaligned{iprobe,g} = popres.ErrXpredPhsAllaligned{iprobe,g} - mean(popres.ErrXpredPhsAllaligned{iprobe,g});
                        popres.ErrXpredNormPhsAllaligned{iprobe,g} = popres.ErrXpredNormPhsAllaligned{iprobe,g}./std(popres.ErrXpredNormPhsAllaligned{iprobe,g});
                        
                        
                        nfold = size(popres.ErrXPhsAllCV,3);
                        if nfold > 1
                            PostErrPhsAll_std = 0;
                            PostErrpredPhsAll_std = 0;
                            PostErrpredNormPhsAll_std = 0;
                            PostErrPhsAllaligned_std = 0;
                            PostErrpredPhsAllaligned_std = 0;
                            PostErrpredNormPhsAllaligned_std = 0;
                            PostErrpredPhsShiftAll_std = 0;
                            PostBiasPhsAll_std = 0;
                            PostBiaspredPhsAll_std = 0;
                            for kiter = 1:nfold
                                [PostErrPhsAll_iter,~,...
                                    PostErrpredPhsAll_iter,~,~,~]...
                                    = ComputePostAve(popres.ErrXPhsAllCV{iprobe,g,kiter}, popres.XsumPhsAllCV{iprobe,g,kiter}, popres, popres.FsmoothOutput, [], phsSmthWin, phsSmthWin, true, true);%NaN, NaN, true, true);%
                                PostErrpredPhsAll_iter = PostErrpredPhsAll_iter - mean(PostErrpredPhsAll_iter);
                                PostErrpredNormPhsAll_iter = PostErrpredPhsAll_iter./std(PostErrpredPhsAll_iter);
                                
                                PostErrPhsAll_std = PostErrPhsAll_std + (nfold - 1)/nfold*(PostErrPhsAll_iter - popres.ErrXPhsAll{iprobe,g}).^2;
                                PostErrpredPhsAll_std = PostErrpredPhsAll_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpredPhsAll_iter,2*pi/Prange*(popres.ErrXpredPhsAll{iprobe,g}-mean(popres.ErrXpredPhsAll{iprobe,g})))).^2;
                                PostErrpredNormPhsAll_std = PostErrpredNormPhsAll_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpredNormPhsAll_iter,2*pi/Prange*popres.ErrXpredNormPhsAll{iprobe,g})).^2;
                                
                                [PostErrPhsAll_aligned_iter,...
                                    PostErrpredPhsAll_aligned_iter,...
                                    ~,~,PostErrpredPhsShiftAll_iter]...
                                    = PhaseAlignement(ref,...
                                    PostErrpredPhsAll_iter,...
                                    PostErrPhsAll_iter,[],[]);
                                PostErrpredPhsAll_aligned_iter = PostErrpredPhsAll_aligned_iter - mean(PostErrpredPhsAll_aligned_iter);
                                PostErrpredNormPhsAll_aligned_iter = PostErrpredPhsAll_aligned_iter./std(PostErrpredPhsAll_aligned_iter);
                                
                                PostErrPhsAllaligned_std = PostErrPhsAllaligned_std + (nfold - 1)/nfold*(PostErrPhsAll_aligned_iter - popres.ErrXPhsAllaligned{iprobe,g}).^2;
                                PostErrpredPhsAllaligned_std = PostErrpredPhsAllaligned_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpredPhsAll_aligned_iter,2*pi/Prange*(popres.ErrXpredPhsAllaligned{iprobe,g}-mean(popres.ErrXpredPhsAllaligned{iprobe,g})))).^2;
                                PostErrpredNormPhsAllaligned_std = PostErrpredNormPhsAllaligned_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpredNormPhsAll_aligned_iter,2*pi/Prange*popres.ErrXpredNormPhsAllaligned{iprobe,g})).^2;
                                PostErrpredPhsShiftAll_std = PostErrpredPhsShiftAll_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpredPhsShiftAll_iter,2*pi/Prange*popres.ErrXpredPhsShiftAll{iprobe,g})).^2;
                                
%                                 [PostBiasPhsAll_iter,~,...
%                                     PostBiaspredPhsAll_iter,~,~,~]...
%                                     = ComputePostAve(popres.BiasXPhsAllCV{iprobe,g,kiter}, popres.XsumPhsAllCV{iprobe,g,kiter}, popres, popres.FsmoothOutput, [], phsSmthWin, phsSmthWin, true, true);%NaN, NaN, true, true);%
%                                 PostBiasPhsAll_std = PostBiasPhsAll_std + (nfold - 1)/nfold*(PostBiasPhsAll_iter - popres.BiasXPhsAll{iprobe,g}).^2;
%                                 PostBiaspredPhsAll_std = PostBiaspredPhsAll_std + (nfold - 1)/nfold*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostBiaspredPhsAll_iter,2*pi/Prange*(popres.BiasXpredPhsAll{iprobe,g}-mean(popres.BiasXpredPhsAll{iprobe,g})))).^2;
                                
                            end
                            popres.ErrXPhsAll_SE{iprobe,g} = sqrt(PostErrPhsAll_std);
                            popres.ErrXpredPhsAll_SE{iprobe,g} = sqrt(PostErrpredPhsAll_std);
                            popres.ErrXpredNormPhsAll_SE{iprobe,g} = sqrt(PostErrpredNormPhsAll_std);
                            
                            popres.ErrXPhsAllaligned_SE{iprobe,g} = sqrt(PostErrPhsAllaligned_std);
                            popres.ErrXpredPhsAllaligned_SE{iprobe,g} = sqrt(PostErrpredPhsAllaligned_std);
                            popres.ErrXpredNormPhsAllaligned_SE{iprobe,g} = sqrt(PostErrpredNormPhsAllaligned_std);
                            popres.ErrXpredPhsShiftAll_SE{iprobe,g} = sqrt(PostErrpredPhsShiftAll_std);
                            
%                             popres.BiasXPhsAll_SE{iprobe,g} = sqrt(PostBiasPhsAll_std);
%                             popres.BiasXpredPhsAll_SE{iprobe,g} = sqrt(PostBiaspredPhsAll_std);
                        end
                    end
                end
%                 [popres.ErrXpred_reg{iprobe,g},popres.xreg_interp] = Xpredcorrection(popres.PostXpred{iprobe,g},popres.PostXpred{iprobe,2},size(popres.PostXpred{iprobe,2},1),popres.dx);
            end
    end
end
end

function [PostErrAll_out, meanErrX_Post_out, PostErrAllpred_out, PostCorrAll_out, PostCorrAllpred_out, meanCorrX_Post_out,...
          meanErrX_Postdiff_out, PostErrAllpreddiff_out] = ComputePostAve(PostErrSum0, Xsum0, params, Fsmthoutput, popref, Xsmthbin, Ysmthbin, XFcircular, YFcircular)
    Prange = size(PostErrSum0,1);
    Xrange = size(PostErrSum0,2);
    if nargin < 5
        popref = [];
    end
    if nargin < 6
        Xsmthbin = params.Xsmthwin;
        Ysmthbin = params.Xsmthwin;%NaN;
        XFcircular = true;
        YFcircular = true;
    end
    if nargin < 8
        XFcircular = true;
        YFcircular = true;
    end
    for k = 1:size(PostErrSum0,3)
        PostErrSum = PostErrSum0(:,:,k);
        Xsum = Xsum0(:,:,k);
        if size(Xsum,1) == 1
            if strcmp(params.DecAveType,'DecPost') || strcmp(params.DecAveType,'DecError')
                PostErrAll = special_smooth2D(PostErrSum./repmat(Xsum,[Prange 1]), [Ysmthbin/size(PostErrSum,1) Xsmthbin/size(PostErrSum,2)],[YFcircular XFcircular]);
                if ~Fsmthoutput
                    PostErrAll_o = PostErrSum./repmat(Xsum,[Prange 1]);
                else
                    PostErrAll_o = PostErrAll;
                end
%             elseif strcmp(params.DecAveType,'DecError')
%                 PostErrAll = special_smooth2D(PostErrSum, [Ysmthbin/size(PostErrSum,1) Xsmthbin/size(PostErrSum,2)],[YFcircular XFcircular])...
%                              ./special_smooth2D(repmat(Xsum,[Prange 1]), [Ysmthbin/size(PostErrSum,1) Xsmthbin/size(PostErrSum,2)],[YFcircular XFcircular]);
            end
%             PostErrAll = special_smooth2D(PostErrSum, [Ysmthbin/size(PostErrSum,1) Xsmthbin/size(PostErrSum,2)],[YFcircular XFcircular])./special_smooth2D(repmat(Xsum,[Prange 1]), [Ysmthbin/size(PostErrSum,1) Xsmthbin/size(PostErrSum,2)],[YFcircular XFcircular]);
        else
            if strcmp(params.DecAveType,'DecPost') || strcmp(params.DecAveType,'DecError')
                PostErrAll = special_smooth2D(PostErrSum./Xsum, [Ysmthbin/size(PostErrSum,1) Xsmthbin/size(PostErrSum,2)],[YFcircular XFcircular]);
                if ~Fsmthoutput
                    PostErrAll_o = PostErrSum./Xsum;
                else
                    PostErrAll_o = PostErrAll;
                end
%             elseif strcmp(params.DecAveType,'DecError')
%                 PostErrAll = special_smooth2D(PostErrSum, [Ysmthbin/size(PostErrSum,1) Xsmthbin/size(PostErrSum,2)],[YFcircular XFcircular])...
%                              ./special_smooth2D(Xsum, [Ysmthbin/size(PostErrSum,1) Xsmthbin/size(PostErrSum,2)],[YFcircular XFcircular]);
            end
%             PostErrAll = special_smooth2D(PostErrSum, [Ysmthbin/size(PostErrSum,1) Xsmthbin/size(PostErrSum,2)],[YFcircular XFcircular])./special_smooth2D(Xsum, [Ysmthbin/size(PostErrSum,1) Xsmthbin/size(PostErrSum,2)],[YFcircular XFcircular]);%smooth2D(PostXSum, params.lambdaSmooth)./smooth2D(Xsum, params.lambdaSmooth);
        end
        
        if params.dx < 1
            [xorig,yorig] = meshgrid(1:3*Xrange,1:3*Prange);
            [xinterp,yinterp] = meshgrid(params.dx:params.dx:3*Xrange,params.dx:params.dx:3*Prange);
            Post_interp = interp2(xorig,yorig,repmat(PostErrAll,[3 3]),xinterp,yinterp);
            Post_interp = Post_interp((size(yinterp,1)/3+1):(2*(size(yinterp,1)/3)),(size(xinterp,2)/3+1):(2*(size(xinterp,2)/3)));
            PostErrAll = Post_interp;
        end
        
        if ~isempty(popref) && ~Fsmthoutput
            popref.ErrXAll(:,:,k) = special_smooth2D(popref.ErrXAll(:,:,k), [Ysmthbin/size(popref.ErrXAll,1) Xsmthbin/size(popref.ErrXAll,2)],[YFcircular XFcircular]);
        end
        
        if ~isempty(popref)
            PostCorrAll = zeros(size(PostErrAll));
            errmap = PostErrAll;
            errmapref = popref.ErrXAll(:,:,k);
            for i=1:size(errmap,2)
                errmap(:,i) = errmap(:,i)-mean(errmap(:,i));
                errmap(:,i) = errmap(:,i)/(mean(errmap(:,i).^2))^0.5;
                errmapref(:,i) = errmapref(:,i)-mean(errmapref(:,i));
                errmapref(:,i) = errmapref(:,i)/(mean(errmapref(:,i).^2))^0.5;
            end
            ishift = 0;
            for xshift = -floor(size(errmap,2)/2)+1:floor(size(errmap,2)/2)
                ishift = ishift+1;
                PostCorrAll(ishift,:) = diag(errmap'*circshift(errmapref,xshift,1)/size(errmap,1));
            end
        else
            PostCorrAll = zeros(size(PostErrAll));
            errmap = PostErrAll;
            errmapref = nanmean(PostErrAll,2);
            errmapref = errmapref-nanmean(errmapref);
            errmapref = errmapref./nanmean(errmapref.^2)^0.5;
            for i=1:size(errmap,2)
                errmap(:,i) = errmap(:,i)-nanmean(errmap(:,i));
                errmap(:,i) = errmap(:,i)/(nanmean(errmap(:,i).^2))^0.5;
            end
            ishift = 0;
            for xshift = -floor(size(errmap,1)/2)+1:floor(size(errmap,1)/2)
                ishift = ishift+1;
                PostCorrAll(ishift,:) = errmap'*circshift(errmapref,xshift,1)/size(errmap,1);
            end            
        end
        
        PostErrAllpred = params.dx*(getCircularAverage(PostErrAll,params.amp_th,params.maxtol,params.interpol));
        PostCorrAlltemp = PostCorrAll;
%         PostCorrAlltemp([1:floor(size(PostErrAll,2)/4) floor(size(PostErrAll,2)*3/4):end],:) = 0;
        PostCorrAllpred = params.dx*(getCircularAverage(PostCorrAlltemp,params.amp_thcorr,params.maxtolcorr,params.interpolcorr));
%         PostCorrAllpred = ones(size(PostErrAllpred));
        
        MeanXerr = nanmean(PostErrAll,2);
        if sum(isnan(MeanXerr)) <= 2 && sum(MeanXerr>params.amp_th) > 0 && numel(MeanXerr) > 1
            [~, imax] = max(MeanXerr);
            meanErrX_Post = params.dx*(getCircularAverage(MeanXerr,params.amp_th,params.maxtol,params.interpol));%nanmean(getCircularAverage(popres.s_ErrXAll{ianimal,iseries,iprobe,g},amp_th_post,maxtol_post));%
        else
            meanErrX_Post = NaN;
        end

        MeanXerr = nanmean(PostErrAll,2);
        MeanXerr = MeanXerr-nanmean(MeanXerr);
        MeanXerr = MeanXerr/nanmean(MeanXerr.^2)^0.5;
        if ~isempty(popref)
            MeanXerrref = nanmean(squeeze(popref.ErrXAll(:,:,k)),2);
        else
            MeanXerrref = nanmean(PostErrAll,2);
        end
        MeanXerrref = MeanXerrref-nanmean(MeanXerrref);
        MeanXerrref = MeanXerrref/nanmean(MeanXerrref.^2)^0.5;
        if sum(isnan(MeanXerr)) <= 2 && sum(MeanXerr>params.amp_th) > 0 && numel(MeanXerr) > 1
            meanCorrX_Post = zeros(numel(MeanXerr),1);
            ishift = 0;
            for xshift = -floor(numel(MeanXerr)/2)+1:floor(numel(MeanXerr)/2)
                ishift = ishift+1;
                meanCorrX_Post(ishift) = MeanXerr(:)'*circshift(MeanXerrref(:),xshift);
            end
            meanCorrX_Post = params.dx*(getCircularAverage(meanCorrX_Post(:),params.amp_thcorr,params.maxtolcorr,params.interpolcorr));
        else
            meanCorrX_Post = NaN;
        end
        
        if ~isempty(popref)
            PostErrAllpreddiff = Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrAllpred,2*pi/Prange*popref.ErrXpred(:,k));
        else
            PostErrAllpreddiff = NaN;
        end
        if ~isempty(popref)
            meanErrX_Postdiff = Prange/(2*pi)*circ_dist(2*pi/Prange*meanErrX_Post,2*pi/Prange*popref.meanErrX_Post(k));
        else
            meanErrX_Postdiff = NaN;
        end
        
        if k == 1
            PostErrAll_out = PostErrAll_o;
            PostCorrAll_out = PostCorrAll;
            PostErrAllpred_out = PostErrAllpred(:);
            PostCorrAllpred_out = PostCorrAllpred(:);
            meanErrX_Post_out = meanErrX_Post;
            meanCorrX_Post_out = meanCorrX_Post;
            PostErrAllpreddiff_out = PostErrAllpreddiff;
            meanErrX_Postdiff_out = meanErrX_Postdiff;
        else
            PostErrAll_out = cat(3,PostErrAll_out,PostErrAll_o);
            PostCorrAll_out = cat(3,PostCorrAll_out,PostCorrAll);
            PostErrAllpred_out = cat(2,PostErrAllpred_out,PostErrAllpred(:));
            PostCorrAllpred_out = cat(2,PostCorrAllpred_out,PostCorrAllpred(:));
            meanErrX_Post_out = cat(1,meanErrX_Post_out,meanErrX_Post);
            meanCorrX_Post_out = cat(1,meanCorrX_Post_out,meanCorrX_Post);
            PostErrAllpreddiff_out = cat(2,PostErrAllpreddiff_out,PostErrAllpreddiff(:));
            meanErrX_Postdiff_out = cat(1,meanErrX_Postdiff_out,meanErrX_Postdiff);
        end
    end
end

function [PostXAll_out, PostErrAll_out, PostBiasAll_out,...
          PostXpred_out, PostErrpred_out, PostBiaspred_out,...
          meanErrX_Post_out, meanBiasX_Post_out,...
          PostXpred_std_out, PostErrpred_std_out, PostBiaspred_std_out,...
          meanErrX_Post_std_out, meanBiasX_Post_std_out,...
          PostXpreddiff_out, PostErrpreddiff_out, meanErrX_Postdiff_out,...
          PostXpreddiff_std_out, PostErrpreddiff_std_out, meanErrX_Postdiff_std_out,...
          PostCorrAll_out, PostCorrpred_out, meanCorrX_Post_out, PostCorrpred_std_out, meanCorrX_Post_std_out] = ComputePostGrandAve(s_PostXAll0,s_ErrXAll0,s_BiasXAll0,s_CorrXAll0,popref,params,Fcenter,Fnormalize,Xsmthbin,Ysmthbin,XFcircular,YFcircular)
    if nargin < 7
        Fcenter = false;
    end
    if nargin < 8
        Fnormalize = false;
    end
    if nargin < 9
        Xsmthbin = params.Xsmthwin;
        Ysmthbin = params.Xsmthwin;%NaN;
        XFcircular = true;
        YFcircular = true;
    end
    if nargin < 11
        XFcircular = true;
        YFcircular = true;
    end
    nanimal = size(s_PostXAll0,1);
    nseries = size(s_PostXAll0,2);
    ndim4 = 0;
    for ianimal = 1:nanimal
        for iseries = 1:nseries
            if ~isempty(s_PostXAll0{ianimal,iseries}) && nansum(s_PostXAll0{ianimal,iseries}(:)) > 0
                ndim4 = max(ndim4,size(s_PostXAll0{ianimal,iseries},3));
            end
        end
    end
    for k = 1:ndim4    
        PostXAll_all = [];
        PostErrAll_all = [];
        PostCorrAll_all = [];
        for ianimal = 1:nanimal
            for iseries = 1:nseries
                if ~isempty(s_PostXAll0{ianimal,iseries}) && nansum(s_PostXAll0{ianimal,iseries}(:)) > 0 && numel(s_PostXAll0{ianimal,iseries}) > 1
                    s_PostXAll = s_PostXAll0{ianimal,iseries}(:,:,k);
                    s_ErrXAll = s_ErrXAll0{ianimal,iseries}(:,:,k);
                    s_CorrXAll = s_CorrXAll0{ianimal,iseries}(:,:,k);
%                     s_BiasXAll = s_BiasXAll0{ianimal,iseries}(:,:,k);
                    if isempty(PostXAll_all)
                        PostXAll_all = s_PostXAll;
                        PostErrAll_all = s_ErrXAll;
                        PostCorrAll_all = s_CorrXAll;
%                         PostBiasAll_all = s_BiasXAll;
                    else
                        PostXAll_all = cat(3,PostXAll_all, s_PostXAll);
                        PostErrAll_all = cat(3,PostErrAll_all, s_ErrXAll);
                        PostCorrAll_all = cat(3,PostCorrAll_all,s_CorrXAll);
%                         PostBiasAll_all = cat(3,PostBiasAll_all,s_BiasXAll);
                    end
                end
            end
        end

        PostXAll = special_smooth2D(nanmean(PostXAll_all,3), [Ysmthbin/size(PostXAll_all,1) Xsmthbin/size(PostXAll_all,2)],[YFcircular XFcircular]);
        PostErrAll = special_smooth2D(nanmean(PostErrAll_all,3), [Ysmthbin/size(PostErrAll_all,1) Xsmthbin/size(PostErrAll_all,2)],[YFcircular XFcircular]);
%         PostBiasAll = special_smooth2D(nanmean(PostBiasAll_all,3), [Ysmthbin/size(PostBiasAll_all,1) Xsmthbin/size(PostBiasAll_all,2)],[YFcircular XFcircular]);
        
%         PostCorrAll = nanmean(PostCorrAll_all,3);
        if ~isempty(popref)
            PostCorrAll = zeros(size(PostErrAll));
            errmap = PostErrAll;
            errmapref = popref.ErrXAll(:,:,k);
            for i=1:size(errmap,2)
                errmap(:,i) = errmap(:,i)-mean(errmap(:,i));
                errmap(:,i) = errmap(:,i)/(mean(errmap(:,i).^2))^0.5;
                errmapref(:,i) = errmapref(:,i)-mean(errmapref(:,i));
                errmapref(:,i) = errmapref(:,i)/(mean(errmapref(:,i).^2))^0.5;
            end
            ishift = 0;
            for xshift = -floor(size(errmap,2)/2)+1:floor(size(errmap,2)/2)
                ishift = ishift+1;
                PostCorrAll(ishift,:) = diag(errmap'*circshift(errmapref,xshift,1)/size(errmap,1));
            end
        else
            PostCorrAll = zeros(size(PostErrAll));
            errmap = PostErrAll;
            errmapref = nanmean(PostErrAll,2);
            errmapref = errmapref-nanmean(errmapref);
            errmapref = errmapref./nanmean(errmapref.^2)^0.5;
            for i=1:size(errmap,2)
                errmap(:,i) = errmap(:,i)-nanmean(errmap(:,i));
                errmap(:,i) = errmap(:,i)/(nanmean(errmap(:,i).^2))^0.5;
            end
            ishift = 0;
            for xshift = -floor(size(errmap,1)/2)+1:floor(size(errmap,1)/2)
                ishift = ishift+1;
                PostCorrAll(ishift,:) = errmap'*circshift(errmapref,xshift,1)/size(errmap,1);
            end            
        end
        
        Prange = size(PostXAll_all,1);

%         PostErrAll = zeros(size(PostXAll));
%         for i = 1:size(PostXAll,2)
%             PostErrAll(:,i) = circshift(PostXAll(:,i,:),-i + floor(size(PostXAll,2)/2) + 1);
%         end
        
        PostXpred = params.dx*(getCircularAverage(PostXAll,params.amp_thpred,params.maxtolpred,params.interpolpred));
        if Fnormalize
            PostXpred = PostXpred - nanmean(PostXpred);
            PostXpred = PostXpred./std(PostXpred);
        end
        if ~isempty(popref)
            PostXpreddiff = Prange/(2*pi)*circ_dist(2*pi/Prange*PostXpred,2*pi/Prange*popref.PostXpred(:,k));
        else
            PostXpreddiff = NaN;
        end
        PostErrpred = params.dx*(getCircularAverage(PostErrAll,params.amp_thpred,params.maxtolpred,params.interpolpred));
        if Fnormalize
            PostErrpred = PostErrpred - nanmean(PostErrpred);
            PostErrpred = PostErrpred./std(PostErrpred);
        end
        if ~isempty(popref)
            PostErrpreddiff = Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpred,2*pi/Prange*popref.ErrXpred(:,k));
        else
            PostErrpreddiff = NaN;
        end
        meanErrX_Post = params.dx*(getCircularAverage(nanmean(PostErrAll,2),params.amp_thpred,params.maxtolpred,params.interpolpred));
        if ~isempty(popref)
            meanErrX_Postdiff = Prange/(2*pi)*circ_dist(2*pi/Prange*meanErrX_Post,2*pi/Prange*popref.meanErrX_Post(k));
        else
            meanErrX_Postdiff = NaN;
        end
        
%         PostBiaspred = params.dx*(getCircularAverage(PostBiasAll,params.amp_thpred,params.maxtolpred,params.interpolpred));
%         if Fnormalize
%             PostBiaspred = PostBiaspred - nanmean(PostBiaspred);
%             PostBiaspred = PostBiaspred./std(PostBiaspred);
%         end
%         meanBiasX_Post = params.dx*(getCircularAverage(nanmean(PostBiasAll,2),params.amp_thpred,params.maxtolpred,params.interpolpred));
        
        PostCorrAlltemp = PostCorrAll;
%         PostCorrAlltemp([1:floor(size(PostCorrAll,1)/4) floor(size(PostCorrAll,1)*3/4):end],:) = 0;
        PostCorrpred = params.dx*(getCircularAverage(PostCorrAlltemp,params.amp_thcorr,params.maxtolcorr,params.interpolcorr));
        if Fnormalize
            PostCorrpred = PostCorrpred - nanmean(PostCorrpred);
            PostCorrpred = PostCorrpred./std(PostCorrpred);
        end
        meanCorrX_Post = params.dx*(getCircularAverage(nanmean(PostCorrAlltemp,2),params.amp_thcorr,params.maxtolcorr,params.interpolcorr));

        cutoff = 0;
        maskmat = cutoff*ones(size(PostXAll));
        for i = 1:size(PostXAll,2)
            [~,igood] = findfield(PostXAll(:,i),params.amp_thpred);
            maskmat(igood,i) = 1;
        end
        
        allsessions = 1:size(PostXAll_all,3);
        nfolds = numel(allsessions);
        PostXpred_std = 0;
        PostErrpred_std = 0;
%         PostBiaspred_std = 0;
        meanErrX_Post_std = 0;
%         meanBiasX_Post_std = 0;
        PostCorrpred_std = 0;
        meanCorrX_Post_std = 0;
        PostXpreddiff_std = 0;
        PostErrpreddiff_std = 0;
        meanErrX_Postdiff_std = 0;
        PostXpred_all = PostXpred;
        PostErrpred_all = PostErrpred;
%         PostBiaspred_all = PostBiaspred;
        PostCorrpred_all = PostCorrpred;
        PostXpreddiff_all = PostXpreddiff;
        PostErrpreddiff_all = PostErrpreddiff;
        meanErrX_Postdiff_all = meanErrX_Postdiff;
        if ~isempty(popref)
            PostXpred_allref = popref.PostXpred(:,k);
            PostErrpred_allref = popref.ErrXpred(:,k);
            meanErrX_Post_allref = popref.meanErrX_Post(k);
        else
            PostXpred_allref = NaN;
            PostErrpred_allref = NaN;
            meanErrX_Post_allref = NaN;
        end
        if Fcenter
            PostXpred_all = PostXpred_all - nanmean(PostXpred_all);
            PostErrpred_all = PostErrpred_all - nanmean(PostErrpred_all);
%             PostBiaspred_all = PostBiaspred_all - nanmean(PostBiaspred_all);
            PostCorrpred_all = PostCorrpred_all - nanmean(PostCorrpred_all);
            PostXpred_allref = PostXpred_allref - nanmean(PostXpred_allref);
            PostErrpred_allref = PostErrpred_allref - nanmean(PostErrpred_allref);
            PostXpreddiff_all = PostXpreddiff_all - mean(PostXpreddiff_all);
            PostErrpreddiff_all = PostErrpreddiff_all - mean(PostErrpreddiff_all);
            meanErrX_Postdiff_all = meanErrX_Postdiff_all - mean(meanErrX_Postdiff_all);
        end
        
        for kiter = 1:nfolds
            PostXAll_iter = maskmat.*special_smooth2D(nanmean(PostXAll_all(:,:,~ismember(allsessions,kiter)),3), [Ysmthbin/size(PostXAll_all,1) Xsmthbin/size(PostXAll_all,2)],[YFcircular XFcircular]);
            PostErrAll_iter = special_smooth2D(nanmean(PostErrAll_all(:,:,~ismember(allsessions,kiter)),3), [Ysmthbin/size(PostErrAll_all,1) Xsmthbin/size(PostErrAll_all,2)],[YFcircular XFcircular]);
%             PostBiasAll_iter = special_smooth2D(nanmean(PostBiasAll_all(:,:,~ismember(allsessions,kiter)),3), [Ysmthbin/size(PostBiasAll_all,1) Xsmthbin/size(PostBiasAll_all,2)],[YFcircular XFcircular]);
           
%             PostCorrAll_iter = special_smooth2D(nanmean(PostCorrAll_all(:,:,~ismember(allsessions,kiter)),3), [Ysmthbin/size(PostCorrAll_all,1) Xsmthbin/size(PostCorrAll_all,2)],[YFcircular XFcircular]);
            if ~isempty(popref)
                PostCorrAll_iter = zeros(size(PostErrAll_iter));
                errmap = PostErrAll_iter;
                errmapref = popref.ErrXAll(:,:,k);
                for i=1:size(errmap,2)
                    errmap(:,i) = errmap(:,i)-mean(errmap(:,i));
                    errmap(:,i) = errmap(:,i)/(mean(errmap(:,i).^2))^0.5;
                    errmapref(:,i) = errmapref(:,i)-mean(errmapref(:,i));
                    errmapref(:,i) = errmapref(:,i)/(mean(errmapref(:,i).^2))^0.5;
                end
                ishift = 0;
                for xshift = -floor(size(errmap,2)/2)+1:floor(size(errmap,2)/2)
                    ishift = ishift+1;
                    PostCorrAll_iter(ishift,:) = diag(errmap'*circshift(errmapref,xshift,1)/size(errmap,1));
                end
            else
                PostCorrAll_iter = zeros(size(PostErrAll_iter));
                errmap = PostErrAll_iter;
                errmapref = nanmean(PostErrAll_iter,2);
                errmapref = errmapref-nanmean(errmapref);
                errmapref = errmapref./nanmean(errmapref.^2)^0.5;
                for i=1:size(errmap,2)
                    errmap(:,i) = errmap(:,i)-nanmean(errmap(:,i));
                    errmap(:,i) = errmap(:,i)/(nanmean(errmap(:,i).^2))^0.5;
                end
                ishift = 0;
                for xshift = -floor(size(errmap,1)/2)+1:floor(size(errmap,1)/2)
                    ishift = ishift+1;
                    PostCorrAll_iter(ishift,:) = errmap'*circshift(errmapref,xshift,1)/size(errmap,1);
                end
            end
            
%             PostCorrAll_iter([1:floor(size(PostCorrAll,1)/4) floor(size(PostCorrAll,1)*3/4):end],:) = 0;
            
%             PostErrAll_iter = zeros(size(PostXAll_iter));
%             for i = 1:size(PostXAll_iter,2)
%                 PostErrAll_iter(:,i) = circshift(PostXAll_iter(:,i),-i + floor(size(PostXAll_iter,2)/2) + 1);
%             end
            PostXpred_iter = params.dx*(getCircularAverage(PostXAll_iter,0,params.maxtolpred,params.interpolpred));
            if Fcenter
                PostXpred_iter = PostXpred_iter - nanmean(PostXpred_iter);
            end
            if Fnormalize
                PostXpred_iter = PostXpred_iter - nanmean(PostXpred_iter);
                PostXpred_iter = PostXpred_iter./std(PostXpred_iter);
            end
            PostXpred_std = PostXpred_std + (nfolds - 1)/nfolds*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostXpred_iter,2*pi/Prange*PostXpred_all)).^2;
            
            if ~isempty(popref)
                PostXpreddiff_iter = Prange/(2*pi)*circ_dist(2*pi/Prange*PostXpred_iter,2*pi/Prange*PostXpred_allref);
            else
                PostXpreddiff_iter = NaN;
            end
            PostXpreddiff_std = PostXpreddiff_std + (nfolds - 1)/nfolds*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostXpreddiff_iter,2*pi/Prange*PostXpreddiff_all)).^2;
            
            PostErrpred_iter = params.dx*(getCircularAverage(PostErrAll_iter,0,params.maxtolpred,params.interpolpred));
            if Fcenter
                PostErrpred_iter = PostErrpred_iter - nanmean(PostErrpred_iter);
            end
            if Fnormalize
                PostErrpred_iter = PostErrpred_iter - nanmean(PostErrpred_iter);
                PostErrpred_iter = PostErrpred_iter./std(PostErrpred_iter);
            end
            PostErrpred_std = PostErrpred_std + (nfolds - 1)/nfolds*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpred_iter,2*pi/Prange*PostErrpred_all)).^2;
            
            if ~isempty(popref)
                PostErrpreddiff_iter = Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpred_iter,2*pi/Prange*PostErrpred_allref);
            else
                PostErrpreddiff_iter = NaN;
            end
            PostErrpreddiff_std = PostErrpreddiff_std + (nfolds - 1)/nfolds*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostErrpreddiff_iter,2*pi/Prange*PostErrpreddiff_all)).^2;
            
            meanErrX_Post_iter = params.dx*(getCircularAverage(nanmean(PostErrAll_iter,2),0,params.maxtolpred,params.interpolpred));
            meanErrX_Post_std = meanErrX_Post_std + (nfolds - 1)/nfolds*(meanErrX_Post_iter - meanErrX_Post).^2;
            
            if ~isempty(popref)
                meanErrX_Postdiff_iter = Prange/(2*pi)*circ_dist(2*pi/Prange*meanErrX_Post_iter,2*pi/Prange*meanErrX_Post_allref);
            else
                meanErrX_Postdiff_iter = NaN;
            end
            meanErrX_Postdiff_std = meanErrX_Postdiff_std + (nfolds - 1)/nfolds*(Prange/(2*pi)*circ_dist(2*pi/Prange*meanErrX_Postdiff_iter,2*pi/Prange*meanErrX_Postdiff_all)).^2;
            
%             PostBiaspred_iter = params.dx*(getCircularAverage(PostBiasAll_iter,0,params.maxtolpred,params.interpolpred));
%             if Fcenter
%                 PostBiaspred_iter = PostBiaspred_iter - nanmean(PostBiaspred_iter);
%             end
%             if Fnormalize
%                 PostBiaspred_iter = PostBiaspred_iter - nanmean(PostBiaspred_iter);
%                 PostBiaspred_iter = PostBiaspred_iter./std(PostBiaspred_iter);
%             end
%             PostBiaspred_std = PostBiaspred_std + (nfolds - 1)/nfolds*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostBiaspred_iter,2*pi/Prange*PostBiaspred_all)).^2;
            
%             meanBiasX_Post_iter = params.dx*(getCircularAverage(nanmean(PostBiasAll_iter,2),0,params.maxtolpred,params.interpolpred));
%             meanBiasX_Post_std = meanBiasX_Post_std + (nfolds - 1)/nfolds*(Prange/(2*pi)*circ_dist(2*pi/Prange*meanBiasX_Post_iter,2*pi/Prange*meanBiasX_Post)).^2;
            
            
            PostCorrpred_iter = params.dx*(getCircularAverage(PostCorrAll_iter,0,params.maxtolcorr,params.interpolcorr));
            if Fcenter
                PostCorrpred_iter = PostCorrpred_iter - nanmean(PostCorrpred_iter);
            end
            if Fnormalize
                PostCorrpred_iter = PostCorrpred_iter - nanmean(PostCorrpred_iter);
                PostCorrpred_iter = PostCorrpred_iter./std(PostCorrpred_iter);
            end
            PostCorrpred_std = PostCorrpred_std + (nfolds - 1)/nfolds*(Prange/(2*pi)*circ_dist(2*pi/Prange*PostCorrpred_iter,2*pi/Prange*PostCorrpred_all)).^2;
            meanCorrX_Post_iter = params.dx*(getCircularAverage(nanmean(PostCorrAll_iter,2),0,params.maxtolcorr,params.interpolcorr));
            meanCorrX_Post_std = meanCorrX_Post_std + (nfolds - 1)/nfolds*(Prange/(2*pi)*circ_dist(2*pi/Prange*meanCorrX_Post_iter,2*pi/Prange*meanCorrX_Post)).^2;
        end
        if  k == 1
            PostXAll_out = PostXAll;
            PostErrAll_out = PostErrAll;
            PostBiasAll_out = [];%PostBiasAll;
            PostCorrAll_out = PostCorrAll;
            PostXpred_out = PostXpred;
            PostErrpred_out = PostErrpred;
            PostBiaspred_out = [];%PostBiaspred;
            PostCorrpred_out = PostCorrpred;
            meanErrX_Post_out = meanErrX_Post;
            meanBiasX_Post_out = [];%meanBiasX_Post;
            meanCorrX_Post_out = meanCorrX_Post;
            PostXpreddiff_out = PostXpreddiff;
            PostErrpreddiff_out = PostErrpreddiff;
            meanErrX_Postdiff_out = meanErrX_Postdiff;
            PostXpred_std_out = sqrt(PostXpred_std);
            PostErrpred_std_out = sqrt(PostErrpred_std);
            PostBiaspred_std_out = [];%sqrt(PostBiaspred_std);
            meanErrX_Post_std_out = sqrt(meanErrX_Post_std);
            meanBiasX_Post_std_out = [];%sqrt(meanBiasX_Post_std);
            PostCorrpred_std_out = sqrt(PostCorrpred_std);
            meanCorrX_Post_std_out = sqrt(meanCorrX_Post_std);
            PostXpreddiff_std_out = sqrt(PostXpreddiff_std);
            PostErrpreddiff_std_out = sqrt(PostErrpreddiff_std);
            meanErrX_Postdiff_std_out = sqrt(meanErrX_Postdiff_std);
        else
            PostXAll_out = cat(3, PostXAll_out, PostXAll);
            PostErrAll_out = cat(3, PostErrAll_out, PostErrAll);
            PostBiasAll_out = [];%cat(3, PostBiasAll_out, PostBiasAll);
            PostCorrAll_out = cat(3, PostCorrAll_out, PostCorrAll);
            PostXpred_out = cat(2, PostXpred_out, PostXpred);
            PostErrpred_out = cat(2, PostErrpred_out, PostErrpred);
            PostBiaspred_out = [];%cat(2, PostBiaspred_out, PostBiaspred);
            PostCorrpred_out = cat(2, PostCorrpred_out, PostCorrpred);
            meanErrX_Post_out = cat(1, meanErrX_Post_out, meanErrX_Post);
            meanBiasX_Post_out = [];%cat(1, meanBiasX_Post_out, meanBiasX_Post);
            meanCorrX_Post_out = cat(1, meanCorrX_Post_out, meanCorrX_Post);
            PostXpreddiff_out = cat(2,PostXpreddiff_out,PostXpreddiff);
            PostErrpreddiff_out = cat(2,PostErrpreddiff_out,PostErrpreddiff);
            meanErrX_Postdiff_out = cat(1,meanErrX_Postdiff_out,meanErrX_Postdiff);
            PostXpred_std_out = cat(2, PostXpred_std_out, sqrt(PostXpred_std));
            PostErrpred_std_out = cat(2, PostErrpred_std_out, sqrt(PostErrpred_std));
            PostBiaspred_std_out = [];%cat(2, PostBiaspred_std_out, sqrt(PostBiaspred_std));
            meanErrX_Post_std_out = cat(1, meanErrX_Post_std_out, sqrt(meanErrX_Post_std));
            meanBiasX_Post_std_out = [];%cat(1, meanBiasX_Post_std_out, sqrt(meanBiasX_Post_std));
            PostCorrpred_std_out = cat(2, PostCorrpred_std_out, sqrt(PostCorrpred_std));
            meanCorrX_Post_std_out = cat(1, meanCorrX_Post_std_out, sqrt(meanCorrX_Post_std));
            PostXpreddiff_std_out = cat(2,PostXpreddiff_std_out,sqrt(PostXpreddiff_std));
            PostErrpreddiff_std_out = cat(2,PostErrpreddiff_std_out,sqrt(PostErrpreddiff_std));
            meanErrX_Postdiff_std_out = cat(1,meanErrX_Postdiff_std_out,sqrt(meanErrX_Postdiff_std));
        end
    end
end

function mat_out = smooth2D(mat,lambdaSmooth)
nonvalid = isnan(mat);
mat(nonvalid) = 1;%1
G = smooth1D(repmat(mat,3,3),lambdaSmooth);
H = smooth1D(G',lambdaSmooth)';
mat_out = H(size(mat,1)+1:2*size(mat,1),size(mat,2)+1:2*size(mat,2));
mat_out(nonvalid) = NaN;
end

function mat_oo = special_smooth2D(mat_i,win,Fcircular)
% lambdaSmooth = 2;
% nonvalid = isnan(mat_i);
% mat_i(nonvalid) = 1;%1
% G = smooth1D(repmat(mat_i,3,3),lambdaSmooth);
% H = smooth1D(G',lambdaSmooth)';
% mat_o = H(size(mat_i,1)+1:2*size(mat_i,1),size(mat_i,2)+1:2*size(mat_i,2));
% mat_o(nonvalid) = NaN;

mat_oo = mat_i;
for k = 1:size(mat_i,3)
    mat_o = mat_i(:,:,k);
    mat_o(isnan(mat_i(:,:,k))) = 0;
    if Fcircular(1)
        mat_o = repmat(mat_o,[3 1]);
    end
    if Fcircular(2)
        mat_o = repmat(mat_o,[1 3]);
    end
    if sum(isnan(win)) > 0
        if isnan(win(1)) && ~isnan(win(2))
            for i = 1:size(mat_o,1)
                mat_o(i,:) = special_smooth_1d(mat_o(i,:), win(2), [], size(mat_i,2));
            end
        end
        if isnan(win(2)) && ~isnan(win(1))
            for j = 1:size(mat_o,2)
                mat_o(:,j) = special_smooth_1d(mat_o(:,j), win(1), [], size(mat_i,1));
            end
        end
    else
        mat_o = special_smooth_2d(mat_o, win, [], [], [size(mat_i,1) size(mat_i,2)]);
    end
    if Fcircular(1)
        mat_o = mat_o(floor(size(mat_o,1)/3)+1:2*floor(size(mat_o,1)/3),:);
    end
    if Fcircular(2)
        mat_o = mat_o(:,floor(size(mat_o,2)/3)+1:2*floor(size(mat_o,2)/3));
    end
    
    mat_o(isnan(mat_i(:,:,k))) = NaN;
    mat_oo(:,:,k) = mat_o;
end
end

function [predave_out,x] = Xpredcorrection(predave,predave_ref,Xrange,dx)
predave_ref(predave_ref - (1:Xrange)' > floor(Xrange/2)) = predave_ref(predave_ref - (1:Xrange)' > floor(Xrange/2)) - Xrange;
predave_ref(predave_ref - (1:Xrange)' < -floor(Xrange/2)) = predave_ref(predave_ref - (1:Xrange)' < -floor(Xrange/2)) + Xrange;
predave(predave - (1:Xrange)' > floor(Xrange/2)) = predave(predave - (1:Xrange)' > floor(Xrange/2)) - Xrange;
predave(predave - (1:Xrange)' < -floor(Xrange/2)) = predave(predave - (1:Xrange)' < -floor(Xrange/2)) + Xrange;

predave_interp = interp1(linspace(0,numel(predave),numel(predave)+1), [predave(1);predave], 0:dx:(Xrange-dx));predave_interp(isnan(predave_interp)) = 0;
predaveref_interp = interp1(linspace(0,numel(predave_ref),numel(predave_ref)+1), [predave_ref(1);predave_ref], 0:dx:(Xrange-dx));predaveref_interp(isnan(predaveref_interp)) = 0;
predave_out = zeros(size(predave_interp));
x = 0:dx:(Xrange-dx);
predaveref_interp = [predaveref_interp(1:end-1)-Xrange predaveref_interp predaveref_interp(2:end)+Xrange];
xrep = [x(1:end-1)-Xrange x x(2:end)+Xrange];
for i = 1:numel(x)
    try
    idxmatch = min(find(abs(predaveref_interp-predave_interp(i)) <= min(abs(predaveref_interp-predave_interp(i)))));
    idxmatch = idxmatch(abs(idxmatch-(numel(x)-1) - i) == min(abs(idxmatch-(numel(x)-1) - i)));
    predave_out(i) = xrep(round(idxmatch)) - x(i);
    catch
        keyboard
    end
end
end

function [PostErrPhsAll_aligned,PostErrpredPhsAll_aligned,PostCorrPhsAll_aligned,PostCorrpredPhsAll_aligned,PostErrpredPhsShift] = PhaseAlignement(PostErrpredPhsAllRef,PostErrpredPhsAll,PostErrPhsAll,PostCorrPhsAll,PostCorrpredPhsAll)
nphsbins = numel(PostErrpredPhsAll);
nphsbins_interp = 360/5;
phsbins = 360*(0:nphsbins)/nphsbins;
phsbins_interp = 360*(0:(nphsbins_interp-1))/nphsbins_interp;

if numel(PostErrpredPhsAllRef) > 1
    PostErrpredPhsRef = zeros(size(PostErrpredPhsAllRef,1)+1,1);
    PostErrpredPhsRef(1:end-1) = PostErrpredPhsAllRef;
    PostErrpredPhsRef(end) = PostErrpredPhsRef(1);
    PostErrpredPhsRef = interp1(phsbins,PostErrpredPhsRef,phsbins_interp,'spline');
else
    PostErrpredPhsRef = PostErrpredPhsAllRef;
end

PostErrpredPhs = zeros(size(PostErrpredPhsAll,1)+1,1);
PostErrpredPhs(1:end-1) = PostErrpredPhsAll;
PostErrpredPhs(end) = PostErrpredPhs(1);
PostErrpredPhs = interp1(phsbins,PostErrpredPhs,phsbins_interp,'spline');

if ~isempty(PostCorrpredPhsAll)
    PostCorrpredPhs = zeros(size(PostCorrpredPhsAll,1)+1,1);
    PostCorrpredPhs(1:end-1) = PostCorrpredPhsAll;
    PostCorrpredPhs(end) = PostCorrpredPhs(1);
    PostCorrpredPhs = interp1(phsbins,PostCorrpredPhs,phsbins_interp,'spline');
end

if ~isempty(PostErrPhsAll)
    PostErrPhs = zeros(size(PostErrPhsAll,1),size(PostErrPhsAll,2)+1);
    PostErrPhs(:,1:end-1) = PostErrPhsAll;
    PostErrPhs(:,end) = PostErrPhs(:,1);
    PostErrPhs = (interp1(phsbins,PostErrPhs',phsbins_interp,'spline'))';
end

if ~isempty(PostCorrPhsAll)
    PostCorrPhs = zeros(size(PostCorrPhsAll,1),size(PostCorrPhsAll,2)+1);
    PostCorrPhs(:,1:end-1) = PostCorrPhsAll;
    PostCorrPhs(:,end) = PostCorrPhs(:,1);
    PostCorrPhs = (interp1(phsbins,PostCorrPhs',phsbins_interp,'spline'))';
end

if numel(PostErrpredPhsRef) > 1
    PostErrpredPhscorr = xcorr(PostErrpredPhs-mean(PostErrpredPhs),PostErrpredPhsRef-mean(PostErrpredPhsRef),numel(PostErrpredPhs)/2,'coeff');
    [~,imaxcorr] = max(PostErrpredPhscorr);
    ishift = (imaxcorr-(numel(PostErrpredPhs)/2+1));
    PostErrpredPhsShift = sign(ishift)*phsbins_interp(abs(ishift)+1);
else
    ishift = sign(PostErrpredPhsRef)*(find(phsbins_interp == abs(PostErrpredPhsRef))-1);
    PostErrpredPhsShift = 0;
end

PostErrpredPhs = circshift(PostErrpredPhs,-ishift);
PostErrpredPhsAll_aligned = interp1(phsbins_interp,PostErrpredPhs,phsbins(1:end-1));

if ~isempty(PostCorrpredPhsAll)
    PostCorrpredPhs = circshift(PostCorrpredPhs,-ishift);
    PostCorrpredPhsAll_aligned = interp1(phsbins_interp,PostCorrpredPhs,phsbins(1:end-1));
else
    PostCorrpredPhsAll_aligned = [];
end

if ~isempty(PostErrPhsAll)
    PostErrPhs = circshift(PostErrPhs,-ishift,2);
    PostErrPhsAll_aligned = (interp1(phsbins_interp,PostErrPhs',phsbins(1:end-1)))';
else
    PostErrPhsAll_aligned = [];
end

if ~isempty(PostCorrPhsAll)
    PostCorrPhs = circshift(PostCorrPhs,-ishift,2);
    PostCorrPhsAll_aligned = (interp1(phsbins_interp,PostCorrPhs',phsbins(1:end-1)))';
else
    PostCorrPhsAll_aligned = [];
end
end
                                